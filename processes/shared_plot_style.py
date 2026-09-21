#!/usr/bin/env python3
"""Repository-wide publication figure contract for ASPIRE.

Importing scripts call :func:`install_publication_style` once.  The save hook is
intentional: it catches seaborn/theme calls made later in legacy renderers and
enforces typography and panel labels immediately before every export.
"""

from __future__ import annotations

import string
import subprocess
from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.collections import QuadMesh


FONT_FAMILY = "Times New Roman"
MIN_FONT_SIZE = 22.0
PANEL_FONT_SIZE = 28.0
_INSTALLED = False
_ORIGINAL_SAVEFIG = Figure.savefig


def _register_system_font() -> list[Path]:
    """Register Times files that Fontconfig knows but Matplotlib has not cached.

    Nextflow gives each task a private ``MPLCONFIGDIR``.  A pre-existing cache in
    that directory can omit a font that is correctly installed system-wide, so
    relying on ``findfont`` alone is not sufficient.
    """
    candidates: set[Path] = set()
    try:
        result = subprocess.run(
            ["fc-list", ":family=Times New Roman", "-f", "%{file}\n"],
            check=True, capture_output=True, text=True,
        )
        candidates.update(Path(line.strip()) for line in result.stdout.splitlines() if line.strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    for root in (Path("/usr/share/fonts"), Path("/usr/local/share/fonts")):
        if root.is_dir():
            candidates.update(root.rglob("*Times*Roman*.ttf"))
            candidates.update(root.rglob("times*.ttf"))
    registered = []
    for path in sorted(candidates):
        if not path.is_file():
            continue
        try:
            font_manager.fontManager.addfont(str(path))
            registered.append(path)
        except (OSError, RuntimeError):
            continue
    return registered


def require_font() -> Path:
    """Resolve the exact required font, registering system files if necessary."""
    try:
        resolved = Path(font_manager.findfont(FONT_FAMILY, fallback_to_default=False))
    except ValueError:
        _register_system_font()
        try:
            resolved = Path(font_manager.findfont(FONT_FAMILY, fallback_to_default=False))
        except ValueError as exc:
            raise RuntimeError(
                f"Required publication font '{FONT_FAMILY}' is unavailable. "
                "Install the Microsoft core fonts and rebuild the font cache."
            ) from exc
    if not resolved.is_file():
        raise RuntimeError(f"Resolved publication font does not exist: {resolved}")
    return resolved


def _is_panel_axis(ax) -> bool:
    if not ax.get_visible() or ax.get_label() == "<colorbar>":
        return False
    try:
        return ax.get_subplotspec() is not None
    except AttributeError:
        return False


def _panel_axes(fig: Figure) -> list:
    candidates = [ax for ax in fig.axes if _is_panel_axis(ax)]
    unique, seen = [], set()
    for ax in candidates:
        bounds = tuple(round(value, 6) for value in ax.get_position().bounds)
        if bounds in seen:  # Skip twinned/overlaid axes occupying the same panel.
            continue
        seen.add(bounds)
        unique.append(ax)
    return unique


def apply_publication_contract(fig: Figure) -> None:
    """Enforce fonts, minimum readable sizes, and multi-panel lettering."""
    compact_typography = bool(
        getattr(fig, "_aspire_compact_publication_typography", False)
    )
    disable_panel_labels = bool(
        getattr(fig, "_aspire_disable_panel_labels", False)
    )
    # Saved captions/notes below the canvas are prohibited by the contract.
    for text in fig.texts:
        if text.get_position()[1] < 0.02:
            text.set_visible(False)
    for text in fig.findobj(match=mpl.text.Text):
        text.set_fontfamily(FONT_FAMILY)
        if (
            not compact_typography
            and text.get_visible()
            and text.get_text()
            and text.get_fontsize() < MIN_FONT_SIZE
        ):
            text.set_fontsize(MIN_FONT_SIZE)

    axes = _panel_axes(fig)
    # Contract heatmaps/clustermaps are grayscale; zero maps to white in Greys.
    for artist in fig.findobj(match=lambda item: isinstance(item, (AxesImage, QuadMesh))):
        if not getattr(artist, "_aspire_preserve_categorical_cmap", False):
            artist.set_cmap("Greys")
    if len(axes) > 1 and fig._suptitle is not None and all(ax.get_title().strip() for ax in axes):
        fig._suptitle.set_visible(False)
    if len(axes) > 1 and not disable_panel_labels:
        for index, ax in enumerate(axes):
            letter = string.ascii_uppercase[index] if index < 26 else f"A{index + 1}"
            existing = [item for item in ax.texts if str(item.get_text()).strip() == letter]
            if existing:
                panel_text = existing[0]
                panel_text.set_position((-0.08, 1.02))
                panel_text.set_transform(ax.transAxes)
            else:
                panel_text = ax.text(
                    -0.08, 1.02, letter, transform=ax.transAxes,
                    ha="right", va="bottom", clip_on=False,
                )
            panel_text.set_fontfamily(FONT_FAMILY)
            panel_text.set_fontsize(PANEL_FONT_SIZE)
            panel_text.set_fontweight("bold")


def _contract_savefig(self: Figure, *args, **kwargs):
    apply_publication_contract(self)
    kwargs.setdefault("dpi", 300)
    kwargs.setdefault("bbox_inches", "tight")
    result = _ORIGINAL_SAVEFIG(self, *args, **kwargs)
    if args and isinstance(args[0], (str, Path)):
        requested = Path(args[0])
        if requested.suffix.lower() in {".pdf", ".png", ".svg"}:
            for suffix in (".pdf", ".png"):
                counterpart = requested.with_suffix(suffix)
                if counterpart == requested:
                    continue
                counterpart_kwargs = dict(kwargs)
                counterpart_kwargs.pop("format", None)
                _ORIGINAL_SAVEFIG(self, counterpart, **counterpart_kwargs)
    return result


def install_publication_style() -> None:
    global _INSTALLED
    if _INSTALLED:
        return
    require_font()
    mpl.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.serif": [FONT_FAMILY],
        "font.size": MIN_FONT_SIZE,
        "axes.titlesize": 24,
        "axes.labelsize": MIN_FONT_SIZE,
        "xtick.labelsize": MIN_FONT_SIZE,
        "ytick.labelsize": MIN_FONT_SIZE,
        "legend.fontsize": MIN_FONT_SIZE,
        "legend.title_fontsize": MIN_FONT_SIZE,
        "figure.titlesize": 24,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    Figure.savefig = _contract_savefig
    _INSTALLED = True
