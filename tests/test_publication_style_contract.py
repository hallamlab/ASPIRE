import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROJECT = Path(__file__).parents[1]
PROCESSES = PROJECT / "processes"
SPEC = importlib.util.spec_from_file_location("shared_plot_style", PROCESSES / "shared_plot_style.py")
STYLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(STYLE)


PYTHON_PLOT_MARKERS = ("matplotlib", "pyplot", "seaborn", "savefig", "plotly")
R_PLOT_MARKERS = ("plot(", "pdf(", "cairo_pdf(", "png(", "svg(", "ggplot", "ggsave")


def test_every_plot_renderer_installs_contract():
    missing = []
    for path in sorted(PROCESSES.rglob("*")):
        if path.suffix.lower() not in {".py", ".r"} or path.name == "shared_plot_style.py":
            continue
        text = path.read_text(errors="ignore")
        markers = PYTHON_PLOT_MARKERS if path.suffix.lower() == ".py" else R_PLOT_MARKERS
        if not any(marker in text for marker in markers):
            continue
        if path.suffix.lower() == ".py" and "install_publication_style" not in text:
            missing.append(str(path.relative_to(PROJECT)))
        if path.suffix.lower() == ".r" and "publication_family" not in text:
            missing.append(str(path.relative_to(PROJECT)))
    assert not missing, f"Plot renderers missing publication contract: {missing}"


def test_no_renderer_declares_a_conflicting_font():
    conflicts = []
    for path in sorted(PROCESSES.rglob("*")):
        if path.suffix.lower() not in {".py", ".r"}:
            continue
        text = path.read_text(errors="ignore").lower()
        if "source sans" in text or "font.family'] = 'sans-serif'" in text:
            conflicts.append(str(path.relative_to(PROJECT)))
    assert not conflicts, f"Renderers declare non-contract fonts: {conflicts}"


def test_required_font_resolves_to_times_new_roman_file():
    resolved = STYLE.require_font()
    assert resolved.is_file()
    assert "times" in resolved.name.lower()


def test_save_hook_enforces_fonts_panels_notes_and_exports(tmp_path):
    STYLE.install_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for axis in axes:
        axis.set_title("Facet", fontfamily="DejaVu Sans", fontsize=8)
        axis.set_xlabel("Score", fontsize=7)
        axis.plot([0, 1], [0, 1])
    note = fig.text(0.5, -0.1, "remove this saved caption", fontsize=6)
    fig.savefig(tmp_path / "contract.svg")

    assert (tmp_path / "contract.svg").is_file()
    assert (tmp_path / "contract.png").is_file()
    assert (tmp_path / "contract.pdf").is_file()
    assert not note.get_visible()
    for axis, letter in zip(axes, "AB"):
        panel = [text for text in axis.texts if text.get_text() == letter]
        assert len(panel) == 1
        assert panel[0].get_fontfamily()[0] == STYLE.FONT_FAMILY
        assert panel[0].get_fontsize() >= STYLE.PANEL_FONT_SIZE
    for text in fig.findobj(match=matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontfamily()[0] == STYLE.FONT_FAMILY
            assert text.get_fontsize() >= STYLE.MIN_FONT_SIZE
    plt.close(fig)


def test_heatmaps_are_grayscale_and_redundant_suptitle_is_removed(tmp_path):
    STYLE.install_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for axis in axes:
        axis.set_title("Metric")
        axis.imshow(np.arange(4).reshape(2, 2), cmap="viridis")
    title = fig.suptitle("Redundant whole-figure title")
    fig.savefig(tmp_path / "heatmap.png")
    assert not title.get_visible()
    for axis in axes:
        assert axis.images[0].get_cmap().name == "Greys"
    plt.close(fig)
