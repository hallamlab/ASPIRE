"""Typography and vector heatmap settings for the CAMI vignette figures."""
import matplotlib.pyplot as plt
import numpy as np


def apply_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "text.color": "#292929", "svg.fonttype": "none",
        "svg.hashsalt": "aspire-cami-vignette",
    })


def heatmap(ax, values, **kwargs):
    """Draw one vector cell per value, centered on integer tick positions."""
    rows, columns = np.shape(values)
    mesh = ax.pcolormesh(np.arange(columns + 1) - 0.5, np.arange(rows + 1) - 0.5,
                         values, shading="flat", edgecolors="face", linewidth=0.1, **kwargs)
    ax.invert_yaxis()
    return mesh
