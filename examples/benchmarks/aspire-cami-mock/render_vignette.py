"""Render documentation figures from a completed CAMI mock: python render_vignette.py RESULTS DEST."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

modules = Path(sys.argv[1]) / "modules"
destination = Path(sys.argv[2])
destination.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#777777", "text.color": "#252525", "svg.fonttype": "none"})


def table(module, filename):
    return pd.read_csv(modules / module / "tables" / filename, sep="\t")


def save(fig, name):
    fig.savefig(destination / (name + ".svg"), bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


groups = ["Airways", "Oral"]
fig, ax = plt.subplots(figsize=(6.5, 3.4))
for group, color, style in zip(groups, ["#303030", "#717171"], ["-", "--"]):
    t = table("collectors_curve", group + "_collector_stats.tsv")
    ax.fill_between(t.n_samples, t.lo_2p5, t.hi_97p5, color=color, alpha=0.09)
    ax.plot(t.n_samples, t.mean_asvs, style, color=color, label=group)
    ax.annotate(str(round(t.mean_asvs.iloc[-1])), (t.n_samples.iloc[-1], t.mean_asvs.iloc[-1]),
                xytext=(3, 5), textcoords="offset points")
ax.set(xlabel="Libraries sampled", ylabel="Cumulative ASVs", xlim=(0, 81))
ax.legend(frameon=False, loc="lower right")
save(fig, "collector")

coords = table("diversity", "umap_coordinates_bray.tsv")
fig, axes = plt.subplots(1, 2, figsize=(8.1, 3.6), layout="constrained")
for group, marker, color in zip(groups, ["o", "^"], ["#333333", "#898989"]):
    t = coords[coords.Type_Group == group]
    axes[1].scatter(t.UMAP1_bray, t.UMAP2_bray, label=group, marker=marker, color=color,
                    edgecolors="white", linewidth=0.4, s=23)
for i, (group, case) in enumerate([(g, c) for g in groups for c in ["Control", "Cancer"]]):
    vals = coords.loc[(coords.Type_Group == group) & (coords.Case == case), "Shannon"]
    axes[0].boxplot(vals, positions=[i], widths=0.58, patch_artist=True,
                    boxprops={"facecolor": "#dedede" if case == "Control" else "#a0a0a0"},
                    medianprops={"color": "#222222"})
axes[0].set_xticks(range(4), ["Airways\ncontrol", "Airways\ncase", "Oral\ncontrol", "Oral\ncase"], fontsize=9)
axes[0].set(ylabel="Shannon diversity", title="Within-library diversity")
axes[1].set(xlabel="UMAP 1", ylabel="UMAP 2", title="Bray–Curtis ordination")
axes[1].legend(frameon=False, fontsize=9)
save(fig, "community")

voc = table("voc_correlation", "patient_asv_voc_permutation_long.tsv")
voc = voc[(voc.normalization == "relative_abundance") & (voc.status == "tested")]
indicators = table("indicator_analysis", "Type_Group_indicator_species_summary.tsv")
comparison = table("voc_correlation", "patient_asv_voc_normalization_comparison.tsv")
airway = indicators.loc[indicators.significant & indicators["s.Airways"].eq(1)
                        & indicators["s.Oral"].eq(0), "ASV"]
acetone = comparison[comparison.voc.eq("acetone") & comparison.support.eq("both_normalizations")
                     & comparison.same_direction].sort_values(["rho", "asv"], ascending=[False, True])
candidates = acetone[acetone.asv.isin(airway) & acetone.rho.gt(0)]
pair_columns = ["asv_label", "voc", "n_patients", "n_nonzero_patients", "rho", "q_value",
                "rho_clr", "q_value_clr"]
candidates.head(1)[pair_columns].to_csv(destination / "patient-isa-acetone.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(6.1, 3.6), layout="constrained")
examples = [*acetone.head(2).asv, *acetone.tail(1).asv]
if candidates.empty:
    ax.axis("off")
    ax.text(0.5, 0.5, "No positive airway indicator–acetone pair supported\nby both normalizations.",
            ha="center", va="center", transform=ax.transAxes)
else:
    pair = candidates.iloc[0]
    abundance = table("voc_correlation", "patient_asv_relative_abundance.tsv")
    chemistry = table("voc_correlation", "patient_voc_matrix.tsv")
    points = abundance[["patient_id", pair.asv]].merge(chemistry[["patient_id", "acetone"]],
                                                       on="patient_id", validate="one_to_one").dropna()
    ax.scatter(100 * points[pair.asv], points.acetone, s=31, color="#484848", alpha=0.7)
    ax.set(xlabel=f"{pair.asv} mean relative abundance (%)", ylabel="Mean acetone (synthetic units)",
           title=f"Airway indicator and acetone across {len(points)} patients")
    examples.insert(0, pair.asv)
save(fig, "patient-isa-acetone")
if not examples:
    examples = voc.groupby("asv").rho.count().nlargest(3).index.tolist()
examples = list(dict.fromkeys(examples))
matrix = voc.pivot(index="asv", columns="voc", values="rho").loc[examples]
fig, ax = plt.subplots(figsize=(7.8, 2.7), layout="constrained")
palette = "Greys"
im = ax.imshow(matrix, cmap=palette, vmin=-1, vmax=1, aspect="auto")
ax.set_yticks(range(len(matrix)), matrix.index)
ax.set_xticks(range(len(matrix.columns)), [v.replace("_", "\n", 1) if v == "dimethyl_sulfide"
                  else v.replace("2_3_butanedione", "2,3-butanedione") for v in matrix.columns], fontsize=9)
for (i, j), value in np.ndenumerate(matrix.to_numpy()):
    ax.text(j, i, f"{value:+.2f}", ha="center", va="center", color="white" if value > 0.1 else "#202020")
fig.colorbar(im, ax=ax, ticks=[-1, 0, 1], label="Spearman ρ")
save(fig, "patient-voc")

members = table("network_analysis", "spieceasi_modules_all.tsv")
sizes = members.groupby(["module_id", "module_label"]).size().sort_index()
fig, ax = plt.subplots(figsize=(8, 2.9))
ax.bar(range(len(sizes)), sizes, color="#646464")
ax.set_xticks(range(len(sizes)), sizes.index.get_level_values("module_label"), rotation=55, fontsize=9)
ax.set(xlabel="Leiden consensus module", ylabel="ASVs")
save(fig, "network-modules")

taxa = table("taxonomy", "cancer_vs_control/taxonomic_abundance_observed.tsv")
taxa = taxa[(taxa.tax_level == "Phylum") & (taxa.sample_type == "Airways")].copy()
taxa["peak"] = taxa[["median_control", "median_cancer"]].max(axis=1)
taxa = taxa.nlargest(5, "peak").sort_values("peak")
fig, ax = plt.subplots(figsize=(6.7, 3.2), layout="constrained")
positions = np.arange(len(taxa))
for column, offset, marker, color, label in [
    ("median_control", 0.11, "o", "#777777", "Control"),
    ("median_cancer", -0.11, "^", "#242424", "Case"),
]:
    ax.scatter(taxa[column] * 100, positions + offset, marker=marker, color=color, label=label)
ax.set_yticks(positions, taxa.taxon)
ax.set(xlabel="Median relative abundance (%)")
ax.set_xlim(left=-1.5)
ax.legend(frameon=False, loc="lower right")
save(fig, "airway-phyla")
