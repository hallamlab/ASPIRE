"""Generate stand-alone SVGs and their display tables from a completed CAMI vignette."""
import argparse
import re
import sys
from pathlib import Path
from shutil import copyfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle, Patch
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "processes" / "graph_network"))
from network_layout import spaced_network_positions
from vignette_style import apply_style, heatmap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import false_discovery_control

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
modules = args.results / "modules"
out = args.output
out.mkdir(parents=True, exist_ok=True)
apply_style()


def read(module, name):
    return pd.read_csv(modules / module / "tables" / name, sep="\t")


def save(fig, name):
    fig.savefig(out / (name + ".svg"), metadata={"Date": None}, bbox_inches="tight")
    plt.close(fig)


def export(frame, name):
    frame.to_csv(out / (name + ".tsv"), sep="\t", index=True)


meta = read("metadata_plots", "metadata_updated_micro.tsv").set_index("sample_id")
meta = meta.sort_values(["Case", "Participant_ID", "Type_Group", "lung_status"])
sites = ["Airways", "Oral"]
palette = read("metadata_plots", "Type_Group_palette.tsv").set_index("value").color
colors = [palette[s] for s in sites]
for rank in ("Family", "Phylum"):
    m = read("clustermaps", f"clustermap_{rank}_plot.tsv").set_index(rank + "_plot")
    m = m.loc[:, meta.index]
    rel = m.div(m.sum(axis=0)) * 100
    name = rank.lower() + "-patient-heatmap"
    mat = np.log10(1 + m)
    fig = plt.figure(figsize=(14, max(4.5, len(m)*0.23 + 2)))
    grid = fig.add_gridspec(3, 3, height_ratios=[1, 0.16, max(2, len(m)*0.19)],
                            width_ratios=[1.2, 8, 0.16], wspace=0.02, hspace=0.05,
                            left=0.02, right=0.70, bottom=0.1)
    top_tree = fig.add_subplot(grid[0, 1])
    side_tree = fig.add_subplot(grid[2, 0])
    order = []
    for profiles, axis, orientation in [(mat, side_tree, "left"), (mat.T, top_tree, "top")]:
        leaves = list(range(len(profiles)))
        if len(profiles) > 1:
            tree = dendrogram(linkage(profiles, method="ward", optimal_ordering=True),
                              ax=axis, orientation=orientation, no_labels=True,
                              link_color_func=lambda _: "#555555")
            leaves = tree["leaves"]
        order.append(leaves)
        axis.axis("off")
    side_tree.invert_yaxis()
    mat = mat.iloc[order[0], order[1]]
    export(mat, name)
    ordered = meta.loc[mat.columns, ["Participant_ID", "Type_Group", "Case", "lung_status"]]
    export(ordered, rank.lower() + "-heatmap-column-order")
    strip = fig.add_subplot(grid[1, 1])
    heatmap(strip, [pd.Categorical(ordered.Type_Group, categories=sites).codes],
                 cmap=ListedColormap(colors), vmin=0, vmax=1)
    strip.axis("off")
    fig.legend(handles=[Patch(color=c, label=s) for s,c in zip(sites, colors)],
               loc="upper center", ncol=2, frameon=False)
    ax = fig.add_subplot(grid[2, 1])
    img = heatmap(ax, mat, cmap="BuPu")
    ax.yaxis.tick_right()
    ax.set_yticks(range(len(mat)), mat.index, fontsize=10)
    ticks = np.arange(0, len(ordered), 20)
    ax.set_xticks(ticks, ticks + 1)
    ax.set_xlabel("Library position in the clustered matrix")
    bar = fig.add_axes([0.935, 0.31, 0.012, 0.30])
    fig.colorbar(img, cax=bar, label="log₁₀(count + 1)").solids.set_rasterized(False)
    save(fig, name)
    if rank == "Family":
        top = rel.mean(axis=1).nlargest(10).index
        patient = rel.loc[top].T.join(meta[["Participant_ID","Type_Group"]]).groupby(
            ["Participant_ID","Type_Group"]).mean()
        export(patient, "family-abundance")
        fig, ax = plt.subplots(figsize=(9,4), layout="constrained")
        for i,site in enumerate(sites):
            vals = patient.xs(site, level="Type_Group")
            ax.boxplot(vals.to_numpy(), positions=np.arange(len(top))+(i-0.5)*0.36, widths=0.32,
                       patch_artist=True, boxprops={"facecolor":colors[i]}, showfliers=False,
                       medianprops={"color":"black"})
        ax.set_xticks(range(len(top)),top,rotation=45,ha="right")
        ax.set_yscale("symlog", linthresh=0.01, linscale=0.6)
        ticks = [0, 0.01, 0.1, 1, 10, 100]
        ax.set_yticks(ticks, [str(t) for t in ticks])
        ax.set_ylim(0, 110)
        ax.set_ylabel("Patient mean relative abundance (%)")
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.legend(handles=[Patch(color=c,label=s) for s,c in zip(sites,colors)],frameon=False)
        save(fig, "family-abundance")
    else:
        comp = rel.T.join(meta.Type_Group).groupby("Type_Group").mean().loc[sites].T
        comp = comp.sort_index()
        export(comp, "phylum-composition")
        fig, ax = plt.subplots(figsize=(8,2.7),layout="constrained")
        base=np.zeros(len(sites))
        hues = plt.get_cmap("tab20").colors
        for index, taxon in enumerate(comp.index):
            color = hues[index % len(hues)]
            ax.barh(sites,comp.loc[taxon],left=base,color=color,label=taxon)
            base+=comp.loc[taxon].to_numpy()
        ax.set(xlim=(0,100),xlabel="Mean library composition (%)")
        for label in ax.get_yticklabels():
            label.set_color(palette[label.get_text()])
        ax.legend(loc="center left",bbox_to_anchor=(1,0.5),frameon=False)
        save(fig,"phylum-composition")


partitions = {}
for mode, suffix in [("Combinations allowed",""),("Single group only","_DULEG")]:
    t=read("indicator_analysis",f"Type_Group_indicator_species{suffix}_summary.tsv")
    t=t[t["significant"].astype(str).str.lower().eq("true")]
    flags=t[["s."+s for s in sites]]
    labels=flags.apply(lambda row:"+".join(s for s in sites if row["s."+s]==1),axis=1)
    partitions[mode]=labels.value_counts()
partition=pd.DataFrame(partitions).reindex(sites+["+".join(sites)]).fillna(0).astype(int)
export(partition,"indicator-partitions")
ax=partition.T.plot.bar(color=[*colors, "#aaaaaa"], rot=0, figsize=(6,3.2))
ax.set_ylabel("Significant ASVs")
ax.legend(frameon=False)
save(ax.figure,"indicator-partitions")

# Each contrast uses one library per patient per anatomical group.
counts = read("non_target_filtering", "ASV_target.tsv").set_index("ASV_ID").loc[:, meta.index]
proportions = counts.div(counts.sum(axis=0)).T
airway = meta[meta.Type_Group.eq("Airways")]
oral = meta[meta.Type_Group.eq("Oral")]
arms = {
    "Tumor-side": airway[airway.lung_status.eq("TumorSide") & airway.Case.eq("Cancer")],
    "Contralateral": airway[airway.lung_status.eq("Contralateral") & airway.Case.eq("Cancer")],
    "Healthy": airway[airway.lung_status.eq("Healthy") & airway.Case.eq("Control")],
    "Cancer": oral[oral.Case.eq("Cancer")],
    "Control": oral[oral.Case.eq("Control")],
}
contrasts = [
    ("Tumor-side", "Healthy", False, "tumor-side-vs-healthy"),
    ("Contralateral", "Healthy", False, "contralateral-vs-healthy"),
    ("Tumor-side", "Contralateral", True, "tumor-side-vs-contralateral"),
    ("Cancer", "Control", False, "oral-case"),
]
edges=pd.read_csv(modules/"network_analysis/tables/spieceasi_edge_list.csv")
nodes=pd.read_csv(modules/"network_analysis/tables/spieceasi_node_features.csv").Taxon
graph=nx.Graph()
graph.add_nodes_from(nodes)
graph.add_weighted_edges_from(edges[["Taxon1","Taxon2","Weight"]].itertuples(index=False,name=None))
pos=nx.spring_layout(graph,seed=42,weight=None)
# Each contrast uses the same positions, including room for its largest indicator.
pos = spaced_network_positions(pos, dict.fromkeys(graph, 5.2), drawing_span=180, gap=1.5)
fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
cohorts = []
for ax, (first, second, paired, slug) in zip(axes.flat, contrasts):
    left, right = [arms[name].sort_values("Participant_ID") for name in (first, second)]
    if left.empty or right.empty or any(t.Participant_ID.duplicated().any() for t in (left, right)):
        raise ValueError(f"{slug} requires one library per patient in each populated arm")
    shared = set(left.Participant_ID) & set(right.Participant_ID)
    if paired:
        if not shared:
            raise ValueError(f"No complete pairs for {slug}")
        left = left[left.Participant_ID.isin(shared)]
        right = right[right.Participant_ID.isin(shared)]
        assert left.Participant_ID.tolist() == right.Participant_ID.tolist()
    elif shared:
        raise ValueError(f"Patient overlap in independent contrast {slug}")
    for label, arm in [(first, left), (second, right)]:
        cohorts.append(arm[["Participant_ID"]].assign(contrast=slug, arm=label, paired=paired))
    frame = proportions.loc[left.index.append(right.index)]
    frame = frame.loc[:, frame.gt(0).any()]
    X = frame.to_numpy()
    present = X > 0
    case = np.arange(len(frame)) < len(left)
    n_case, n_control = case.sum(), (~case).sum()
    totals = X.sum(axis=0)
    prevalence = present.sum(axis=0)
    binary = present.astype(float)

    def indicator(masks):
        case_sum = masks @ X
        case_present = masks @ binary
        means = np.stack([case_sum / n_case, (totals - case_sum) / n_control])
        fidelity = np.stack([case_present / n_case, (prevalence - case_present) / n_control])
        return np.sqrt(np.maximum(0, means / means.sum(axis=0) * fidelity))

    observed = indicator(case.astype(float))
    stat = observed.max(axis=0)
    rng = np.random.default_rng(42)
    exceed = np.zeros(len(stat), dtype=int)
    permutations = 99_999
    for start in range(0, permutations, 1000):
        size = min(1000, permutations - start)
        if paired:
            # Rows are aligned by patient across the two halves of the matrix.
            swaps = rng.integers(0, 2, size=(size, n_case))
            masks = np.concatenate([swaps, 1 - swaps], axis=1).astype(float)
        else:
            masks = np.stack([rng.permutation(case) for _ in range(size)]).astype(float)
        exceed += (indicator(masks).max(axis=0) >= stat - 1e-12).sum(axis=0)
    p = (1 + exceed) / (1 + permutations)
    result=pd.DataFrame({"group":np.where(observed.argmax(axis=0)==0,first,second),
                         "stat":stat,"p":p,"q":false_discovery_control(p)},index=frame.columns)
    result["highlighted"]=(result.q<0.05)&(result.stat>=0.25)
    result["in_network"]=result.index.isin(graph)
    export(result, slug + "-indicators")
    ax.add_collection(LineCollection([(pos[a],pos[b]) for a,b in graph.edges],
                                    colors="#dddddd",linewidths=0.35,zorder=0))
    for group,color in [("Background","#c9c9c9"),(second,"#858585"),(first,"#181818")]:
        hits=result[result.highlighted & result.in_network]
        ids=[n for n in graph if n not in hits.index] if group=="Background" else hits.index[hits.group==group]
        if len(ids):
            xy=np.array([pos[n] for n in ids])
            ax.scatter(*xy.T,marker="^" if group==first else "o",s=8 if group=="Background" else 90*result.loc[ids,"stat"]**2,
                       color=color,edgecolors="white",linewidths=0.25)
    title = "Oral case vs control" if slug == "oral-case" else f"{first} vs {second.lower()}"
    ax.set_title(f"{title} · {len(hits)} network indicators",
                 color=palette["Oral" if slug == "oral-case" else "Airways"])
    ax.legend(handles=[Patch(color=c, label=g) for g,c in [(first,"#181818"),(second,"#858585"),
                                                          ("Background","#c9c9c9")]],
              loc="lower center", bbox_to_anchor=(0.5,-0.06), ncol=3, frameon=False, fontsize=10)
    ax.set_aspect("equal")
    ax.axis("off")
export(pd.concat(cohorts), "indicator-contrast-cohorts")
save(fig,"site-case-networks")

# Show every defined case-versus-control statistic, including nonsignificant tests.
case_fit = read("indicator_analysis", "Case_indicator_species_DULEG_summary.tsv").set_index("ASV")
tested = case_fit.dropna(subset=["stat", "q.value"])
fig, ax = plt.subplots(figsize=(6.2, 4.1), layout="constrained")
criteria = tested["q.value"] < 0.05
for label, keep, marker, shade in [
    ("Nonsignificant", ~criteria, "o", "#bcbcbc"),
    ("Case", criteria & tested["s.Cancer"].eq(1), "^", "#161616"),
    ("Control", criteria & tested["s.Control"].eq(1), "s", "#525252"),
]:
    sub = tested.loc[keep]
    ax.scatter(sub.stat, -np.log10(sub["q.value"]), marker=marker, c=shade,
               s=27, alpha=0.8, label=f"{label} ({len(sub)})")
ax.axhline(-np.log10(0.05), color="#777777", linestyle=":", linewidth=1)
ax.set(xlabel="Indicator statistic", ylabel="−log₁₀(adjusted p-value)",
       title=f"Case vs control · {len(tested)} scored, {len(case_fit)-len(tested)} unscored")
ax.legend(frameon=False)
save(fig, "case-isa-single-group")
export(case_fit, "case-isa-single-group")

def vignette_svg(source, target, keep_site=False):
    """Desaturate exported vector paint while keeping site labels and legend keys consistent."""
    preserve_hues = "phylum" in source.stem.lower()
    def neutral(match):
        paint, code = match.groups()
        if preserve_hues or (keep_site and code.lower() in {"009e73", "6a3d9a"}):
            return match.group()
        rgb = [int(code[i:i+2], 16) for i in (0, 2, 4)]
        gray = round(np.dot(rgb, [0.2126, 0.7152, 0.0722]))
        return paint + "#" + f"{gray:02x}" * 3
    svg = re.sub(r"((?:fill|stroke):\s*)#([0-9a-fA-F]{6})\b", neutral, source.read_text())
    if source.stem == "power_curves_cancer_vs_control_simple":
        for site in sites:
            svg = re.sub(r'(<text style="[^"]*?fill: )#[a-fA-F0-9]{6}("[^>]*>' + site + r'</text>)',
                         lambda m: m[1] + palette[site] + m[2], svg)
    svg = re.sub(r"font-size: ([0-9.]+)px",
                 lambda m: f"font-size: {min(12, max(10, round(float(m[1]))))}px", svg)
    svg = re.sub(r'font-family: [^;"]+', "font-family: 'DejaVu Sans', sans-serif", svg)
    target.write_text(svg)

# Collect the matching native exports without modifying the published run.
for module, plot in [

    ("contamination_filtering", "three_tier_results/read_depth_by_class"),
    ("upset", "final_micro_upset"),
    ("bubbleplotter", "bubble_plot_asv_summary_type_group"),
    ("lung_status_analysis", "Airways/Distance_contrasts_ABC"),
    ("lung_status_analysis", "Airways/PCoA_lung_status"),
    ("diversity", "patient_aware/Bray_Boxplot_sample_type_within_patient"),
    ("power_analysis", "power_curves_cancer_vs_control_simple"),
    ("indicator_analysis", "Type_Group__Case/type_group_ISA_plot"),
    ("indicator_analysis", "Type_Group__Case/type_group_ISA_plot_Phylum"),
]:
    source = modules / module / "plots" / (plot + ".svg")
    vignette_svg(source, out / source.name, keep_site=plot.split("/")[-1] in {
        "type_group_ISA_plot", "bubble_plot_asv_summary_type_group", "final_micro_upset", "read_depth_by_class"})

for filename in (
    "network_Type_Group_ISA", "network_Case_ISA", "network_betweenness_POS_ALL",
    "network_modules_POS_ALL", "network_modules_POS_SUB", "network_phylum_ISA",
    "network_phylum_ABUND_MEAN", "network_phylum_ABUND_MEDIAN",
):
    source = modules / "network_analysis/plots" / (filename + ".svg")
    vignette_svg(source, out / source.name, keep_site=filename == "network_Type_Group_ISA")
for filename in ("case_ISA_plot", "Combined_ISA_plot", "Combined_ISA_plot_Phylum"):
    source = modules / "indicator_analysis/plots/Type_Group__Case" / (filename + ".svg")
    vignette_svg(source, out / source.name, keep_site=filename == "Combined_ISA_plot")
sankey_dir = modules / "sankey/plots"
copyfile(sankey_dir / "data_loss_sankey.label.html", out / "data_loss_sankey.label.html")
all_library_sankey = sankey_dir / "data_loss_sankey_raw.label.html"
if all_library_sankey.is_file():
    copyfile(all_library_sankey, out / all_library_sankey.name)
presence = read("upset", "final_micro_final_venn_presence_table.tsv")
overlap = presence.groupby("grouping").ASV_ID.nunique().reindex(
    ["Airways", "Airways + Oral", "Oral"], fill_value=0)
export(overlap, "site-overlap")
fig, ax = plt.subplots(figsize=(5.6, 2.8))
for center, color in zip([-0.55, 0.55], colors):
    ax.add_patch(Circle((center, 0), 0.95, facecolor=color, alpha=0.27, edgecolor=color))
for x, count in zip([-0.85, 0, 0.85], overlap):
    ax.text(x, 0, str(count), ha="center", va="center", fontsize=12)
for x, group in zip([-0.95, 0.95], sites):
    ax.text(x, 1, group, ha="center")
ax.set(xlim=(-1.6, 1.6), ylim=(-1, 1.25), aspect="equal")
ax.axis("off")
save(fig, "site-overlap")
