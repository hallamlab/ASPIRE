#!/usr/bin/env python3
"""Prepare independent prevalence cohorts and apply their union without subtraction.

The depth cutoff is a biological-sample inclusion criterion. Negative controls
are deliberately exempt; zero-depth controls are reported but not modelled.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd


CLASSES = ("biological", "technical", "bio_control", "positive")


def read_counts(path):
    with open(path, newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header[1:]) != len(set(header[1:])):
        raise ValueError("Count table must have unique sample IDs")
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    if counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError("Count table must have unique ASV and sample IDs")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Expected finite, nonnegative integer ASV counts")
    return counts.astype("int64")


def save_counts(counts, path):
    counts.to_csv(path, sep="\t", index_label="ASV_ID")


def record(records, step, counts, **extra):
    entry = dict(step=step, samples=counts.shape[1], asvs=counts.shape[0],
                 nonzero_asvs=int((counts.sum(axis=1) > 0).sum()),
                 reads=int(counts.to_numpy().sum()), **extra)
    records.append(entry)
    print(json.dumps(entry), flush=True)


def align_metadata(counts, metadata, config):
    sample_col = config["sample_col"]
    class_col = config["class_col"]
    if sample_col not in metadata or class_col not in metadata:
        raise ValueError(f"Metadata requires {sample_col!r} and {class_col!r}")
    # Existing ASPIRE convention permits hyphen/underscore differences, but an
    # ambiguous normalization must never silently assign the wrong patient.
    normalize = lambda value: str(value).replace("-", "_")
    if metadata[sample_col].isna().any():
        raise ValueError("Missing metadata sample IDs")
    metadata = metadata.copy()
    metadata.index = metadata[sample_col].map(normalize)
    keys = pd.Index([normalize(s) for s in counts.columns])
    if metadata.index.has_duplicates or keys.has_duplicates:
        raise ValueError("Ambiguous/duplicate normalized sample IDs")
    missing = keys.difference(metadata.index)
    if len(missing):
        raise ValueError(f"Count samples missing metadata: {missing.tolist()}")
    # Preserve all metadata rows, including zero-read samples omitted by the
    # count-matrix writer, and every procedure/pairing column.
    lookup = dict(zip(keys, counts.columns))
    metadata[sample_col] = [lookup.get(k, s) for k, s in zip(metadata.index, metadata[sample_col])]
    metadata = metadata.set_index(sample_col, drop=False)
    mapping = {}
    for cls in CLASSES:
        labels = config["labels"][cls]
        if not isinstance(labels, list):
            raise ValueError(f"{cls}_labels must be a list of metadata values")
        for label in labels:
            label = str(label).strip()
            if label in mapping:
                raise ValueError(f"Overlapping sample class label: {label}")
            mapping[label] = cls
    classes = metadata[class_col].astype(str).str.strip().map(mapping)
    if classes.isna().any():
        bad = metadata.loc[classes.isna(), class_col].unique().tolist()
        raise ValueError(f"Unmapped sample classes in {class_col}: {bad}")
    metadata["decontam_class"] = classes
    return metadata


def prepare(counts_path, metadata_path, config, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    counts = read_counts(counts_path)
    metadata = align_metadata(counts, pd.read_csv(metadata_path, sep="\t", dtype=str), config)
    cutoff = config["min_biological_reads"]
    if not isinstance(cutoff, (int, float)) or cutoff < 0 or not np.isfinite(cutoff):
        raise ValueError("min_biological_reads must be finite and nonnegative")
    metadata["post_qc_reads"] = counts.sum().reindex(metadata.index, fill_value=0)
    metadata["in_count_table"] = metadata.index.isin(counts.columns)
    metadata["biological_pass"] = ((metadata.decontam_class == "biological") &
                                    metadata.in_count_table & (metadata.post_qc_reads > 0) &
                                    (metadata.post_qc_reads >= cutoff))
    metadata["qc_status"] = "retained_control"
    metadata.loc[metadata.post_qc_reads == 0, "qc_status"] = "zero_usable_reads"
    metadata.loc[metadata.decontam_class == "positive", "qc_status"] = "positive_qc_only"
    metadata.loc[metadata.decontam_class == "biological", "qc_status"] = "below_biological_depth"
    metadata.loc[metadata.biological_pass, "qc_status"] = "biological_pass"
    metadata.to_csv(outdir / "sample_qc.tsv", sep="\t", index=False)
    records = []
    record(records, "raw_post_qc", counts)
    bio_ids = metadata.index[metadata.biological_pass]
    if not len(bio_ids):
        raise ValueError("No biological samples pass the depth criterion")
    biological = counts.loc[:, bio_ids]
    save_counts(biological, outdir / "biological_qc_counts.tsv")
    record(records, "biological_depth_qc", biological)
    positive_ids = metadata.index[(metadata.decontam_class == "positive") & metadata.in_count_table]
    save_counts(counts.loc[:, positive_ids], outdir / "positive_qc_counts.tsv")
    record(records, "positive_qc_only", counts.loc[:, positive_ids])
    prevalence = pd.DataFrame(index=counts.index)
    for cls in ("biological", "technical", "bio_control"):
        ids = bio_ids if cls == "biological" else metadata.index[
            (metadata.decontam_class == cls) & metadata.in_count_table & (metadata.post_qc_reads > 0)]
        prevalence[f"n_{cls}"] = len(ids)
        prevalence[f"prevalence_{cls}"] = ((counts.loc[:, ids] > 0).mean(axis=1)
                                           if len(ids) else np.nan)
    prevalence.to_csv(outdir / "prevalence.tsv", sep="\t", index_label="ASV_ID")
    for arm, cls in (("TECH", "technical"), ("BIO", "bio_control")):
        options = config[cls]
        threshold = options["threshold"]
        if not 0 < threshold < 1:
            raise ValueError(f"{arm} threshold must be between 0 and 1")
        if not options["enabled"]:
            continue
        neg_ids = metadata.index[(metadata.decontam_class == cls) &
                                 metadata.in_count_table & (metadata.post_qc_reads > 0)]
        if not len(neg_ids):
            raise ValueError(f"{arm} enabled but no nonzero usable controls; supply controls or explicitly disable this arm")
        subset = counts.loc[:, list(bio_ids) + list(neg_ids)]
        save_counts(subset, outdir / f"{arm}_counts.tsv")
        pd.DataFrame({"Sample": subset.columns,
                      "is_negative": subset.columns.isin(neg_ids)}).to_csv(
                          outdir / f"{arm}_metadata.tsv", sep="\t", index=False)
        record(records, f"{arm}_input", subset, negatives=len(neg_ids), biological=len(bio_ids))
    pd.DataFrame(records).to_csv(outdir / "step_counts.tsv", sep="\t", index=False)
    (outdir / "settings.json").write_text(json.dumps(config, indent=2) + "\n")


def finalize(taxonomy_path, outdir):
    config = json.loads((outdir / "settings.json").read_text())
    calls = pd.read_csv(outdir / "prevalence.tsv", sep="\t", index_col=0)
    calls.index = calls.index.astype(str)
    taxonomy = pd.read_csv(taxonomy_path, sep="\t", dtype=str)
    id_col = "Feature ID" if "Feature ID" in taxonomy else taxonomy.columns[0]
    if "Taxon" not in taxonomy or taxonomy[id_col].duplicated().any():
        raise ValueError("Taxonomy requires unique ASV IDs and a Taxon column")
    calls.insert(0, "taxonomy", taxonomy.set_index(id_col).Taxon.reindex(calls.index).fillna("Unassigned"))
    for arm, cls in (("TECH", "technical"), ("BIO", "bio_control")):
        if config[cls]["enabled"]:
            scores = pd.read_csv(outdir / f"{arm}_scores.tsv", sep="\t", index_col=0)
            scores.index = scores.index.astype(str)
            if scores.index.has_duplicates or set(scores.index) != set(calls.index):
                raise ValueError(f"{arm} scores must cover every original ASV exactly once")
            calls[f"{arm}_score"] = scores["score"].reindex(calls.index)
            flag = scores["contaminant"].astype(str).str.lower()
            if not flag.isin(["true", "false"]).all():
                raise ValueError(f"Invalid {arm} contaminant calls")
            calls[f"{arm}_contaminant"] = flag.eq("true").reindex(calls.index)
            calls[f"{arm}_status"] = scores["status"].reindex(calls.index)
        else:
            calls[f"{arm}_score"] = np.nan
            calls[f"{arm}_contaminant"] = False
            calls[f"{arm}_status"] = "disabled"
    tech, bio_control = calls.TECH_contaminant, calls.BIO_contaminant
    calls["final_category"] = np.select([tech & bio_control, tech, bio_control], ["TECH+BIO", "TECH", "BIO"], default="CLEAN")
    calls["removal_reason"] = calls.final_category.where(tech | bio_control, "")
    calls.to_csv(outdir / "contamination_calls.tsv", sep="\t", index_label="ASV_ID")
    calls.loc[tech | bio_control].to_csv(outdir / "removed_asvs_with_taxonomy.tsv", sep="\t", index_label="ASV_ID")
    counts = read_counts(outdir / "biological_qc_counts.tsv")
    cleaned = counts.loc[~counts.index.isin(calls.index[tech | bio_control])]
    # ASV removal only: retained counts are byte-for-byte integer values from
    # the QC-passing biological table. Never subtract or rescale control counts.
    save_counts(cleaned, outdir / "ASV_cleaned.tsv")
    records = []
    record(records, "union_removal", cleaned, removed_asvs=int((tech | bio_control).sum()))
    pd.concat([pd.read_csv(outdir / "step_counts.tsv", sep="\t"), pd.DataFrame(records)],
              ignore_index=True).to_csv(outdir / "step_counts.tsv", sep="\t", index=False)
    summary = calls.final_category.value_counts().reindex(["CLEAN", "TECH", "BIO", "TECH+BIO"], fill_value=0)
    (outdir / "filtered").mkdir(exist_ok=True)
    (outdir / "filtered/filter_summary.txt").write_text(summary.to_string() + "\n" + json.dumps(records[0]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "finalize"])
    parser.add_argument("--counts", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.counts, args.metadata, json.loads(args.config.read_text()), args.outdir)
    else:
        finalize(args.taxonomy, args.outdir)
