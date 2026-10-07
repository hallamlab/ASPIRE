#!/usr/bin/env python3
"""Create a portable ASPIRE configuration for a supplied mock dataset."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import yaml


REQUIRED_DATASET_FILES = (
    "fastq_manifest.tsv",
    "sample_metadata.tsv",
    "chemistry.tsv",
    "references/mitochondria.fasta",
    "references/contaminants.fasta",
    "ground_truth_reference_filters.tsv",
    "ground_truth_feature_registry.tsv",
    "ground_truth_group_effects.tsv",
    "ground_truth_asv_chem.tsv",
    "ground_truth_network_modules.tsv",
    "ground_truth_extraction_controls.tsv",
)
REQUIRED_METADATA_COLUMNS = {
    "sample_id", "Participant_ID", "Case", "Type_Group", "lung_status", "batch",
    "DNA_conc", "is_negative_control", "is_positive_control",
}
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
THREAD_KEYS = {"threads", "ncores", "num_core", "num_cores", "n_core", "cpus", "conqur_num_core"}
CONFIG_TIERS = ("core", "standard", "optional")


def config_section(config: dict, name: str) -> dict:
    """Return a named section from tiered or legacy-flat configuration."""
    matches = []
    if name in config:
        matches.append(config[name])
    for tier in CONFIG_TIERS:
        tier_config = config.get(tier, {})
        if isinstance(tier_config, dict) and name in tier_config:
            matches.append(tier_config[name])
    if len(matches) != 1 or not isinstance(matches[0], dict):
        raise SystemExit(
            f"Configuration section '{name}' must be defined exactly once as a mapping"
        )
    return matches[0]


def write_if_changed(path: Path, content: str) -> None:
    """Preserve timestamps of identical inputs so Nextflow can resume their tasks."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file() or path.read_text() != content:
        path.write_text(content)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_checksums(dataset: Path) -> None:
    checksum_file = dataset / "checksums.tsv"
    if not checksum_file.is_file():
        return
    rows = checksum_file.read_text().splitlines()
    failures = []
    for index, line in enumerate(rows):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if index == 0 and fields[0].strip().lower() in {"relative_path", "path"}:
            continue
        if len(fields) < 2:
            failures.append(f"invalid row: {line}")
            continue
        path = dataset / fields[0]
        expected = fields[1].strip().lower()
        if not path.is_file() or sha256(path) != expected:
            failures.append(fields[0])
    if failures:
        raise SystemExit("Dataset checksum validation failed:\n  " + "\n  ".join(failures))


def read_manifest(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment="#")
    if "sample_id" not in frame.columns:
        frame = pd.read_csv(path, sep="\t", comment="#", header=None)
        frame.columns = ["sample_id", "fastq_r1", "fastq_r2"][: len(frame.columns)]
    if "fastq_r1" not in frame.columns:
        raise SystemExit(f"Manifest must contain sample_id and fastq_r1 columns: {path}")
    if "fastq_r2" not in frame.columns:
        frame["fastq_r2"] = ""
    frame = frame[["sample_id", "fastq_r1", "fastq_r2"]].fillna("")
    frame["sample_id"] = frame["sample_id"].astype(str)
    if frame["sample_id"].eq("").any() or frame["sample_id"].duplicated().any():
        raise SystemExit(f"Manifest contains empty or duplicate sample IDs: {path}")
    return frame


def resolve_fastq(raw_path: str, manifest_path: Path, dataset: Path, sample_id: str, read: str) -> Path:
    candidates = []
    if raw_path:
        supplied = Path(raw_path).expanduser()
        candidates.extend([supplied, manifest_path.parent / supplied, dataset / supplied])
    fastq_dir = dataset / "fastq"
    candidates.extend(
        path for path in fastq_dir.glob(f"{sample_id}*{read}*")
        if path.name.endswith(FASTQ_SUFFIXES)
    )
    existing = []
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate.is_file() and candidate not in existing:
            existing.append(candidate)
    if not existing:
        raise SystemExit(f"Could not resolve {read} FASTQ for sample '{sample_id}'")
    return existing[0]


def write_portable_manifest(dataset: Path, destination: Path) -> tuple[Path, set[str]]:
    source = dataset / "fastq_manifest.tsv"
    frame = read_manifest(source)
    rows = ["sample_id\tfastq_r1\tfastq_r2"]
    for row in frame.itertuples(index=False):
        r1 = resolve_fastq(str(row.fastq_r1), source, dataset, row.sample_id, "R1")
        r2 = resolve_fastq(str(row.fastq_r2), source, dataset, row.sample_id, "R2")
        rows.append(f"{row.sample_id}\t{r1}\t{r2}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_if_changed(destination, "\n".join(rows) + "\n")
    return destination.resolve(), set(frame["sample_id"])


CAMI_COLORS = {"Airways": "#009E73", "Oral": "#6A3D9A", "Skin": "#CC79A7", "Control": "#8C8C8C"}


def is_cami_metadata(metadata):
    return {"source_sample", "body_site", "study_role"} <= set(metadata.columns) and {
        "Airways", "Oral", "Skin"
    } <= set(metadata["body_site"].dropna())


def has_synthetic_clinical_design(metadata):
    return bool('synthetic_clinical' in metadata and metadata['synthetic_clinical'].fillna(False).astype(str).str.lower().isin(['true','1']).any())


def configure_cami_profile(config, project_dir, clinical=False):
    """Configure the airway/oral comparison and explicit control classes."""
    def adapt(value):
        if isinstance(value, dict):
            return {("batch" if key == "Case" and not clinical else key): adapt(item) for key,item in value.items()}
        if isinstance(value, list):
            return [adapt(item) for item in value]
        if isinstance(value, str):
            if value == "Case" and not clinical: return "batch"
            if value == "Participant_ID" and not clinical: return "source_sample"
            return value.replace("Bronchial Brush", "Airways").replace("BAL", "Oral")
        return value
    config = adapt(config)
    palette = str((project_dir / 'examples/cami_palette.tsv').resolve())
    primary = 'Airways=#009E73,Oral=#6A3D9A,Airways+Oral=#B8E186,not_indicator=#D3D3D3'
    secondary = 'Case' if clinical else 'batch'
    batches = 'Control=#8C8C8C,Cancer=#D55E00' if clinical else 'plate_1=#8C8C8C,plate_2=#0072B2'
    for section in ('metadata_plots', 'sankey'):
        config_section(config,section)['palette_file'] = palette
    config_section(config,'table_filter')['min_sample_reads'] = 5000
    # CAMI classes come from metadata: skin is the biological/sample control;
    # oral and airways are biological samples; extraction blanks are technical.
    config_section(config,'control_decontam').update(
        class_col='Type_Group', biological_labels=['Airways','Oral'],
        technical_labels=['Control'], bio_control_labels=['Skin'],
        positive_labels=['Positive'], technical_enabled=True, bio_control_enabled=True)
    config_section(config,'metadata_plots').update(
        keep_types=['Airways','Oral'], group_order=['Airways','Oral'], input_table='filtered')
    config_section(config,'sankey')['vertical_order'] = ['Airways','Oral','Skin','Control']
    config_section(config,'batch_correction')['biological_covariates'] = 'Type_Group,Case' if clinical else 'Type_Group'
    for section in ('umap_clustering','diversity','clustermaps'):
        config_section(config,section).update(group1_palette=primary, group2_palette=batches)
    for section in ('indicspecies','spieceasi'):
        config_section(config,section)['group_palettes'] = {'Type_Group':primary, secondary:batches}
        config_section(config,section)['group_orders'] = {'Type_Group':['Airways','Oral'], secondary:['Control','Cancer'] if clinical else ['plate_1','plate_2']}
    # Patient-aware branches use the explicit synthetic clinical design.
    config_section(config,'diversity')['patient_aware']['enabled'] = clinical
    for section in ('taxonomy_patient_aware','lung_status_analysis','power_analysis'):
        config_section(config,section)['enabled'] = clinical
    config_section(config,'voc_correlation').update(patient_inference=clinical, sample_types='Airways,Oral',
        isa_focus_groups=['Airways'], isa_all_type_groups=['Airways','Oral'])
    return config


def validate_tabular_inputs(dataset: Path, manifest_ids: set[str]) -> None:
    metadata = pd.read_csv(dataset / "sample_metadata.tsv", sep="\t")
    required = REQUIRED_METADATA_COLUMNS
    if is_cami_metadata(metadata) and not has_synthetic_clinical_design(metadata):
        required = required - {'Participant_ID', 'Case', 'lung_status'}
    missing_columns = sorted(required - set(metadata.columns))
    if missing_columns:
        raise SystemExit(f"Metadata is missing required columns: {', '.join(missing_columns)}")
    metadata_ids = set(metadata["sample_id"].dropna().astype(str))
    if len(metadata_ids) != len(metadata) or metadata_ids != manifest_ids:
        raise SystemExit(
            f"Manifest/metadata sample mismatch: manifest={len(manifest_ids)}, metadata={len(metadata_ids)}"
        )
    chemistry = pd.read_csv(dataset / "chemistry.tsv", sep="\t")
    if "sample_id" not in chemistry.columns:
        raise SystemExit("chemistry.tsv must contain a sample_id column")
    chemistry_ids = set(chemistry["sample_id"].dropna().astype(str))
    if not chemistry_ids or not chemistry_ids <= metadata_ids:
        raise SystemExit(
            f"Chemistry sample IDs must be a non-empty subset of metadata IDs: chemistry={len(chemistry_ids)}"
        )
    numeric = chemistry.drop(columns=["sample_id"]).select_dtypes(include="number")
    if numeric.shape[1] == 0:
        raise SystemExit("chemistry.tsv must contain at least one numeric VOC column")


def absolute_project_paths(value, project_dir: Path):
    if isinstance(value, dict):
        return {key: absolute_project_paths(item, project_dir) for key, item in value.items()}
    if isinstance(value, list):
        return [absolute_project_paths(item, project_dir) for item in value]
    if isinstance(value, str) and value.startswith("processes/"):
        return str((project_dir / value).resolve())
    return value


def default_thread_count() -> int:
    available = os.cpu_count() or 1
    return max(1, int(available * 0.8))


def apply_thread_defaults(value, threads: int):
    if isinstance(value, dict):
        updated = {}
        for key, item in value.items():
            if str(key) in THREAD_KEYS and isinstance(item, int):
                updated[key] = threads
            else:
                updated[key] = apply_thread_defaults(item, threads)
        return updated
    if isinstance(value, list):
        return [apply_thread_defaults(item, threads) for item in value]
    return value


def build_config(
    template: Path,
    dataset: Path,
    output: Path,
    runtime: Path,
    project_dir: Path,
    threads: int | None = None,
):
    threads = default_thread_count() if threads is None else threads
    if threads < 1:
        raise ValueError("threads must be at least 1")
    config = yaml.safe_load(template.read_text())
    config = absolute_project_paths(config, project_dir)
    config = apply_thread_defaults(config, threads)

    metadata = str((dataset / "sample_metadata.tsv").resolve())
    config_section(config, "paths").update(
        input_dir=str((dataset / "fastq").resolve()),
        manifest=str((dataset / "fastq_manifest.tsv").resolve()),
        output_dir=str(output.resolve()),
        runtime_dir=str(runtime.resolve()),
        keep_runtime_dir=True,
    )
    config_section(config, "table_filter")["script"] = str(
        (project_dir / "processes/filter_table/filter_ASV_table.py").resolve()
    )
    config_section(config, "mito")["mito_fasta"] = str((dataset / "references/mitochondria.fasta").resolve())
    config_section(config, "mito")["contaminant_fasta"] = str(
        (dataset / "references/contaminants.fasta").resolve()
    )
    for section in ("filter_counts", "control_decontam", "sankey", "metadata_plots", "spieceasi"):
        config_section(config, section)["metadata"] = metadata
    config_section(config, "voc_correlation")["voc_table"] = str((dataset / "chemistry.tsv").resolve())
    palette = str((project_dir / "examples/metadata_palette.template.tsv").resolve())
    config_section(config, "metadata_plots")["palette_file"] = palette
    config_section(config, "sankey")["palette_file"] = palette
    metadata_file = dataset / 'sample_metadata.tsv'
    if metadata_file.is_file() and is_cami_metadata(pd.read_csv(metadata_file, sep='\t')):
        config = configure_cami_profile(config, project_dir, has_synthetic_clinical_design(pd.read_csv(metadata_file, sep='\t')))
    manifest_file = dataset / "manifest.json"
    if manifest_file.is_file() and json.loads(manifest_file.read_text()).get("name") == "aspire-quickstart":
        # Exercise every power-analysis path with a short installation grid.
        # Truth-recovery thresholds and all other analysis settings are shared.
        config_section(config, "power_analysis").update(n_simulations=10, n_perm=49)
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path, help="Supplied ASPIRE mock dataset directory")
    parser.add_argument("--output", required=True, type=Path, help="Final ASPIRE output directory")
    parser.add_argument("--runtime", type=Path, help="Runtime directory; defaults to <output>/.aspire")
    parser.add_argument("--config-out", required=True, type=Path, help="Generated YAML path")
    parser.add_argument(
        "--threads",
        type=int,
        default=default_thread_count(),
        help="Thread/core count for generated mock config; defaults to 80%% of detected CPUs",
    )
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1")

    dataset = args.dataset.expanduser().resolve()
    output = args.output.expanduser().resolve()
    runtime = args.runtime.expanduser().resolve() if args.runtime else output / ".aspire"
    project_dir = Path(__file__).resolve().parents[2]
    template = project_dir / "examples/mock.local.yml"

    missing = [str(dataset / relative) for relative in REQUIRED_DATASET_FILES if not (dataset / relative).is_file()]
    fastq_dir = dataset / "fastq"
    if not fastq_dir.is_dir() or not any(fastq_dir.glob("*.f*q*")):
        missing.append(f"{fastq_dir} (directory containing FASTQs)")
    if missing:
        raise SystemExit("Dataset is incomplete; missing:\n  " + "\n  ".join(missing))
    validate_checksums(dataset)
    if output == runtime:
        raise SystemExit("--output and --runtime must be different directories")

    portable_manifest = args.config_out.with_suffix(".manifest.tsv").resolve()
    portable_manifest, manifest_ids = write_portable_manifest(dataset, portable_manifest)
    validate_tabular_inputs(dataset, manifest_ids)
    config = build_config(template, dataset, output, runtime, project_dir, args.threads)
    metadata_frame = pd.read_csv(dataset / 'sample_metadata.tsv', sep='\t')
    if is_cami_metadata(metadata_frame):
        metadata_frame['Color'] = metadata_frame['Type_Group'].map(CAMI_COLORS)
        if metadata_frame['Color'].isna().any():
            raise SystemExit('Unrecognized CAMI Type_Group; cannot assign a consistent palette')
        metadata_frame['source_sample'] = metadata_frame['source_sample'].fillna(metadata_frame['sample_id'])
        colored_metadata = args.config_out.with_suffix('.metadata.tsv').resolve()
        colored_metadata.parent.mkdir(parents=True, exist_ok=True)
        write_if_changed(colored_metadata, metadata_frame.to_csv(sep='\t', index=False))
        for section in ('filter_counts','control_decontam','sankey','metadata_plots','spieceasi'):
            config_section(config,section)['metadata'] = str(colored_metadata)
        print('Study profile: airway–oral comparison; TECH blanks and BIO skin prevalence filtering')
    config_section(config, "paths")["manifest"] = str(portable_manifest)
    args.config_out.parent.mkdir(parents=True, exist_ok=True)
    write_if_changed(args.config_out, yaml.safe_dump(config, sort_keys=False))
    tier_counts = ", ".join(
        f"{tier} ({len(config.get(tier, {}))} sections)" for tier in CONFIG_TIERS
    )
    enabled_optional = [
        name for name, section in config.get("optional", {}).items()
        if isinstance(section, dict) and section.get("enabled") is True
    ]
    disabled_optional = [
        name for name, section in config.get("optional", {}).items()
        if isinstance(section, dict) and section.get("enabled") is False
    ]
    print(f"Wrote mock-run configuration: {args.config_out.resolve()}")
    print(f"Wrote mock-run manifest: {portable_manifest} ({len(manifest_ids)} samples)")
    print(f"Configured mock-run threads: {args.threads}")
    print(f"Configuration tiers: {tier_counts}")
    print(f"Enabled optional benchmark modules: {', '.join(enabled_optional)}")
    print(f"Disabled (input not supplied): {', '.join(disabled_optional)}")
    print(f"Run: ./run_asv_pipeline.sh {args.config_out.resolve()} --no-resume")


if __name__ == "__main__":
    main()
