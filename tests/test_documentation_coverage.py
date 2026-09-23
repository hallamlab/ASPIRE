import re
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def config_leaves(value, prefix=""):
    if isinstance(value, dict):
        if not value:
            yield prefix, value
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            yield from config_leaves(child, path)
    else:
        yield prefix, value


def rendered(value):
    if value is None:
        text = "null"
    elif isinstance(value, bool):
        text = str(value).lower()
    elif isinstance(value, (list, dict)):
        text = json.dumps(value, separators=(",", ":"))
    else:
        text = str(value)
    return text.replace("|", "\\|").replace("\n", " ")


def test_every_template_parameter_is_catalogued():
    config = yaml.safe_load((ROOT / "asv_pipeline_nextflow.yml").read_text())
    catalogue = (ROOT / "docs" / "CONFIG_PARAMETERS.md").read_text()
    missing = [path for path, _ in config_leaves(config) if f"`{path}`" not in catalogue]
    assert not missing

    stale = [
        path for path, value in config_leaves(config)
        if not re.search(
            rf"^\| `{re.escape(path)}` \| [^|]+ \| `{re.escape(rendered(value))}` \|",
            catalogue,
            re.M,
        )
    ]
    assert not stale


def test_every_registered_stage_is_documented():
    launcher = (ROOT / "run_asv_pipeline.sh").read_text()
    block = re.search(r"PROCESS_ORDER=\(\n(.*?)\n\)", launcher, re.S)
    assert block
    stages = block.group(1).split()
    process_reference = (ROOT / "docs" / "PROCESS_REFERENCE.md").read_text()
    missing = [stage for stage in stages if f"`{stage}`" not in process_reference]
    assert not missing


def test_public_docs_do_not_reference_private_study_configs():
    docs = "\n".join(
        path.read_text()
        for path in [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
    )
    for private_name in (
        "SPARK_compatible.yml",
        "SPARK_host_filtered_full_run.yml",
        "set1-2.local.yml",
        "set1-2.reproduction.yml",
    ):
        assert private_name not in docs
