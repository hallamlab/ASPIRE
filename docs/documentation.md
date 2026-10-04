# Maintain and publish the documentation

The README holds installation basics, the reviewer test and links into the
complete guides. Markdown pages under `docs/` are shared between GitHub and
Read the Docs; `docs/index.md` defines navigation. The reviewer guide includes
`examples/MOCK_DATASET_TESTING.md` so its detailed test procedure has one source.

## Build locally

From the repository root, create a separate documentation environment:

```bash
mamba create -n aspire-docs --override-channels -c conda-forge python=3.11 pip -y
conda activate aspire-docs
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs /tmp/aspire-docs-html
```

Open `/tmp/aspire-docs-html/index.html`. No workflow environment, references or
sample data are needed to build the documentation. Workflow previews are committed SVGs and display offline without Mermaid JavaScript. The original publication SVGs and PDFs remain available as downloads.

## Connect Read the Docs

1. Sign in to [Read the Docs](https://app.readthedocs.org/) using your GitHub account.
2. Add a project from `hallamlab/ASPIRE`. Grant the GitHub integration access to this repository; organization access may need an administrator.
3. Select the branch containing these documentation files for the first build. Use `docs/main-user-guide` for the initial review, then `main` after merging.
4. Keep `.readthedocs.yaml` at the repository root as the build configuration. It selects `docs/conf.py`, installs `docs/requirements.txt`, and treats Sphinx warnings as errors.
5. Start a build and inspect its log. Open the resulting site; check the reviewer page, workflow figures, Mermaid diagrams, parameter tables and PDF downloads.
6. In the project versions/settings, activate the versions you want to publish and choose the default version visitors should see. After merging, make sure the production branch is the source of the default documentation build.
7. Confirm that a subsequent push triggers a new build. The GitHub integration normally handles update notifications automatically; inspect the integration settings if it does not.

Read the Docs documents the [GitHub integration and automatic notifications](https://docs.readthedocs.com/platform/stable/reference/git-integration.html).
The GitHub App subscribes to the required events; it does not require you to
create an additional repository webhook manually.

The public guide is hosted at [hallamlab-aspire.readthedocs.io](https://hallamlab-aspire.readthedocs.io/). Keep README guide links pointed at the hosted pages; retain only installation basics, a runnable reviewer test, the main workflow figure, and support links in the README.

## Keep documentation aligned with the code

- Update the YAML template, configuration guide and process reference together.
- Regenerate the parameter catalogue with `python docs/generate_config_parameters.py`.
- Check registered-stage and parameter coverage with `python -m pytest tests/test_documentation_coverage.py` in a development environment containing pytest and PyYAML.
- Build with warnings treated as errors and inspect diagrams in a browser.
- Keep the SVG/PDF workflow figure synchronized when the analysis graph changes.

`docs/conf.py` explicitly lists published pages. Add new public pages to both
that list and the index navigation. Local manuscript drafts, audits and private
study configurations are not part of the public documentation build.

## Shared diagram scale

ASPIRE and MetaPathways use a shared 2,240-unit-wide white canvas for documentation previews. Smaller diagrams are centered without stretching; Mermaid diagrams use a common 1.5× scale to bring their 16-pixel labels close to the publication figures’ typography. Preview width is responsive, but relative scale stays consistent across pages. Click a diagram to open its SVG for closer inspection. Original publication SVG/PDF downloads stay tightly cropped.

Edit Mermaid sources under `docs/diagrams/`, or the original publication SVGs under `docs/assets/`. To rebuild the committed previews from the repository root, in the documentation environment:

```bash
python -m pip install -r docs/diagram-requirements.txt
python -m playwright install chromium
python scripts/render_workflow_diagrams.py
```

Chromium requires its usual Linux system libraries. Rendering downloads the pinned Mermaid bundle; ordinary Sphinx builds need neither Chromium nor network access for diagrams. `docs/diagrams/figures.json` records the source mapping and shared canvas size. If a future diagram needs a wider canvas, update both projects together. Do not hand-edit generated files in `docs/assets/diagrams/`.
