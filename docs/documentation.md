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
sample data are needed to build the documentation. For diagrams, the generated
HTML loads Mermaid JavaScript; an offline browser may require that asset to be
cached. The SVG and PDF workflow figures are standalone files.

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

After the first successful deployment, add the assigned public site URL to the
README and repository description. The site URL is intentionally not guessed.

## Keep documentation aligned with the code

- Update the YAML template, configuration guide and process reference together.
- Regenerate the parameter catalogue with `python docs/generate_config_parameters.py`.
- Check registered-stage and parameter coverage with `python -m pytest tests/test_documentation_coverage.py` in a development environment containing pytest and PyYAML.
- Build with warnings treated as errors and inspect diagrams in a browser.
- Keep the SVG/PDF workflow figure synchronized when the analysis graph changes.

`docs/conf.py` explicitly lists published pages. Add new public pages to both
that list and the index navigation. Local manuscript drafts, audits and private
study configurations are not part of the public documentation build.
