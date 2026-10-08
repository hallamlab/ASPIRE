"""Build user documentation without importing analysis dependencies."""
import os

project = 'ASPIRE'
author = 'Hallam Lab and ASPIRE contributors'
copyright = '2026, ASPIRE contributors'
extensions = ['myst_parser', 'sphinxcontrib.mermaid']
source_suffix = {'.md': 'markdown'}
root_doc = 'index'
# Explicit publication boundary: local manuscript drafts and audits never enter
# the documentation build, even when they exist in a developer checkout.
include_patterns = [
    'index.md', 'installation.md', 'test.md', 'cami-mock.md', 'reviewer-test.md', 'getting-started.md',
    'primer-trimming.md', 'inputs.md', 'workflow.md', 'decontamination.md', 'analyses.md', 'outputs.md',
    'CONFIGURATION.md', 'CONFIG_PARAMETERS.md', 'PROCESS_REFERENCE.md',
    'EXPERT_GUIDE.md', 'VOC_STATISTICS.md', 'SPARK_FUNCTIONALITY.md',
    'troubleshooting.md', 'documentation.md', 'WORKFLOW_STYLE.md',
]
myst_heading_anchors = 4
myst_fence_as_directive = ['mermaid']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {'collapse_navigation': False, 'navigation_depth': 2}
html_title = 'ASPIRE: Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems'
html_baseurl = os.environ.get('READTHEDOCS_CANONICAL_URL', '')
html_static_path = ['assets']
html_css_files = ['docs.css']
mermaid_version = '11.12.1'
mermaid_init_config = {'startOnLoad': False, 'theme': 'neutral', 'flowchart': {'htmlLabels': False}}
mermaid_light_theme = 'neutral'
mermaid_dark_theme = 'neutral'
mermaid_fullscreen = True
mermaid_height = 'auto'
