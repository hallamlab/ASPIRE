#!/usr/bin/env python3
"""Check the shared MP nodal workflow contract without third-party dependencies."""
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
NS = '{http://www.w3.org/2000/svg}'

def check(path):
    root = ET.parse(path).getroot()
    source = path.read_text()
    errors = []
    if path.name == 'workflow-brief.svg':
        approved = ROOT / 'docs/diagrams/workflow-brief.source.svg'
        if path.read_bytes() != approved.read_bytes():
            errors.append('brief workflow must exactly reproduce the approved Illustrator SVG')
        return errors
    if 'Trebuchet MS' not in source:
        errors.append('ASPIRE workflow text must use the approved Trebuchet MS substitute typography')
    if any(font in source for font in ('Times New Roman', 'Arial', 'Helvetica')):
        errors.append('remove legacy workflow font declarations')
    circles = list(root.iter(NS + 'circle'))
    for color, role in [('#DAE8F6', 'input'), ('#D5E7D4', 'output')]:
        count = sum(n.get('fill', '').upper() == color for n in circles)
        if count < 2:
            errors.append(f'keep circular {role} nodes beyond the legend')
    if root.get('data-workflow-style') != 'aspire-approved-branched':
        errors.append('keep the detailed shared-chain / ASV / ASV–VOC branch structure')
    diamonds = [n for n in root.iter(NS+'path')
                if len(re.findall('[Ll]', n.get('d', ''))) >= 3
                and n.get('fill', '').upper() == '#FFFFFF']
    if len(diamonds) < 3:
        errors.append('keep diamond compute nodes beyond the legend')
    if not list(root.iter(NS+'marker')):
        errors.append('keep directional connectors')
    return errors

def main():
    paths = [ROOT / 'docs/assets' / name for name in
             ('workflow.svg', 'workflow-main.svg', 'workflow-brief.svg')]
    paths = [p for p in paths if p.exists()]
    if not paths:
        print('No primary workflow SVGs in this checkout; apply docs/WORKFLOW_STYLE.md to new figures.')
        return 0
    failed = False
    for path in paths:
        errors = check(path)
        failed |= bool(errors)
        print(('FAIL' if errors else 'PASS') + ' ' + str(path.relative_to(ROOT)))
        for error in errors:
            print('  ' + error)
    return int(failed)

if __name__ == '__main__':
    sys.exit(main())
