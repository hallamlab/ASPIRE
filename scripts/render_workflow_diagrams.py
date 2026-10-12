#!/usr/bin/env python3
"""Rebuild shared-scale SVG previews from Mermaid and publication SVG sources.

Run from the repository root after installing docs/diagram-requirements.txt
and running `python -m playwright install chromium`. Rendering downloads the
pinned Mermaid bundle; documentation builds use the committed SVGs offline.
"""
from pathlib import Path
import json
import argparse
import xml.etree.ElementTree as ET
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
NS = '{http://www.w3.org/2000/svg}'
ET.register_namespace('', NS[1:-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mermaid-js', type=Path, help='Local Mermaid 11.12.1 bundle for offline rendering')
    args = parser.parse_args()
    config = json.loads((ROOT / 'docs/diagrams/figures.json').read_text())
    rendered = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto('about:blank')
        page.set_default_timeout(45000)
        if args.mermaid_js:
            page.add_script_tag(path=str(args.mermaid_js.resolve()))
        else:
            page.add_script_tag(url='https://cdn.jsdelivr.net/npm/mermaid@11.12.1/dist/mermaid.min.js')
        page.evaluate('mermaid.initialize({startOnLoad:false,securityLevel:"strict"})')
        for item in config['figures']:
            code = (ROOT / 'docs/diagrams' / (item['name'] + '.mmd')).read_text()
            svg = page.evaluate('''async ({code,id}) => {
                const result = await mermaid.render(id, code);
                const div = document.createElement('div');
                div.innerHTML = result.svg;
                return new XMLSerializer().serializeToString(div.firstElementChild);
            }''', {'code': code, 'id': ROOT.name.lower() + '-' + item['name']})
            rendered.append((item['name'], ET.fromstring(svg), config['mermaid_scale']))
        browser.close()
    for name in ['workflow', 'workflow-brief']:
        source = ROOT / 'docs/assets' / (name + '.svg')
        if source.exists():
            rendered.append((name, ET.parse(source).getroot(), 1))
    width = config['canvas_width']
    destination = ROOT / 'docs/assets/diagrams'
    destination.mkdir(exist_ok=True)
    # Check every diagram before replacing any previews. Keep the common width
    # synchronized between ASPIRE and MetaPathways if a future figure needs more.
    for name, element, scale in rendered:
        natural_width = float(element.get('viewBox').split()[2]) * scale
        if natural_width > width - 40:
            raise SystemExit(f'{name} needs {natural_width + 40:g} canvas units; update the shared canvas width.')
    for name, element, scale in rendered:
        _, _, w, h = map(float, element.get('viewBox').split())
        w *= scale
        h *= scale
        canvas = ET.Element(NS + 'svg', {
            'viewBox': f'0 0 {width} {h + 40:g}', 'width': str(width),
            'height': f'{h + 40:g}', 'role': 'img',
        })
        ET.SubElement(canvas, NS + 'title').text = name.replace('-', ' ').title()
        ET.SubElement(canvas, NS + 'rect', {
            'width': str(width), 'height': f'{h + 40:g}', 'fill': '#ffffff',
        })
        element.set('x', f'{(width - w) / 2:g}')
        element.set('y', '20')
        element.set('width', f'{w:g}')
        element.set('height', f'{h:g}')
        element.attrib.pop('style', None)
        canvas.append(element)
        (destination / (name + '.svg')).write_text(ET.tostring(canvas, encoding='unicode') + '\n')
        print(f'Rendered {name} on {width}-unit canvas')


if __name__ == '__main__':
    main()
