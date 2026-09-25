import os, re, glob

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'vienna_ps_capabilities.md')

out_lines = ['# ViennaPS Current Capabilities Inventory\n']
out_lines.append('\nGenerated from `include/viennaps/*` directories.\n')

sections = [
    ('Models', 'include/viennaps/models/*.hpp'),
    ('Process / Strategy', 'include/viennaps/process/*.hpp'),
    ('Fields / Multiphysics', 'include/viennaps/fields/*.hpp'),
    ('GDS / Layout', 'include/viennaps/gds/*.hpp'),
    ('Geometries', 'include/viennaps/geometries/*.hpp'),
    ('Materials', 'include/viennaps/materials/*.hpp'),
    ('Compact', 'include/viennaps/compact/*.hpp'),
]

for section, pattern in sections:
    out_lines.append(f"\n## {section}\n")
    files = sorted(glob.glob(os.path.join(REPO_ROOT, pattern)))
    if not files:
        out_lines.append('_No files._\n')
        continue
    for f in files:
        name = os.path.basename(f)
        try:
            with open(f, 'r', encoding='utf-8', errors='ignore') as fh:
                content = fh.read()
        except Exception as e:
            out_lines.append(f"- {name} (read error: {e})\n")
            continue
        classes = re.findall(r'class\s+(\w+)', content)
        brief = ''
        for line in content.splitlines():
            m = re.search(r'@brief\s+(.+)', line)
            if m:
                brief = m.group(1).strip()
                break
        # fallback: first doxygen-ish comment line
        if not brief:
            for line in content.splitlines()[:10]:
                if line.strip().startswith('///') and not line.strip().startswith('////'):
                    brief = line.strip().lstrip('/').strip()
                    if brief:
                        break
        out_lines.append(f"\n### {name}\n")
        out_lines.append(f"- Classes: {', '.join(classes[:5])}\n")
        if brief:
            out_lines.append(f"- Brief: {brief}\n")


with open(OUT_PATH, 'w', encoding='utf-8') as f:
    f.write(''.join(out_lines))
print(f'Wrote {OUT_PATH}')
