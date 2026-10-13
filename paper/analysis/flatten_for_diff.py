"""Flatten \\input{tables/...} and expand the number macros (numbers_*.tex) so latexdiff can mark
changes made inside tables and generated text too. Usage: python3 -I flatten_expand.py <paper_dir> <out.tex>"""
import re, sys
from pathlib import Path

root, out = Path(sys.argv[1]), Path(sys.argv[2])
s = (root / 'main.tex').read_text(encoding='utf8')


def braced(text, i):
    """text[i] == '{': return (content, index after the closing brace)."""
    depth, j = 0, i
    while True:
        c = text[j]
        if c == '\\':
            j += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return text[i + 1:j], j + 1
        j += 1


macros = {}
for f in sorted(root.glob('numbers_*.tex')):
    t = f.read_text(encoding='utf8')
    for m in re.finditer(r'\\newcommand\{\\([A-Za-z]+)\}', t):
        body, _ = braced(t, m.end())
        macros[m.group(1)] = body

# drop the \input of number files; inline every other \input (tables)
s = re.sub(r'\\input\{numbers_[a-z]+\.tex\}[^\n]*\n', '', s)
s = re.sub(r'\\input\{(tables/[^}]+)\}', lambda m: (root / (m.group(1) if m.group(1).endswith('.tex') else m.group(1) + '.tex')).read_text(encoding='utf8'), s)

names = sorted(macros, key=len, reverse=True)
pat = re.compile(r'\\(' + '|'.join(names) + r')(?![A-Za-z])(\{\})?')
for _ in range(3):                       # definitions may use other macros
    s = pat.sub(lambda m: macros[m.group(1)], s)   # a function's return value is inserted literally
out.write_text(s, encoding='utf8')
print(f'{len(macros)} macros expanded; left:', sorted(set(re.findall(pat, s))))
