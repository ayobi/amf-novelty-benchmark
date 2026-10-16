#!/usr/bin/env python3
"""Word copy of the manuscript: flatten main.tex (tables inlined, number macros expanded, cross-references and
citations resolved from main.aux, figures as PNG) and convert it with pandoc.

Run after a full LaTeX build (main.aux must be current):
  python3 analysis/make_docx.py --paper . --out amf_stage3_paper_v12.docx
Needs pandoc (3.x) and pdftoppm.
"""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path


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


def sub_command(s, name, fn):
    """Replace every \\name{arg} (balanced braces) by fn(arg)."""
    out, i = [], 0
    pat = re.compile(r'\\' + name + r'(?![A-Za-z])\s*\{')
    while True:
        m = pat.search(s, i)
        if not m:
            out.append(s[i:])
            return ''.join(out)
        arg, j = braced(s, m.end() - 1)
        out.append(s[i:m.start()])
        out.append(fn(arg))
        i = j


def fit_tables(path, total_in=6.3):
    """Column widths in proportion to their longest cell text (capped), so no column collapses."""
    import docx
    from docx.shared import Inches
    d = docx.Document(path)
    for t in d.tables:
        ncol = max(len(r.cells) for r in t.rows)
        text = [0] * ncol                       # longest cell text (capped) per column
        word = [0] * ncol                       # longest single word per column
        for r in t.rows:
            for k, c in enumerate(r.cells):
                text[k] = max(text[k], min(len(c.text), 60))
                word[k] = max(word[k], max((len(w) for w in c.text.split()), default=0))
        base = [0.085 * w + 0.2 for w in word]  # inches: the longest word must fit on one line
        if sum(base) > total_in:
            base = [b * total_in / sum(base) for b in base]
        extra = total_in - sum(base)
        need = [max(0.0, 0.085 * x + 0.2 - b) for x, b in zip(text, base)]
        tot = sum(need) or 1.0
        widths = [b + extra * n / tot for b, n in zip(base, need)]
        t.autofit = False
        for r in t.rows:
            for k, c in enumerate(r.cells):
                c.width = Inches(widths[k])
        # the grid and the table width too (LibreOffice and Google Docs read these)
        from docx.oxml.ns import qn
        tbl = t._tbl
        grid = tbl.tblGrid
        cols = grid.findall(qn('w:gridCol'))
        if len(cols) == ncol:
            for k, gc in enumerate(cols):
                gc.set(qn('w:w'), str(int(widths[k] * 1440)))
        pr = tbl.tblPr
        tw = pr.find(qn('w:tblW'))
        if tw is None:
            from docx.oxml import OxmlElement
            tw = OxmlElement('w:tblW')
            pr.append(tw)
        tw.set(qn('w:w'), str(int(total_in * 1440)))
        tw.set(qn('w:type'), 'dxa')
        lay = pr.find(qn('w:tblLayout'))
        if lay is None:
            from docx.oxml import OxmlElement
            lay = OxmlElement('w:tblLayout')
            pr.append(lay)
        lay.set(qn('w:type'), 'fixed')
    d.save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paper', default='.')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    root = Path(a.paper).resolve()
    s = (root / 'main.tex').read_text(encoding='utf8')
    aux = (root / 'main.aux').read_text(encoding='utf8')

    # number macros
    macros = {}
    for f in sorted(root.glob('numbers_*.tex')):
        t = f.read_text(encoding='utf8')
        for m in re.finditer(r'\\newcommand\{\\([A-Za-z]+)\}', t):
            macros[m.group(1)] = braced(t, m.end())[0]
    s = re.sub(r'\\input\{numbers_[a-z_]+\.tex\}[^\n]*\n', '', s)
    s = re.sub(r'\\input\{(tables/[^}]+)\}', lambda m: (root / m.group(1)).read_text(encoding='utf8')
               if m.group(1).endswith('.tex') else (root / (m.group(1) + '.tex')).read_text(encoding='utf8'), s)
    names = sorted(macros, key=len, reverse=True)
    pat = re.compile(r'\\(' + '|'.join(names) + r')(?![A-Za-z])(\{\})?')
    for _ in range(3):
        s = pat.sub(lambda m: macros[m.group(1)], s)

    # cross-references and citations from main.aux
    labels = {m.group(1): m.group(2) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', aux)}
    cites = {m.group(1): int(m.group(2)) for m in re.finditer(r'\\bibcite\{([^}]+)\}\{(\d+)\}', aux)}
    s = re.sub(r'\\ref\{([^}]+)\}', lambda m: labels.get(m.group(1), '??'), s)

    def cite(arg):
        ns = sorted(cites[k.strip()] for k in arg.split(','))
        return '[' + ', '.join(map(str, ns)) + ']'
    s = sub_command(s, 'cite', cite)

    # preamble and title block: replace by plain paragraphs
    body = s[s.index('\\begin{document}') + len('\\begin{document}'):s.index('\\end{document}')]
    title = braced(s, s.index('\\title{') + 6)[0].replace('\\textbf', '')
    title = re.sub(r'[{}]', '', title)
    head = ('\\documentclass{article}\n\\begin{document}\n'
            f'\\section*{{{title}}}\n\n'
            "Aaron O'Brien\\textsuperscript{1,*}, C\\'esar Mar\\'in\\textsuperscript{2}, Leho Tedersoo\\textsuperscript{3}\n\n"
            "\\textsuperscript{1}Centro de Biotecnolog\\'ia de Sistemas (CSB), Universidad Andr\\'es Bello, Santiago, Chile. "
            "\\textsuperscript{2}Centro de Investigaci\\'on e Innovaci\\'on para el Cambio Clim\\'atico (CiiCC), Facultad de "
            "Ciencias, Universidad Santo Tom\\'as, Valdivia, Chile. "
            '\\textsuperscript{3}Mycology and Microbiology Center, Institute of Technology, University of Tartu, Tartu, '
            'Estonia.\n\n'
            '\\textsuperscript{*}Corresponding author: obrien.aaron.g@gmail.com\n\n'
            )
    body = body.replace('\\maketitle', '')
    body = body.replace('\\begin{abstract}', '\\subsection*{Abstract}').replace('\\end{abstract}', '')

    # simple commands
    body = body.replace('{\\small\\begin{longtable}', '\\begin{longtable}').replace('\\end{longtable}}', '\\end{longtable}')
    body = sub_command(body, 'spn', lambda x: '\\emph{' + x + '}')
    body = sub_command(body, 'paragraph', lambda x: '\\textbf{' + x + '}')
    body = sub_command(body, 'verify', lambda x: '[VERIFY: ' + x + ']')
    body = sub_command(body, 'todo', lambda x: '[TODO: ' + x + ']')
    body = sub_command(body, 'path', lambda x: '\\texttt{' + x.replace('_', '\\_') + '}')
    for cmd in ('clearpage', 'noindent', 'centering', 'small', 'footnotesize', 'raggedright'):
        body = re.sub(r'\\' + cmd + r'(?![A-Za-z])', '', body)
    body = re.sub(r'\\vspace\{[^}]*\}', '', body)
    body = re.sub(r'\\setlength\{\\tabcolsep\}\{[^}]*\}', '', body)
    body = re.sub(r'\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}', '', body)
    body = re.sub(r'^\\renewcommand\{\\the[A-Za-z]+\}.*$', '', body, flags=re.M)
    body = re.sub(r'\\setcounter\{[a-z]+\}\{\d+\}', '', body)
    body = re.sub(r'[Lp]\{[0-9.]+\\linewidth\}', 'l', body)       # widths are set after conversion

    # tables: caption with its number, notes after the table
    def number(t):
        m = re.search(r'\\label\{([^}]+)\}', t)
        return labels[m.group(1)]

    def fix_table(m):
        t = m.group(0)
        cap = ''
        if '\\caption{' in t:
            i = t.index('\\caption{')
            cap, j = braced(t, i + 8)
            t = t[:i] + t[j:]
        notes = ''
        if '\\begin{tablenotes}' in t:
            i, j = t.index('\\begin{tablenotes}'), t.index('\\end{tablenotes}')
            notes = t[i + len('\\begin{tablenotes}'):j].replace('\\item', '').strip()
            t = t[:i] + t[j + len('\\end{tablenotes}'):]
        t = t.replace('\\begin{threeparttable}', '').replace('\\end{threeparttable}', '')
        t = re.sub(r'\\begin\{table\}(\[[^\]]*\])?', '', t).replace('\\end{table}', '')
        n = number(m.group(0))
        t = re.sub(r'\\label\{[^}]+\}', '', t)
        out = f'\n\n\\textbf{{Table {n}.}} {cap}\n\n{t}\n\n'
        if notes:
            out += f'{{\\small {notes}}}\n\n'
        return out
    body = re.sub(r'\\begin\{table\}.*?\\end\{table\}', fix_table, body, flags=re.S)

    def fix_longtable(m):
        t = m.group(0)
        n = number(t)
        t = re.sub(r'\\label\{[^}]+\}', '', t)
        i = t.index('\\caption{')
        cap, j = braced(t, i + 8)
        t = t[:i] + t[j:]
        t = re.sub(r'\\endfirsthead.*?\\endhead', '', t, flags=re.S)          # keep the first header only
        t = t.replace('\\endfoot', '').replace('\\begin{longtable}', '\\begin{tabular}').replace('\\end{longtable}',
                                                                                                    '\\end{tabular}')
        t = re.sub(r'(\\begin\{tabular\}\{[^}]*\})\s*\\\\', r'\1', t)
        return f'\n\n\\textbf{{Table {n}.}} {cap}\n\n{t}\n\n'
    body = re.sub(r'\\begin\{longtable\}.*?\\end\{longtable\}', fix_longtable, body, flags=re.S)

    # figures: PNG files and numbered captions
    tmp = Path(tempfile.mkdtemp())

    def fix_figure(m):
        t = m.group(0)
        g = re.search(r'\\includegraphics(\[[^\]]*\])?\{([^}]+)\}', t)
        pdf = root / g.group(2)
        png = tmp / (pdf.stem + '.png')
        subprocess.run(['pdftoppm', '-r', '200', '-png', '-singlefile', str(pdf), str(png.with_suffix(''))], check=True)
        i = t.index('\\caption{')
        cap = braced(t, i + 8)[0]
        n = number(t)
        return (f'\n\n\\begin{{figure}}\n\\includegraphics[width=6in]{{{png}}}\n\\caption{{Figure {n}. {cap}}}\n'
                f'\\end{{figure}}\n\n')
    body = re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', fix_figure, body, flags=re.S)

    body = re.sub(r'\\label\{[^}]+\}', '', body)

    # bibliography as numbered paragraphs
    def fix_bib(m):
        t = m.group(0)
        items = re.split(r'\\bibitem\{([^}]+)\}', t)[1:]
        out = ['\n\\section*{References}\n']
        for k, txt in zip(items[0::2], items[1::2]):
            txt = txt.replace('\\end{thebibliography}', '').strip()
            out.append(f'[{cites[k]}] {txt}\n')
        return '\n\n'.join(out)
    body = re.sub(r'\\begin\{thebibliography\}\{\d+\}.*?\\end\{thebibliography\}', fix_bib, body, flags=re.S)

    flat = head + body + '\n\\end{document}\n'
    flat_path = root / 'main_flat_for_docx.tex'
    flat_path.write_text(flat, encoding='utf8')
    subprocess.run(['pandoc', str(flat_path), '-f', 'latex', '-t', 'docx', '--number-sections', '-o', a.out], check=True)
    fit_tables(a.out)
    print(f'wrote {a.out} ({len(cites)} references, {len(labels)} labels)')


if __name__ == '__main__':
    main()
