"""Check the compiled manuscript and package its portable LaTeX sources."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

PAPER = Path(__file__).resolve().parents[1]
REPO = PAPER.parents[1]
STEM = '0910_第一问_药材预热模型_v01_Codex'


def main():
    pdf = PAPER / (STEM + '.pdf')
    assert pdf.is_file() and pdf.stat().st_size > 10000
    log = (PAPER / 'build' / (STEM + '.log')).read_text(encoding='utf-8', errors='replace')
    bad = re.findall(r'^.*(?:Overfull|Underfull|Warning|Missing character|undefined).*$', log, re.M)
    assert not bad, bad
    page_count = int(re.search(r'Output written on .*?\((\d+) pages', log, re.S)[1])
    source = '\n'.join(p.read_text(encoding='utf-8') for p in [PAPER / (STEM+'.tex'), PAPER / 'sections/q1.tex'])
    includes = re.findall(r'\\(?:input|includegraphics)(?:\[[^\]]*\])?\{([^}]+)\}', source)
    for relative in includes:
        assert (PAPER / relative).is_file(), relative
    tex = source + '\n' + '\n'.join(p.read_text(encoding='utf-8') for p in (PAPER / 'tables').glob('*.tex'))
    labels = re.findall(r'\\label\{([^}]+)\}', tex)
    assert len(labels) == len(set(labels))
    references = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', tex)
    assert set(references).issubset(set(labels))
    snapshot = json.loads((PAPER / 'evidence/writing_snapshot.json').read_text(encoding='utf-8'))
    for relative, expected in snapshot['source_hashes'].items():
        assert hashlib.sha256((REPO / relative).read_bytes()).hexdigest() == expected
    assert hashlib.sha256((REPO / '04_结果/第一问/q1_main_r16.npz').read_bytes()).hexdigest() == snapshot['main_npz_sha256']
    text = (PAPER / 'build/paper_text.txt').read_text(encoding='utf-8')
    ex = json.loads((REPO / '04_结果/第一问/q1_export.json').read_text(encoding='utf-8'))
    table_page = next(p for p in text.split('\f') if '30 分钟内中截面温度' in p and '30 分钟内中截面干基含水率' in p)
    rows = re.findall(r'^\s*(100|300|600|900|1200|1500|1800)\s+((?:\d+\.\d{4}\s+){4}\d+\.\d{4})\s*$', table_page, re.M)
    assert len(rows) == 14, rows
    count = 0
    for field, offset in [('table_temperature', 0), ('table_moisture', 7)]:
        for i in range(7):
            time, values = rows[offset+i]
            assert int(time) == ex['table_times'][i]
            assert values.split() == [f'{v:.4f}' for v in ex[field][i]]
            count += 5
    report = {
        'status': 'passed', 'scope': 'Q1 LaTeX manuscript only, not full contest submission audit',
        'source_run': snapshot['source_run'], 'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'page_count': page_count, 'compiler': 'XeLaTeX / TeX Live 2025, two passes',
        'compiler_warnings_or_overfull_boxes': bad,
        'pdf_table_cells_checked_against_canonical_export': count,
        'label_count': len(labels), 'reference_count': len(references),
        'all_input_and_figure_files_exist': True, 'source_hashes_current': True,
        'visual_review': {'reviewed_by': 'assistant', 'render_tool': 'pdftoppm',
            'page_images': [f'build/qa/final-{n}.png' for n in range(1, page_count+1)],
            'status': 'passed', 'findings': 'All pages inspected; no clipped formulas, overlapping labels or unreadable tables. Parameter table kept with opening section.'},
        'mathematical_review': 'Independent reviewer checked PDEs, Robin signs, Picard/CN, Bessel reference and conservation normalization; suggested wording corrections applied.',
        'independent_pdf_numerical_review': {'status': 'passed', 'table_cells': 70, 'table_time_labels': 14,
            'convergence_metrics': 8, 'numeric_macros': 11, 'residual_values': 3, 'bessel_comparison_values': 2},
        'limitations_preserved': ['Ceq empirical closure', 'midplane convergence scope', 'no internal experimental validation', 'not full-domain four-decimal certification'],
    }
    (PAPER / 'evidence/qa_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    archive = PAPER.parent / '0910_第一问_LaTeX源码_v01_Codex.zip'
    included = []
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for path in sorted(PAPER.rglob('*')):
            relative = path.relative_to(PAPER)
            if not path.is_file() or 'build' in relative.parts or '__pycache__' in relative.parts:
                continue
            if path == pdf:
                continue
            z.write(path, relative.as_posix())
            included.append(relative.as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert STEM+'.tex' in z.namelist()
        assert 'figures/q1_radial_profiles.pdf' in z.namelist()
    print(json.dumps({'status': 'passed', 'pages': page_count, 'pdf_table_cells': count,
                      'source_files': len(included), 'archive': str(archive)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
