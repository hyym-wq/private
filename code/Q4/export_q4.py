from pathlib import Path
import numpy as np, json, zipfile

root = Path(__file__).resolve().parents[2]
out = root / 'results/Q4'
z = np.load(out / 'q4_ALE_N640_tol1e-09_full.npz')
times, xi, Rs, C = z['times'], z['xi'], z['R'], z['C']
targets = np.arange(0, 2.0001, .1)
assert C.shape[0] == len(times) and times[-1] == 182968

rows_data = []
for i, t in enumerate(times):
    row = []
    for q in targets:
        row.append(float(np.interp(q/100/Rs[i], xi, C[i])) if q/100 <= Rs[i] + 1e-14 else None)
    row.append(float(C[i, -1]))
    rows_data.append(row)

def col(n):
    s = ''; n += 1
    while n:
        n, rem = divmod(n-1, 26); s = chr(65+rem) + s
    return s
def cell(ref, val):
    if val is None: return ''
    if isinstance(val, str): return f'<c r="{ref}" t="inlineStr"><is><t>{val}</t></is></c>'
    return f'<c r="{ref}"><v>{float(val):.12g}</v></c>'

header = ['时间\\到药材中心的距离'] + [f'{x:.1f}' for x in targets] + ['药材表面']
rows = ['<row r="1">' + ''.join(cell(col(i), v) for i, v in enumerate(header)) + '</row>']
for ii, (t, row) in enumerate(zip(times, rows_data), 2):
    rows.append(f'<row r="{ii}">' + cell('A'+str(ii), t) + ''.join(cell(col(i+1)+str(ii), v) for i, v in enumerate(row)) + '</row>')
sheet = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + ''.join(rows) + '</sheetData></worksheet>'
ct = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>'
rels = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'
wb = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
wr = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
path = out / 'result4.xlsx'
with zipfile.ZipFile(path, 'w') as f:
    f.writestr('[Content_Types].xml', ct); f.writestr('_rels/.rels', rels); f.writestr('xl/workbook.xml', wb); f.writestr('xl/_rels/workbook.xml.rels', wr); f.writestr('xl/worksheets/sheet1.xml', sheet)

table_times = np.r_[np.arange(21600, 182968, 21600), 182968]
ti = [int(np.flatnonzero(times == t)[0]) for t in table_times]
table = [[times[i], *rows_data[i][::5], rows_data[i][-1]] for i in ti]
csv = ['时间/h,0 cm,0.5 cm,1 cm,1.5 cm,2 cm,药材表面']
csv += [f'{r[0]/3600:.6f},' + ','.join('' if x is None else f'{x:.8f}' for x in r[1:]) for r in table]
(out / 'table6.csv').write_text('\n'.join(csv)+'\n', encoding='utf-8-sig')
summary = {'method':'Q2 framework + material-coordinate shrinkage + Appendix 4 + Kirchhoff flux + adaptive BDF','N':640,'rtol':1e-9,'max_step_s':180,'continuous_event_s':182967.2685544065,'first_integer_dry_s':182968,'drying_h':182968/3600,'radius_end_cm':float(Rs[-1]*100),'previous_max_C':0.1500001387526868,'terminal_max_C':0.1499996220903369,'output_rows':len(times),'out_of_domain_policy':'blank','table':table,'mesh_convergence_h':{'160':50.8248103873,'320':50.8243523544,'640':50.8242412651}}
(out / 'q4_final_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(path, len(times), len(header))
