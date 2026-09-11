from pathlib import Path
import zipfile
p=Path(__file__).resolve().parents[2]/'results/Q3/result3.xlsx'
with zipfile.ZipFile(p) as z:
    x=z.read('xl/worksheets/sheet1.xml').decode('utf8')
    print(p, 'rows', x.count('<row '), 'cells', x.count('<c '), 'valid_zip', z.testzip() is None)
    assert x.count('<row ') == 3432
    assert x.count('<c ') == 3432*22
