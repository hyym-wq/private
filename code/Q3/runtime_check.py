import importlib.util
import sys
from pathlib import Path
if len(sys.argv) > 1:
    sys.stdout.reconfigure(encoding='utf-8')
    if sys.argv[1] == 'read':
        p = Path(sys.argv[2]); lines = p.read_text(encoding='utf-8').splitlines()
        lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        hi = int(sys.argv[4]) if len(sys.argv) > 4 else len(lines)
        print('\n'.join(lines[lo:hi])); raise SystemExit
    if sys.argv[1] == 'paths':
        base = Path('C:/Users/29945/.codex/plugins/cache/openai-primary-runtime/spreadsheets/26.909.12148/skills/spreadsheets')
        print(list(base.rglob('mark_artifact_operation_started.mjs')))
        print(Path('C:/Users/29945/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool').exists())
        raise SystemExit
print(sys.executable, flush=True)
for m in ['numpy','scipy','numba','openpyxl','psutil']:
    print(m, bool(importlib.util.find_spec(m)), flush=True)
try:
    import psutil
    for p in psutil.process_iter(['pid','name','cmdline']):
        if (p.info['name'] or '').lower() == 'python.exe': print(p.info,flush=True)
except ImportError:
    pass
