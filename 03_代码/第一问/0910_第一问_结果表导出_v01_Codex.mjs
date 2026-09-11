/**
 * 2026-09-10 第一问结果表导出 v01 Codex
 * Reads the canonical solver JSON and authors result1.xlsx with @oai/artifact-tool.
 * Runtime: bundled Node.js from load_workspace_dependencies.
 * Set Q1_ARTIFACT_RUNTIME to a writable temporary directory containing a
 * node_modules junction to the bundled runtime packages; default is
 * %TEMP%/codex-q1-workbook-0910. No repository dependency installation is needed.
 * Usage: node <this file> [--inspect-template]
 * Mark the first workbook authoring operation using the skill's marker before
 * first execution. --inspect-template is strictly read-only for source files.
 */
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath, pathToFileURL } from 'node:url';

const scriptDirectory = path.dirname(fileURLToPath(import.meta.url));
const privateDirectory = path.resolve(scriptDirectory, '..', '..');
const workspaceDirectory = path.dirname(privateDirectory);
const outputDirectory = path.join(privateDirectory, '04_结果', '第一问');
const runtimeDirectory = process.env.Q1_ARTIFACT_RUNTIME
  || path.join(os.tmpdir(), 'codex-q1-workbook-0910');
const requireRuntime = createRequire(path.join(runtimeDirectory, 'bootstrap.cjs'));
const { FileBlob, Workbook, SpreadsheetFile } = await import(
  pathToFileURL(requireRuntime.resolve('@oai/artifact-tool')).href
);

async function main() {
if (process.argv.includes('--inspect-template')) {
  const templatePath = path.join(workspaceDirectory,
    'CUMCM2026Problems', 'A题', '附件', '附件3', 'result1.xlsx');
  const sourceWorkbook = await SpreadsheetFile.importXlsx(await FileBlob.load(templatePath));
  console.log((await sourceWorkbook.inspect({
    kind: 'workbook,sheet,table', maxChars: 4500, tableMaxRows: 6,
    tableMaxCols: 6, tableMaxCellChars: 80,
  })).ndjson);
  for (const sheetName of ['温度', '水分浓度']) {
    const preview = await sourceWorkbook.render({ sheetName, range: 'A1:F5', scale: 1.5 });
    await fs.writeFile(path.join(runtimeDirectory, `template_${sheetName}.png`),
      new Uint8Array(await preview.arrayBuffer()));
  }
  return;
}

const sourcePath = path.join(outputDirectory, 'q1_export.json');
const data = JSON.parse(await fs.readFile(sourcePath, 'utf8'));
const assert = (condition, message) => { if (!condition) throw new Error(message); };
assert(data.times?.length === 1800, 'Expected 1800 output times.');
assert(data.times.every((t, i) => t === i + 1), 'Output times must be 1..1800 seconds.');
assert(data.radii_cm?.length === 21, 'Expected 21 radii.');
assert(data.radii_cm.every((r, i) => Math.abs(r - i / 10) < 1e-12),
  'Radii must span 0..2 cm in steps of 0.1 cm.');
for (const field of ['temperature', 'moisture']) {
  assert(data[field]?.length === 1800, `${field}: expected 1800 rows.`);
  assert(data[field].every(row => row.length === 21
    && row.every(v => typeof v === 'number' && Number.isFinite(v))),
  `${field}: every row must contain 21 finite numerical values.`);
}
for (const field of ['temperature', 'moisture']) {
  const table = data[`table_${field}`];
  assert(table?.length === data.table_times?.length, `${field}: invalid summary rows.`);
  data.table_times.forEach((t, i) => data.table_radii_cm.forEach((r, j) => {
    const column = Math.round(10 * r);
    assert(Math.abs(table[i][j] - data[field][t - 1][column]) < 1e-9,
      `${field}: summary and detailed result disagree at t=${t}, r=${r}.`);
  }));
}

const workbook = Workbook.create();
const fields = [
  { name: '温度', key: 'temperature', unit: '温度(°C)', color: '#194C6B' },
  { name: '水分浓度', key: 'moisture', unit: '水分浓度(kg/kg)', color: '#346049' },
];
for (const field of fields) {
  const sheet = workbook.worksheets.add(field.name);
  sheet.showGridLines = false;
  sheet.tabColor = field.color;
  const rows = [
    [`时间(s)\\距中心距离(cm)\n${field.unit}`, ...data.radii_cm],
    ...data.times.map((t, i) => [t, ...data[field.key][i]]),
  ];
  const fullRange = sheet.getRange('A1:V1801');
  fullRange.values = rows;
  fullRange.format.font = { name: 'Microsoft YaHei', size: 10, color: '#202C35' };
  fullRange.format.verticalAlignment = 'center';
  fullRange.format.rowHeight = 19;
  fullRange.format.horizontalAlignment = 'right';
  sheet.getRange('A1:A1801').format.columnWidthPx = 215;
  sheet.getRange('B1:V1801').format.columnWidthPx = 96;
  sheet.getRange('A2:A1801').setNumberFormat('0');
  sheet.getRange('B2:V1801').setNumberFormat('0.0000');
  const headers = sheet.getRange('A1:V1');
  headers.format = {
    fill: field.color,
    font: { name: 'Microsoft YaHei', size: 10, bold: true, color: '#FFFFFF' },
    rowHeight: 38,
    verticalAlignment: 'center',
    horizontalAlignment: 'center',
    borders: { insideVertical: { style: 'thin', color: '#FFFFFF' } },
  };
  sheet.getRange('A1').format.wrapText = true;
  sheet.getRange('B1:V1').setNumberFormat('0.0');
  sheet.getRange('A2:A1801').format.fill = '#F0F3F5';
  sheet.getRange('A2:A1801').format.font = { name: 'Microsoft YaHei', size: 10, color: '#526372' };
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}
workbook.recalculate();
const inspections = [];
for (const field of fields) {
  for (const range of ['A1:D4', 'S1799:V1801']) {
    inspections.push((await workbook.inspect({
      kind: 'table', range: `'${field.name}'!${range}`, include: 'values,formulas',
      tableMaxRows: 4, tableMaxCols: 4, maxChars: 2500,
    })).ndjson);
  }
}
const errors = await workbook.inspect({
  kind: 'match',
  searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',
  options: { useRegex: true, maxResults: 10 }, maxChars: 1500,
  summary: 'Q1 final error scan',
});
inspections.push(errors.ndjson);
await fs.mkdir(outputDirectory, { recursive: true });
await fs.writeFile(path.join(runtimeDirectory, 'artifact_inspection.ndjson'), inspections.join('\n'));
for (const field of fields) {
  const preview = await workbook.render({
    sheetName: field.name, range: 'A1:V10', scale: 1.25, format: 'png',
  });
  await fs.writeFile(path.join(runtimeDirectory, `result1_${field.name}_preview.png`),
    new Uint8Array(await preview.arrayBuffer()));
}
const outputPath = path.join(outputDirectory, 'result1.xlsx');
const file = await SpreadsheetFile.exportXlsx(workbook);
await file.save(outputPath);
console.log(JSON.stringify({ output: outputPath, sheets: fields.map(f => f.name),
  shape_per_sheet: [1801, 22], data_values: 75600, number_format: '0.0000',
  previews: runtimeDirectory, error_scan: errors.ndjson }, null, 2));
}

await main();
