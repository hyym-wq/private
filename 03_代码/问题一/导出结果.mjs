import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';

const here = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(here,'../../work/表格数据.json'),'utf8'));
const book = await SpreadsheetFile.importXlsx(await FileBlob.load(data.模板));
if (process.argv.includes('--preview')) {
  console.log((await book.inspect({kind:'workbook,sheet,table',maxChars:2000,tableMaxRows:5,tableMaxCols:6})).ndjson);
  const blob = await book.render({sheetName:'温度',range:'A1:F5',scale:1.5,format:'png'});
  await fs.writeFile(path.join(here,'../../work/模板预览.png'),new Uint8Array(await blob.arrayBuffer()));
} else {
  for (const name of ['温度','水分浓度']) {
    const sheet = book.worksheets.getItem(name);
    // 展开模板中省略号所代表的时间、半径；保留两张原表和A1原标签。
    sheet.getRange('B1:V1').copyFrom(sheet.getRange('B1'),'all');
    sheet.getRange('A2:A1801').copyFrom(sheet.getRange('A2'),'all');
    sheet.getRange('B2:V1801').copyFrom(sheet.getRange('B2'),'all');
    sheet.getRange('B1:V1').values = [data.半径_cm];
    sheet.getRange('A2:V1801').values = data[name];
    sheet.getRange('A2:A1801').setNumberFormat('0');
    sheet.getRange('B1:V1').setNumberFormat('0.0');
    sheet.getRange('B2:V1801').setNumberFormat('0.0000');
    sheet.getRange('B1:V1801').format.columnWidth = 11;
    sheet.freezePanes.freezeRows(1);
  }
  book.recalculate();
  for (const name of ['温度','水分浓度']) {
    console.log((await book.inspect({kind:'table',range:`${name}!A1799:V1801`,include:'values,formulas',tableMaxRows:3,tableMaxCols:5,maxChars:1800})).ndjson);
    const blob = await book.render({sheetName:name,range:'A1:G8',scale:1.5,format:'png'});
    await fs.writeFile(path.join(here,`../../work/${name}结果预览.png`),new Uint8Array(await blob.arrayBuffer()));
  }
  console.log((await book.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!',options:{useRegex:true,maxResults:10},maxChars:1000})).ndjson);
  const out = path.join(here,'../../04_结果/问题一/result1.xlsx');
  await (await SpreadsheetFile.exportXlsx(book)).save(out);
  console.log('已保存',out);
}
