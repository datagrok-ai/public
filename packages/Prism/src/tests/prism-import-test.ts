import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {parsePrismFile} from '../prism/prism-parser';
import {prismSheetToDataFrame, prismAnalysisToDataFrame} from '../prism/prism-to-dataframe';

import JSZip from 'jszip';


/** Creates a .prism file with y_sd_n format for testing. */
async function createSdnPrismFile(): Promise<Uint8Array> {
  const zip = new JSZip();

  const sheetId = 'sheet-sdn';
  const tableId = 'tbl-sdn';
  const dsX = 'ds-x';
  const dsY = 'ds-y';

  zip.file('document.json', JSON.stringify({
    sheets: {data: [sheetId]},
    sheetAttributesMap: {[sheetId]: {title: 'SD+N Data'}},
  }));

  zip.file(`data/sheets/${sheetId}/sheet.json`, JSON.stringify({
    title: 'SD+N Data',
    table: {
      uid: tableId,
      dataFormat: 'y_sd_n',
      dataSets: [dsY],
      xDataSet: dsX,
    },
  }));

  // X, Mean, SD, N
  zip.file(`data/tables/${tableId}/data.csv`, '1,10,1.5,5\n2,20,2.1,5\n3,30,3.0,5\n');
  zip.file(`data/sets/${dsX}.json`, JSON.stringify({uid: dsX, title: 'Time'}));
  zip.file(`data/sets/${dsY}.json`, JSON.stringify({uid: dsY, title: 'Response'}));

  return await zip.generateAsync({type: 'uint8array'});
}

/** Creates a .prism file with analysis results. */
async function createAnalysisPrismFile(): Promise<Uint8Array> {
  const zip = new JSZip();

  const dataSheetId = 'sheet-data';
  const analysisSheetId = 'sheet-analysis';
  const dataTableId = 'tbl-data';
  const resultRefId = 'result-ref';
  const resultDataSheetId = 'result-data-sheet';
  const resultTableId = 'result-tbl';
  const dsX = 'ds-x';
  const dsY = 'ds-y';
  const dsRes1 = 'ds-res1';
  const dsRes2 = 'ds-res2';

  zip.file('document.json', JSON.stringify({
    sheets: {
      data: [dataSheetId],
      analyses: [analysisSheetId],
    },
    sheetAttributesMap: {
      [dataSheetId]: {title: 'Raw Data'},
      [analysisSheetId]: {title: 'Curve Fit'},
    },
  }));

  // Data sheet
  zip.file(`data/sheets/${dataSheetId}/sheet.json`, JSON.stringify({
    title: 'Raw Data',
    table: {uid: dataTableId, dataFormat: 'y_single', dataSets: [dsY], xDataSet: dsX},
  }));
  zip.file(`data/tables/${dataTableId}/data.csv`, '1,10\n2,20\n3,30\n');
  zip.file(`data/sets/${dsX}.json`, JSON.stringify({uid: dsX, title: 'X'}));
  zip.file(`data/sets/${dsY}.json`, JSON.stringify({uid: dsY, title: 'Y'}));

  // Analysis sheet (at analyses/<id>/sheet.json)
  zip.file(`analyses/${analysisSheetId}/sheet.json`, JSON.stringify({
    title: 'Curve Fit',
    preferredResultSheet: resultRefId,
  }));

  // Result sheet reference (at analyses/<id>/result_sheets/<resultId>.json)
  zip.file(`analyses/${analysisSheetId}/result_sheets/${resultRefId}.json`, JSON.stringify({
    dataSheet: resultDataSheetId,
  }));

  // Result data sheet (at data/sheets/<id>/sheet.json)
  zip.file(`data/sheets/${resultDataSheetId}/sheet.json`, JSON.stringify({
    title: 'Fit Results',
    table: {uid: resultTableId, dataFormat: 'text', dataSets: [dsRes1, dsRes2]},
  }));
  zip.file(`data/tables/${resultTableId}/data.csv`, 'IC50,0.5\nHill Slope,1.2\nR squared,0.98\n');
  zip.file(`data/sets/${dsRes1}.json`, JSON.stringify({uid: dsRes1, title: 'Parameter'}));
  zip.file(`data/sets/${dsRes2}.json`, JSON.stringify({uid: dsRes2, title: 'Value'}));

  return await zip.generateAsync({type: 'uint8array'});
}


category('Prism Import', () => {
  test('y_sd_n format: correct column structure', async () => {
    const bytes = await createSdnPrismFile();
    const prism = await parsePrismFile(bytes);
    const df = prismSheetToDataFrame(prism.sheets[0]);

    expect(df.rowCount, 3);
    // X + Mean + SD + N = 4 columns
    expect(df.columns.length, 4);
    expect(df.columns.byIndex(0).name, 'Time');
    expect(df.columns.byIndex(1).name, 'Response_Mean');
    expect(df.columns.byIndex(2).name, 'Response_SD');
    expect(df.columns.byIndex(3).name, 'Response_N');
  });

  test('analysis results: parsed correctly', async () => {
    const bytes = await createAnalysisPrismFile();
    const prism = await parsePrismFile(bytes);

    expect(prism.analyses.length, 1);
    expect(prism.analyses[0].title, 'Curve Fit');
    expect(prism.analyses[0].resultData.length, 3);

    const df = prismAnalysisToDataFrame(prism.analyses[0]);
    expect(df.rowCount, 3);
    expect(df.columns.length >= 2, true);
  });

  test('full import: returns all DataFrames', async () => {
    const bytes = await createAnalysisPrismFile();
    const prism = await parsePrismFile(bytes);

    const results = [];
    for (const sheet of prism.sheets)
      results.push(prismSheetToDataFrame(sheet));
    for (const analysis of prism.analyses)
      results.push(prismAnalysisToDataFrame(analysis));

    expect(results.length, 2); // 1 data sheet + 1 analysis
    expect(results[0].name, 'Raw Data');
  });
});
