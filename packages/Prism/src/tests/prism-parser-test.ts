import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {parsePrismFile} from '../prism/prism-parser';
import {prismSheetToDataFrame} from '../prism/prism-to-dataframe';
import {hasXYData, prismSheetToFitChartData, prismToFitDataFrame} from '../prism/prism-curves';

import JSZip from 'jszip';


/** Creates a synthetic .prism file matching the real format structure. */
async function createTestPrismFile(): Promise<Uint8Array> {
  const zip = new JSZip();

  const sheetId = 'sheet-001';
  const tableId = 'table-001';
  const xDatasetId = 'ds-x-001';
  const yDatasetId = 'ds-y-001';

  // document.json (sheetAttributesMap values are {title: "..."} objects)
  zip.file('document.json', JSON.stringify({
    sheets: {
      data: [sheetId],
      analyses: [],
    },
    sheetAttributesMap: {
      [sheetId]: {title: 'Dose Response'},
    },
  }));

  // data/sheets/<sheetId>/sheet.json
  zip.file(`data/sheets/${sheetId}/sheet.json`, JSON.stringify({
    title: 'Dose Response',
    table: {
      uid: tableId,
      dataFormat: 'y_replicates',
      replicatesCount: 2,
      dataSets: [yDatasetId],
      xDataSet: xDatasetId,
    },
  }));

  // data/tables/<tableId>/data.csv (X, Y_rep1, Y_rep2)
  zip.file(`data/tables/${tableId}/data.csv`, [
    '0.001,5.1,4.8',
    '0.01,25.3,23.1',
    '0.1,52.7,55.2',
    '1,85.4,83.9',
    '10,97.2,96.8',
  ].join('\n'));

  // data/sets/<datasetId>.json
  zip.file(`data/sets/${xDatasetId}.json`, JSON.stringify({
    uid: xDatasetId,
    title: 'Concentration',
  }));
  zip.file(`data/sets/${yDatasetId}.json`, JSON.stringify({
    uid: yDatasetId,
    title: 'Compound X',
  }));

  return await zip.generateAsync({type: 'uint8array'});
}

/** Creates a multi-sheet .prism file for testing. */
async function createMultiSheetPrismFile(): Promise<Uint8Array> {
  const zip = new JSZip();

  const sheet1Id = 'sheet-xy';
  const sheet2Id = 'sheet-col';
  const table1Id = 'tbl-xy';
  const table2Id = 'tbl-col';
  const dsX = 'ds-x';
  const dsY1 = 'ds-y1';
  const dsY2 = 'ds-y2';

  zip.file('document.json', JSON.stringify({
    sheets: {
      data: [sheet1Id, sheet2Id],
      analyses: [],
    },
    sheetAttributesMap: {
      [sheet1Id]: {title: 'XY Data'},
      [sheet2Id]: {title: 'Column Data'},
    },
  }));

  // XY sheet
  zip.file(`data/sheets/${sheet1Id}/sheet.json`, JSON.stringify({
    title: 'XY Data',
    table: {
      uid: table1Id,
      dataFormat: 'y_single',
      dataSets: [dsY1],
      xDataSet: dsX,
    },
  }));
  zip.file(`data/tables/${table1Id}/data.csv`, '1,10\n2,20\n3,30\n');
  zip.file(`data/sets/${dsX}.json`, JSON.stringify({uid: dsX, title: 'X'}));
  zip.file(`data/sets/${dsY1}.json`, JSON.stringify({uid: dsY1, title: 'Y'}));

  // Column sheet (with row titles)
  zip.file(`data/sheets/${sheet2Id}/sheet.json`, JSON.stringify({
    title: 'Column Data',
    table: {
      uid: table2Id,
      dataFormat: 'y_single',
      dataSets: [dsY2],
      rowTitlesDataSet: 'ds-rowtitles',
    },
  }));
  zip.file(`data/tables/${table2Id}/data.csv`, 'GroupA,5\nGroupB,6\nGroupC,7\n');
  zip.file(`data/sets/${dsY2}.json`, JSON.stringify({uid: dsY2, title: 'Values'}));

  return await zip.generateAsync({type: 'uint8array'});
}


category('Prism Parser', () => {
  test('parsePrismFile: parses XY sheet', async () => {
    const bytes = await createTestPrismFile();
    const prism = await parsePrismFile(bytes);
    expect(prism.sheets.length, 1);
    expect(prism.sheets[0].title, 'Dose Response');
    expect(prism.sheets[0].dataFormat, 'y_replicates');
    expect(prism.sheets[0].replicatesCount, 2);
    expect(prism.sheets[0].xTitle, 'Concentration');
    expect(prism.sheets[0].yColumnTitles.length, 1);
    expect(prism.sheets[0].yColumnTitles[0], 'Compound X');
    expect(prism.sheets[0].data.length, 5);
  });

  test('parsePrismFile: parses multiple sheets', async () => {
    const bytes = await createMultiSheetPrismFile();
    const prism = await parsePrismFile(bytes);
    expect(prism.sheets.length, 2);
    expect(prism.sheets[0].title, 'XY Data');
    expect(prism.sheets[1].title, 'Column Data');
  });

  test('parsePrismFile: handles empty archive', async () => {
    const zip = new JSZip();
    const bytes = await zip.generateAsync({type: 'uint8array'});
    const prism = await parsePrismFile(bytes);
    expect(prism.sheets.length, 0);
  });

  test('prismSheetToDataFrame: converts XY sheet', async () => {
    const bytes = await createTestPrismFile();
    const prism = await parsePrismFile(bytes);
    const df = prismSheetToDataFrame(prism.sheets[0]);
    expect(df.rowCount, 5);
    expect(df.columns.length, 3); // X + 2 replicates
    expect(df.columns.byIndex(0).name, 'Concentration');
  });

  test('prismSheetToDataFrame: converts column sheet', async () => {
    const bytes = await createMultiSheetPrismFile();
    const prism = await parsePrismFile(bytes);
    const df = prismSheetToDataFrame(prism.sheets[1]);
    expect(df.rowCount, 3);
    expect(df.columns.length, 2); // Row title + Values
  });

  test('hasXYData: detects numeric X data', async () => {
    const bytes = await createTestPrismFile();
    const prism = await parsePrismFile(bytes);
    expect(hasXYData(prism.sheets[0]), true);
  });

  test('prismSheetToFitChartData: converts to fit chart', async () => {
    const bytes = await createTestPrismFile();
    const prism = await parsePrismFile(bytes);
    const chartData = prismSheetToFitChartData(prism.sheets[0]);
    expect(chartData !== null, true);
    expect(chartData!.series!.length, 1);
    expect(chartData!.series![0].name, 'Compound X');
    expect(chartData!.series![0].points.length, 10); // 5 x-values * 2 replicates
    expect(chartData!.chartOptions!.logX, true);
  });

  test('prismToFitDataFrame: creates fit DataFrame', async () => {
    const bytes = await createTestPrismFile();
    const prism = await parsePrismFile(bytes);
    const df = prismToFitDataFrame(prism.sheets);
    expect(df.rowCount, 1);
    expect(df.col('Fitted Curve') !== null, true);
    expect(df.col('Fitted Curve')!.semType, 'fit');
  });
});
