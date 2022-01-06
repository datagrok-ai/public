/**
 * @jest-environment jsdom
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as utils from '@datagrok-libraries/utils/src/test-utils';
// @ts-ignore
// import * as publish from 'datagrok-tools/bin/commands/publish';
import puppeteer from 'puppeteer';

const P_START_TIMEOUT: number = 100000;
let browser: puppeteer.Browser;
let page: puppeteer.Page;

beforeAll(async () => {
  const out = await utils.getBrowserPage(puppeteer);
  browser = out.browser;
  page = out.page;
}, P_START_TIMEOUT);

afterAll(async () => {
  await browser.close();
});

const fs = require('fs');
function requireTextRequire(name: string, require: any) {
  return fs.readFileSync(require.resolve(name)).toString();
};
function requireText(name: string) {
  return requireTextRequire(name, require);
}

// 5-cyclohexyl-1,3-dihydro-2H-benzo[e][1,4]diazepin-2-one
// 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'

async function _testSearchSubstructure(csv: string, pattern: string, trueIndices: number[], params: any | null = null) {
  const testOutput = await page.evaluate(async (csv: string, pattern: string, trueIndices: number[], params: any) => {
    // await grok.functions.call('Chem:initChem');
    const df = DG.DataFrame.fromCsv(csv);
    const col = df.columns[0];
    const bitset: DG.BitSet = (params !== null) ?
      (await grok.chem.searchSubstructure(col, pattern, params)) :
      (await grok.chem.searchSubstructure(col, pattern));
    const bitsetString = bitset.toBinaryString();
    return {bitsetString};
  }, csv, pattern, trueIndices, params);
  const bitset = [...testOutput.bitsetString];
  for (let k = 0; k < trueIndices.length; ++k) {
    expect(bitset[trueIndices[k]] === '1').toBe(true);
    bitset[trueIndices[k]] = '0';
  }
  for (let i = 0; i < bitset.length; ++i) {
    expect(bitset[i] === '0').toBe(true);
  }
}

async function _testSearchSubstructureSARSmall(params: any | null = null) {
  const fileSARSmall = requireText('./sar_small.csv');
  const testOutput = await page.evaluate(async (file: string, params: any) => {
    // await grok.functions.call('Chem:initChem');
    const df = DG.DataFrame.fromCsv(file);
    const col = df.columns[0];
    const bitset: DG.BitSet = (params !== null) ?
      (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', params)) :
      (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'));
    const countDataframe = col.length;
    const countResult = bitset.trueCount;
    return {countDataframe, countResult};
  }, fileSARSmall, params);
  expect(testOutput.countDataframe).toBe(200);
  expect(testOutput.countResult).toBe(90);
}

async function _testSearchSubstructureAllParameters(foo: any) {
  await foo({substructLibrary: false});
  await foo({substructLibrary: true});
  await foo({});
  await foo();
}

// 90 hits out of 200 on sar_small.csv
it('chem.searchSubstructure.sar_small', async () => {
  await _testSearchSubstructureAllParameters(
    _testSearchSubstructureSARSmall);
});

// Number of molecules is smaller than a number of threads
it('chem.searchSubstructure.5_rows', async () => {
  const targetSmiles = [
    'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
    'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
  ];
  const substructureSmiles = 'C1CC1';
  const csv =
`smiles
${targetSmiles[0]}
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
${targetSmiles[1]}
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5
`;
  await _testSearchSubstructureAllParameters(
    async (params: any) => await _testSearchSubstructure(
      csv, substructureSmiles, [0, 2], params,
    ),
  );
});

it('chem.findSimilar.sar_small', async () => {
  const fileSARSmall = requireText('./sar_small.csv');
  const testOutput = await page.evaluate(async (file: string) => {
    //await grok.functions.call('Chem:initChem');
    const dfInput = DG.DataFrame.fromCsv(file);
    const colInput = dfInput.columns[0];
    const dfResult: DG.DataFrame = // shouldn't be null
      (await grok.chem.findSimilar(colInput, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))!;
    const numRowsOriginal = dfInput.rowCount;
    const numRows = dfResult.rowCount;
    const columnNames = [
      dfResult.columns[0].name,
      dfResult.columns[1].name,
      dfResult.columns[2].name,
    ];
    const first5Rows: any[] = [];
    for (let i = 0; i < 5; ++i) {
      const molecule: string = dfResult.columns[0].get(i);
      const score: number = dfResult.columns[1].get(i);
      const index: number = dfResult.columns[2].get(i);
      const obj = {molecule, score, index};
      first5Rows[i] = obj;
    }
    return {numRowsOriginal, numRows, columnNames, first5Rows};
  }, fileSARSmall);
  expect(testOutput.numRows).toBe(testOutput.numRowsOriginal);
  expect(testOutput.columnNames[0]).toBe('molecule');
  expect(testOutput.columnNames[1]).toBe('score');
  expect(testOutput.columnNames[2]).toBe('index');
  const areEqualFloat = (a: number, b: number) => Math.abs(a-b) < 0.001;
  const arr = testOutput.first5Rows;
  expect(arr[0].molecule).toBe('O=C1CN=C(c2ccccc2N1)C3CCCCC3');
  expect(arr[1].molecule).toBe('O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3');
  expect(arr[2].molecule).toBe('O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3');
  expect(arr[3].molecule).toBe('CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3');
  expect(arr[4].molecule).toBe('O=C1CN=C(c2cc(Br)ccc2N1)C3CCCCC3');
  expect(areEqualFloat(arr[0].score, 1.0000)).toBe(true);
  expect(areEqualFloat(arr[1].score, 0.7813)).toBe(true);
  expect(areEqualFloat(arr[2].score, 0.7429)).toBe(true);
  expect(areEqualFloat(arr[3].score, 0.7429)).toBe(true);
  expect(areEqualFloat(arr[4].score, 0.7429)).toBe(true);
  expect(arr[0].index).toBe(0);
  expect(arr[1].index).toBe(30);
  expect(arr[2].index).toBe(5);
  expect(arr[3].index).toBe(15);
  expect(arr[4].index).toBe(25);
});
