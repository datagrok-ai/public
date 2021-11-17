/**
 * @jest-environment jsdom
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as utils from 'datagrok-tools/src/test-utils';
// @ts-ignore
import * as publish from 'datagrok-tools/bin/commands/publish';
import puppeteer from 'puppeteer';

const P_START_TIMEOUT: number = 100000;
let browser: puppeteer.Browser;
let page: puppeteer.Page;

beforeAll(async () => {
  let out = await utils.getBrowserPage(puppeteer);
  browser = out.browser;
  page = out.page;
}, P_START_TIMEOUT);

afterAll(async () => {
  await browser.close();
});

let fs = require('fs');
function requireText(name: string, require: any) {
   return fs.readFileSync(require.resolve(name)).toString();
};

async function testSearchSubstructure(csv: string, pattern: string, trueIndices: number[], params: any | null = null) {
  let testOutput = await page.evaluate(async (csv: string, pattern: string, trueIndices: number[], params: any) => {
    const df = DG.DataFrame.fromCsv(csv);
    const col = df.columns[0];
    const bitset = (params !== null) ?
      (await grok.chem.searchSubstructure(col, pattern, params)) :
      (await grok.chem.searchSubstructure(col, pattern));
    let bitsetString = bitset.toBinaryString();
    return {bitsetString};
  }, csv, pattern, trueIndices, params);
  let bitset = [...testOutput.bitsetString];
  for (let k = 0; k < trueIndices.length; ++k) {
    expect(bitset[trueIndices[k]] === '1').toBe(true);
    bitset[trueIndices[k]] = '0';
  }
  for (let i = 0; i < bitset.length; ++i) {
    expect(bitset[i] === '0').toBe(true);
  }
}

async function testSearchSubstructureSARSmall(params: any | null = null) {
  let fileSARSmall = requireText("./sar_small.csv", require);
  const testOutput = await page.evaluate(async (file: string, params: any) => {
    const df = DG.DataFrame.fromCsv(file);
    const col = df.col('smiles')!;
    const bitset = (params !== null) ?
      (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', params)) :
      (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'));
    const countDataframe = col.length;
    const countResult = bitset.trueCount;
    return {countDataframe, countResult};
  }, fileSARSmall, params);
  expect(testOutput.countDataframe).toBe(200);
  expect(testOutput.countResult).toBe(90);
}

async function testSearchSubstructureAllParameters(foo: any) {
  await foo({substructLibrary: false});
  await foo({substructLibrary: true});
  await foo({});
  await foo(); 
}

// 90 hits out of 200 on sar_small.csv
it('chem.searchSubstructure.SARSmall', async () => {
  await testSearchSubstructureAllParameters(
    testSearchSubstructureSARSmall);
});

// Number of molecules is smaller than a number of threads
it('chem.searchSubstructure.5rows', async () => {
  const targetSmiles = [
    'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
    'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5'
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
  await testSearchSubstructureAllParameters(
    async (params: any) => await testSearchSubstructure(
      csv, substructureSmiles, [0, 2], params
    )
  );
});