/**
 * @jest-environment jsdom
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as utils from 'datagrok-tools/src/test-utils';
// @ts-ignore
import * as publish from 'datagrok-tools/bin/commands/publish';
import puppeteer from 'puppeteer';
// import * as pkg from '../package';

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

async function testSearchSubstructure(params: any | null = null) {
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

/*
it('My test', async () => {
  const testOutput = await page.evaluate(async (_package: any) => {
    const df = await _package.testFunc();
    const countResult = df.rowCount;
    return {countResult};
  }, pkg._package);
  expect(testOutput.countResult).toBe(3);
});
*/

it('searchSubstructure graph', async () => { await testSearchSubstructure({substructLibrary: false}); });
it('searchSubstructure library 1 (explicit parameter)', async () => { await testSearchSubstructure({substructLibrary: true}); });
it('searchSubstructure library 2 (implicit parameter)', async () => { await testSearchSubstructure({}); });
it('searchSubstructure library 3 (no parameter)', async () => { await testSearchSubstructure(); });