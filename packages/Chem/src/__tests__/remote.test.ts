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

it('searchSubstructure graph', async () => {
  let fileSARSmall = requireText("./sar_small.csv", require);
  const testOutput = await page.evaluate(async (file: string) => {
    const df = DG.DataFrame.fromCsv(file);
    const col = df.col('smiles')!;
    let bitset = await grok.chem.searchSubstructure(
      col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', {substructLibrary: false});
    const countDataframe = col.length;
    const countResult = bitset.trueCount;
    return {countDataframe, countResult};
  }, fileSARSmall);
  expect(testOutput.countDataframe).toBe(200);
  expect(testOutput.countResult).toBe(90);
});