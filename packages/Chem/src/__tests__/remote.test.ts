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

it('searchSubstructure graph', async () => {
  const testOutput = await page.evaluate(async () => {
    const df = await grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv'); // TODO: localize
    const col = df.col('smiles')!;
    let bitset = await grok.chem.searchSubstructure(
      col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', {substructLibrary: false});
    const countDataframe = col.length;
    const countResult = bitset.trueCount;
    return {countDataframe, countResult};
  });
  expect(testOutput.countDataframe).toBe(200);
  expect(testOutput.countResult).toBe(90);
});