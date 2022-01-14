/**
 * @jest-environment jsdom
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as utils from '@datagrok-libraries/utils/src/test-node';
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

it('sum', async () => {
  let r = await page.evaluate(async () => {
    return await grok.functions.eval('PowerPack:sum(1, 2)')
  });
  expect(r).toBe(3);
});

it('current user', async () => {
  let name = await page.evaluate(async () => {
    let u = await grok.dapi.users.current();
    return u.firstName;
  });
  expect(name).toBe('Alexander');
});

it('dataFrame construction test', async () => {
  const testOutput = await page.evaluate(() => {
    const df = DG.DataFrame.fromCsv(
      `columnA,columnB
        valueA-1,valueB-1
        valueA-2,valueB-2`);
    const exist = df.columns.byName('columnB') != null;
    const notExist = df.columns.byName('BlackSabbath') != null;

    return {exist, notExist};
  });
  expect(testOutput.exist).toBe(true);
  expect(testOutput.notExist).toBe(false);
});