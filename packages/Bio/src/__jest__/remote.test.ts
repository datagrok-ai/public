/**
 * @jest-environment jsdom
 */

import * as utils from './test-node';
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

it('TEST', async () => {
  const targetPackage: string = process.env.TARGET_PACKAGE ?? 'Bio';
  console.log(`Testing ${targetPackage} package`);

  //console.log(require('root-require')('package.json').version);
  const r = await page.evaluate((targetPackage): Promise<object> => {
    return new Promise<object>((resolve, reject) => {
      (<any>window).grok.functions.eval(targetPackage + ':test()').then((df: any) => {
        const cStatus = df.columns.byName('success');
        const cMessage = df.columns.byName('result');
        const cCat = df.columns.byName('category');
        const cName = df.columns.byName('name');
        let failed = false;
        let report = '';
        for (let i = 0; i < df.rowCount; i++) {
          if (!cStatus.get(i)) {
            report += `${cCat.get(i)}.${cName.get(i)}: ${cMessage.get(i)}\n`;
            failed = true;
          }
        }
        resolve({report, failed});
      }).catch((e: any) => reject(e));
    });
  }, targetPackage);
  // @ts-ignore
  console.log(r.report);
  // @ts-ignore
  expect(r.failed).toBe(false);
}, 100000);

// it('WebLogo.getAlphabetSimilarity', () => {
//
// });
