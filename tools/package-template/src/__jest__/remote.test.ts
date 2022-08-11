/**
 * @jest-environment jsdom
 */

import * as utils from './test-node';
import puppeteer from 'puppeteer';

const P_START_TIMEOUT: number = 3600000;
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

it('TEST', async () => {
  const target_package:string = process.env.TARGET_PACKAGE;
  console.log(`Testing ${target_package} package`);

  //console.log(require('root-require')('package.json').version);
  let r = await page.evaluate((target_package):Promise<object> => {
    return new Promise<object>((resolve, reject) => {
      (<any>window).grok.functions.eval(target_package + ':test()').then((df: any) => {
        let cStatus = df.columns.byName('success');
        let cMessage = df.columns.byName('result');
        let cCat = df.columns.byName('category');
        let cName = df.columns.byName('name');
        let failed = false;
        let report = '';
        for (let i = 0; i < df.rowCount; i++)
          if (!cStatus.get(i)) {
            report += `${cCat.get(i)}.${cName.get(i)}: ${cMessage.get(i)}\n`;
            failed = true;
          }
        resolve({report, failed});
      }).catch((e: any) => reject(e));
    });
  }, target_package);
  // @ts-ignore
  console.log(r.report);
  // @ts-ignore
  expect(r.failed).toBe(false);
}, 3600000);
