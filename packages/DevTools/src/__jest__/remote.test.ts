/**
 * @jest-environment jsdom
 */

import * as utils from './test-node';
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
  await browser?.close();
});

expect.extend({
  checkOutput(received, expected, context) {
    if (received === expected) {
      return {
        message: () => context,
        pass: true
      };
    } else {
      return {
        message: () => context,
        pass: false
      };
    }
  }
});

it('TEST', async () => {
  const target_package:string = process.env.TARGET_PACKAGE ?? 'DevTools';
  console.log(`Testing ${target_package} package`);

  let r = await page.evaluate((target_package):Promise<object> => {
    return new Promise<object>((resolve, reject) => {
      (<any>window).grok.functions.eval(target_package + ':test()').then((df: any) => {
        let cStatus = df.columns.byName('success');
        let cMessage = df.columns.byName('result');
        let cCat = df.columns.byName('category');
        let cName = df.columns.byName('name');
        let failed = false;
        let passReport = '';
        let failReport = '';
        for (let i = 0; i < df.rowCount; i++) {
          if (cStatus.get(i)) {
            passReport += `Test result : ${target_package}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
          } else {
            failed = true;
            failReport += `Test result : ${target_package}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
          }
        }
        resolve({failReport, passReport, failed});
      }).catch((e: any) => reject(e));
    });
  }, target_package);
  // @ts-ignore
  console.log(r.passReport);
  // @ts-ignore
  expect(r.failed).checkOutput(false, r.failReport);
}, 100000);
