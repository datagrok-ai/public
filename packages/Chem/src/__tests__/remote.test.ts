/**
 * @jest-environment jsdom
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as utils from '@datagrok-libraries/utils/src/test-node';
// @ts-ignore
// import * as publish from 'datagrok-tools/bin/commands/publish';
import puppeteer from 'puppeteer';

const START_TIMEOUT: number = 100000;
const DEFAULT_NAVIGATION_TIMEOUT = 60000;
let browser: puppeteer.Browser;
let page: puppeteer.Page;

beforeAll(async () => {
  const out = await utils.getBrowserPage(puppeteer);
  browser = out.browser;
  page = out.page;
  page.setDefaultNavigationTimeout(DEFAULT_NAVIGATION_TIMEOUT);
}, START_TIMEOUT);

afterAll(async () => {
  await browser.close();
});

//todo: write jest adapter to run internal tests