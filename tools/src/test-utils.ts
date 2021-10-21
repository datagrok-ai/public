import puppeteer from 'puppeteer';
import path from "path";
import os from "os";
import fs from "fs";
// @ts-ignore
import * as yaml from 'js-yaml';
// @ts-ignore
import * as utils from '../bin/utils'

export async function getToken(devKey: string) {
  return 'GGuDtdZ1Qu3sHOwscchQdl2ozmzXwzcjJ52VByCyYqOQS6v5XzMHcLtOs8ehgqJeWlMkU48AXWkLwEkx2L516gRoRzXBQd5oT3784YyeysOP1za9BebKpJ3ZxsQAb1SJ';
}

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

function mapURL(conf: object):object {
  let urls = {};
  // @ts-ignore
  for (let server in conf.servers) {
    // @ts-ignore
    urls[conf['servers'][server]['url']] = conf['servers'][server];
  }
  return urls;
}

export function getDevKey(hostKey: string): {url: string, key: string} {
  let config = yaml.safeLoad(fs.readFileSync(confPath));
  let host = hostKey == '' ? config.default : hostKey;
  host = host.trim();
  let urls = mapURL(config);
  let key = '';
  let url = '';
  try {
    let url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls)
      { // @ts-ignore
        key = config['servers'][urls[url]]['key'];
      }
  } catch (error) {
    if (config['servers'][host] == null)
      throw `Unknown server alias. Please add it to ${confPath}`;
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
  }
  return {url, key};
}

export async function getBrowserPage(): Promise<{browser: puppeteer.Browser, page: puppeteer.Page}> {
  let url:string = process.env.HOST ?? '';
  let cfg = getDevKey(url);
  url = cfg.url;
  if (url.endsWith('/api'))
    url = url.slice(0, -4);
  let key = cfg.key;
  let browser = await puppeteer.launch({
    args: ['--disable-dev-shm-usage', '--disable-features=site-per-process'],
    ignoreHTTPSErrors: true,
  });
  let page = await browser.newPage();
  let token = await getToken(key);
  await page.goto(`${url}/api/info/server`);
  await page.setCookie({name: 'auth', value: token});
  await page.evaluate((token: any) => {
    window.localStorage.setItem('auth', token);
  }, token);
  await page.goto(url);
  try {
    await page.waitForSelector('.grok-view');
  } catch (error) {
    throw 'Can`t load the page';
  }
  return {browser, page};
}