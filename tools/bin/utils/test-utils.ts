import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
const fetch = require('node-fetch');


const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

export async function getToken(url: string, key: string) {
  let response = await fetch(`${url}/users/login/dev/${key}`, {method: 'POST'});
  let json = await response.json();
  if (json.isSuccess == true)
    return json.token;
  else
    throw 'Unable to login to server. Check your dev key';
}

export async function getWebUrl(url: string, token: string) {
  let response = await fetch(`${url}/admin/plugins/admin/settings`, {headers: {Authorization: token}});
  let json = await response.json();
  return json.settings.webRoot;
}

export function getDevKey(hostKey: string): {url: string, key: string} {
  let config = yaml.load(fs.readFileSync(confPath, 'utf8')) as utils.Config;
  let host = hostKey == '' ? config.default : hostKey;
  host = host.trim();
  let urls = utils.mapURL(config);
  let key = '';
  let url = '';
  try {
    let url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls) key = config['servers'][urls[url]]['key'];
  } catch (error) {
    if (config['servers'][host] == null)
      throw `Unknown server alias. Please add it to ${confPath}`;
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
  }
  return {url, key};
}

export async function getBrowserPage(puppeteer: any): Promise<{browser: any, page: any}> {
  let url:string = process.env.HOST ?? '';
  let cfg = getDevKey(url);
  url = cfg.url;

  let key = cfg.key;
  let token = await getToken(url, key);
  url = await getWebUrl(url, token);
  console.log(`Using web root: ${url}`);

  let browser = await puppeteer.launch({
    args: ['--disable-dev-shm-usage', '--disable-features=site-per-process'],
    ignoreHTTPSErrors: true,
  });

  let page = await browser.newPage();
  await page.setDefaultNavigationTimeout(0);
  await page.goto(`${url}/oauth/`);
  await page.setCookie({name: 'auth', value: token});
  await page.evaluate((token: string) => {
    window.localStorage.setItem('auth', token);
  }, token);
  await page.goto(url);
  try {
//    await page.waitForSelector('.grok-preloader', { timeout: 1800000 });
    await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, {timeout: 3600000});
  } catch (error) {
    throw error;
  }
  return {browser, page};
}

export function runWithTimeout(timeout: number, f: () => any): Promise<any> {
  return new Promise(async (resolve, reject) => {
    const timeoutId = setTimeout(() => reject(`Timeout exceeded: ${timeout} ms`), timeout);
    const resolveValue = await f();
    clearTimeout(timeoutId);
    resolve(resolveValue);
  });
}
