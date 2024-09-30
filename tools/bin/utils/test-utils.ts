import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import { PuppeteerNode } from 'puppeteer';
import { spaceToCamelCase } from '../utils/utils';

const fetch = require('node-fetch');

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

export const defaultLaunchParameters: utils.Indexable = {
  args: [
    '--disable-dev-shm-usage',
    '--disable-features=site-per-process',
    '--window-size=1920,1080',
  ],
  ignoreHTTPSErrors: true,
  headless: 'new',
  protocolTimeout: 0,
};

export async function getToken(url: string, key: string) {
  const response = await fetch(`${url}/users/login/dev/${key}`, { method: 'POST' });
  const json = await response.json();
  if (json.isSuccess == true)
    return json.token;
  else
    throw new Error('Unable to login to server. Check your dev key');
}

export async function getWebUrl(url: string, token: string) {
  const response = await fetch(`${url}/admin/plugins/admin/settings`, { headers: { Authorization: token } });
  const json = await response.json();
  return json.settings.webRoot;
}

export function getDevKey(hostKey: string): { url: string, key: string } {
  const config = yaml.load(fs.readFileSync(confPath, 'utf8')) as utils.Config;
  let host = hostKey == '' ? config.default : hostKey;
  host = host.trim();
  const urls = utils.mapURL(config);
  let key = '';
  let url = '';
  try {
    let url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    if (url in urls) key = config['servers'][urls[url]]['key'];
  } catch (error) {
    if (config['servers'][host] == null)
      throw new Error(`Unknown server alias. Please add it to ${confPath}`);
    url = config['servers'][host]['url'];
    key = config['servers'][host]['key'];
  }
  return { url, key };
}

export async function getBrowserPage(puppeteer: PuppeteerNode, params: {} = defaultLaunchParameters): Promise<{ browser: any, page: any }> {
  let url: string = process.env.HOST ?? '';
  const cfg = getDevKey(url);
  url = cfg.url;

  const key = cfg.key;
  const token = await getToken(url, key);
  url = await getWebUrl(url, token);
  console.log(`Using web root: ${url}`);

  const browser = await puppeteer.launch(params);

  const page = await browser.newPage();
  await page.setViewport({
    width: 1920,
    height: 1080,
  });
  page.setDefaultNavigationTimeout(0);
  await page.goto(`${url}/oauth/`);
  await page.setCookie({ name: 'auth', value: token });
  await page.evaluate((token: string) => {
    window.localStorage.setItem('auth', token);
  }, token);
  await page.goto(url);
  try {
    //    await page.waitForSelector('.grok-preloader', { timeout: 1800000 });
    await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, { timeout: 3600000 });
  } catch (error) {
    throw error;
  }
  return { browser, page };
}

export function runWithTimeout(timeout: number, f: () => any): Promise<any> {
  return new Promise(async (resolve, reject) => {
    const timeoutId = setTimeout(() => reject(new Error(`Timeout exceeded: ${timeout} ms`)), timeout);
    try {
      const resolveValue = await f();
      clearTimeout(timeoutId);
      resolve(resolveValue);
    } catch (e) {
      reject(e);
    }
  });
}

export async function timeout(func: () => Promise<any>, testTimeout: number, timeoutReason: string = 'EXECUTION TIMEOUT'): Promise<any> {
  let timeout: any = null;
  const timeoutPromise = new Promise<any>((_, reject) => {
    timeout = setTimeout(() => {
      // eslint-disable-next-line prefer-promise-reject-errors
      reject(timeoutReason);
    }, testTimeout);
  });
  try {
    return await Promise.race([func(), timeoutPromise]);
  } finally {
    if (timeout)
      clearTimeout(timeout);
  }
}

export function exitWithCode(code: number): void {
  console.log(`Exiting with code ${code}`);
  process.exit(code);
}

export class TestContext {
  catchUnhandled = true;
  report = false;

  constructor(catchUnhandled?: boolean, report?: boolean) {
    if (catchUnhandled !== undefined) this.catchUnhandled = catchUnhandled;
    if (report !== undefined) this.report = report;
  };
}

export const recorderConfig = {
  followNewTab: true,
  fps: 25,
  ffmpeg_Path: null,
  videoFrame: {
    width: 1280,
    height: 630,
  },
  videoCrf: 18,
  videoCodec: 'libx264',
  videoPreset: 'ultrafast',
  videoBitrate: 1000,
  autopad: {
    color: 'black',
  },
  // aspectRatio: '16:9',
};

export async function loadPackages(packagesDir: string, packagesToLoad?: string, host?: string, skipPublish?: boolean, skipBuild?: boolean, linkPackage?: boolean): Promise<string[]> {
  let packagesToRun = new Map<string, boolean>();
  let hostString = host === undefined ? `` : `${host}`;
  if (packagesToLoad !== "all") {
    for (let pacakgeName of (packagesToLoad ?? "").split(' ')) {
      packagesToRun.set(spaceToCamelCase(pacakgeName).toLocaleLowerCase(), false);
    }
  }

  for (let dirName of fs.readdirSync(packagesDir)) {
    let packageDir = path.join(packagesDir, dirName);
    if (!fs.lstatSync(packageDir).isFile()) {

      try {
        const packageJsonData = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json'), { encoding: 'utf-8' }));
        const packageFriendlyName = packagesToRun.get(spaceToCamelCase(packageJsonData["friendlyName"] ?? packageJsonData["name"].split("/")[1]).toLocaleLowerCase() ?? "");

        if (utils.isPackageDir(packageDir) && (packageFriendlyName !== undefined || packagesToLoad === "all")) {
          try {
            if (skipPublish != true) {
              await utils.runScript(`npm install`, packageDir);
              if (linkPackage)
                await utils.runScript(`grok link`, packageDir);
              if (skipBuild != true)
                await utils.runScript(`npm run build`, packageDir);
              await utils.runScript(`grok publish ${hostString}`, packageDir);
            }
            packagesToRun.set(dirName, true);
            console.log(`Package published ${dirName}`);
          }
          catch (e: any) {
            console.log(`Package wasn't published ${dirName}`);
          }
        }
      }
      catch (e: any) {
        if (utils.isPackageDir(packageDir) && (packagesToRun.get(spaceToCamelCase(dirName).toLocaleLowerCase()) !== undefined || packagesToLoad === "all"))
          console.log(`Couldn't read package.json  ${dirName}`);
      }
    }
  }
  console.log();
  return Array.from(packagesToRun)
    .filter(([key, value]) => value === true)
    .map(([key]) => key);;
}

export interface WorkerOptions {
  path?: string, catchUnhandled?: boolean, core?: boolean,
  report?: boolean, record?: boolean, verbose?: boolean, benchmark?: boolean, platform?: boolean, category?: string, test?: string,
  stressTest?: boolean, gui?:boolean
}