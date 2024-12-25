
import fs from 'fs';
import os from 'os';
import path from 'path'; 
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import * as Papa from 'papaparse';
import * as testUtils from '../utils/test-utils';
import { getBrowserPage, OrganizedTests as OrganizedTest, timeout, defaultLaunchParameters } from '../utils/test-utils';
import { setAlphabeticalOrder } from '../utils/order-functions'; 
import { PuppeteerNode } from 'puppeteer';
import { PuppeteerScreenRecorder } from 'puppeteer-screen-recorder';
import { spaceToCamelCase } from '../utils/utils';
import puppeteer from 'puppeteer';
import { Browser, Page } from 'puppeteer';

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

export  async function loadNPM(args: LoadNPMArgs) {
    const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;
    utils.setHost(args.host, config);
    const params = Object.assign({}, defaultLaunchParameters);
    const out = await getBrowserPage(puppeteer, params);
    const browser: Browser = out.browser;
    const page: Page = out.page; 

    let testingResults = await page.evaluate((): Promise<void> => { 
      return new Promise<any>((resolve, reject) => {
        (<any>window).grok.functions.call("triggerNpmInstall:test", {}).then((results: any) => {
            resolve(results);
          })
          .catch((e:any)=>{
            resolve(e);
          })
      })
    });
    process.exit(0);
}

interface LoadNPMArgs {
    host?: string;
}