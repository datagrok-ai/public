import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import {PuppeteerNode} from 'puppeteer';
import type * as DG from 'datagrok-api/dg';
import {PuppeteerScreenRecorder} from 'puppeteer-screen-recorder';
import {spaceToCamelCase} from '../utils/utils';
import puppeteer from 'puppeteer';
import {Browser, Page} from 'puppeteer';
import * as color from '../utils/color-utils';
import Papa from 'papaparse';

const fetch = require('node-fetch');

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const testCollectionTimeout = 600000;

export const defaultLaunchParameters: utils.Indexable = {
  args: [
    '--disable-dev-shm-usage',
    '--disable-features=site-per-process',
    '--window-size=1920,1080',
    '--js-flags=--expose-gc',
  ],
  ignoreHTTPSErrors: true,
  headless: 'new',
  protocolTimeout: 0,
};

export async function getToken(url: string, key: string) {
  const response = await fetch(`${url}/users/login/dev/${key}`, {method: 'POST'});
  const json = await response.json();
  if (json.isSuccess == true)
    return json.token;
  else
    throw new Error('Unable to login to server. Check your dev key');
}

export async function getWebUrl(url: string, token: string) {
  const response = await fetch(`${url}/admin/plugins/admin/settings`, {headers: {Authorization: token}});
  const json = await response.json();
  return json.settings.webRoot;
}

export function getWebUrlFromPage(page: Page): string {
  const url = page.url();
  // Extract the base URL (protocol + host)
  const urlObj = new URL(url);
  return `${urlObj.protocol}//${urlObj.host}`;
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
  return {url, key};
}

export async function getBrowserPage(
  puppeteer: PuppeteerNode, 
  params: any = defaultLaunchParameters): Promise<{ browser: Browser, page: Page }> {
  let url: string = process.env.HOST ?? '';
  const cfg = getDevKey(url);
  url = cfg.url;
  const key = cfg.key;
  const token = await getToken(url, key);
  url = await getWebUrl(url, token);
  console.log(`Using web root: ${url}`);

  const browser = await puppeteer.launch(params);
  if (params.debug) {
    const targets = await browser.targets();
    const devtoolsTarget = targets.find((t) => {
      return t.type() === 'other' && t.url().startsWith('devtools://');
    });
    if (devtoolsTarget) {
      const client = await devtoolsTarget.createCDPSession();
      await client.send('Runtime.enable');
      await client.send('Runtime.evaluate', {
        expression: `
      window.UI.viewManager.showView('network');
      window.UI.dockController.setDockSide('bottom')
    `,
      });
    }
  }

  const page = await browser.newPage();
  await page.setViewport({
    width: 1920,
    height: 1080,
  });
  page.setDefaultNavigationTimeout(0);
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

export async function loadPackage(
    packageDir: string,
    dirName: string,
    hostString: string,
    skipPublish?: boolean,
    skipBuild?: boolean,
    linkPackage?: boolean,
    release?: boolean
): Promise<void> {
  if (skipPublish != true) {
    process.stdout.write(`Building and publishing ${dirName}...`);
    await utils.runScript(`npm install`, packageDir);
    if (linkPackage)
      await utils.runScript(`grok link`, packageDir);
    if (skipBuild != true)
      await utils.runScript(`npm run build`, packageDir);
    await utils.runScript(`grok publish ${hostString}${release ? ' --release' : ''}`, packageDir);
    process.stdout.write(` success!\n`);
  }
}

export async function loadPackages(
    packagesDir: string,
    packagesToLoad?: string,
    host?: string,
    skipPublish?: boolean,
    skipBuild?: boolean,
    linkPackage?: boolean,
    release?: boolean
): Promise<string[]> {
  const packagesToRun = new Map<string, boolean>();
  const hostString = host === undefined ? `` : `${host}`;
  if (packagesToLoad && packagesToLoad !== 'all') {
    const packageNames = packagesToLoad
        .split(' ')
        .map((p) => p.trim())
        .filter((p) => p.length > 0);

    if (skipPublish && skipBuild && !linkPackage)
      return packageNames;

    for (const name of packageNames)
      packagesToRun.set(spaceToCamelCase(name).toLowerCase(), false);
  }

  for (const dirName of fs.readdirSync(packagesDir)) {
    const packageDir = path.join(packagesDir, dirName);
    if (!fs.lstatSync(packageDir).isFile()) {
      let shouldLoad = false;
      try {
        const packageJsonData = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json'), {encoding: 'utf-8'}));
        const packageFriendlyName =
            packagesToRun.get(
                spaceToCamelCase(packageJsonData['friendlyName'] ??
                    packageJsonData['name'].split('/')[1] ?? packageJsonData['name'] ?? '').toLocaleLowerCase() ?? ''
            ) ?? packagesToRun.get(dirName);

        shouldLoad = utils.isPackageDir(packageDir) && (packageFriendlyName !== undefined || packagesToLoad === 'all');
      } catch (e: any) {
        if (utils.isPackageDir(packageDir) &&
            (packagesToRun.get(spaceToCamelCase(dirName).toLocaleLowerCase()) !== undefined || packagesToLoad === 'all'))
          console.log(`Couldn't read package.json  ${dirName}`);
      }
      if (shouldLoad) {
        await loadPackage(packageDir, dirName, hostString, skipPublish, skipBuild, linkPackage, release);
        packagesToRun.set(dirName, true);
      }
    }
  }

  console.log();
  return Array.from(packagesToRun)
      .filter(([_, value]) => value)
      .map(([key, _]) => key);
}

export async function loadTestsList(packages: string[], core: boolean = false, record: boolean = false): Promise<Test[]> {
  const packageTestsData = await timeout(async () => {
    const params = Object.assign({}, defaultLaunchParameters);
    const out = await getBrowserPage(puppeteer, params);
    const browser: Browser = out.browser;
    const page: Page = out.page;

    let recorder = null;
    if (record) {
      const suffix = process.env.BACKUP_SIZE && process.env.WORKER_ID && process.env.TOTAL_WORKERS
          ? `_${process.env.BACKUP_SIZE}_${process.env.WORKER_ID}_${process.env.TOTAL_WORKERS}`
          : '';
      const logsDir = `./load-test-console-output${suffix}.log`;
      const recordDir = `./load-test-record${suffix}.mp4`;

      recorder = new PuppeteerScreenRecorder(page, recorderConfig);
      await recorder.start(recordDir);
      await page.exposeFunction('addLogsToFile', addLogsToFile);

      fs.writeFileSync(logsDir, ``);
      page.on('console', (msg) => {addLogsToFile(logsDir, `CONSOLE LOG ENTRY: ${msg.text()}\n`);});
      page.on('pageerror', (error) => {addLogsToFile(logsDir, `CONSOLE LOG ERROR: ${error.message}\n`);});
      page.on('response', (response) => {
        addLogsToFile(logsDir, `CONSOLE LOG REQUEST: ${response.status()}, ${response.url()}\n`);
      });
    }

    const r = await page.evaluate((packages, coreTests): Promise<LoadedPackageData[] | { failReport: string }> => {
      return new Promise<LoadedPackageData[] | { failReport: string }>((resolve, reject) => {
        const promises: any[] = [];
        try {
          packages.map((packageName: string) => {
            const p = (<any>window).DG.Func.find({package: packageName, name: 'test'})[0]?.package;
            if (p) {
              try {
                promises.push(p.getTests(coreTests).catch((e: any) => {
                  console.error('something else went wrong with collecting package tests');
                  console.error(e?.message);
                })
                  .then((ts: any) => ({packageName: packageName, tests: ts})));
              } catch (e: any) {
                console.error('something went wrong while adding test collection promise');
                console.error(e.message);
              }
            }
          });

        } catch (err) {
          console.error('Error during evaluation in browser context:', err);
          // eslint-disable-next-line prefer-promise-reject-errors
          reject();
        }
        Promise.all(promises)
          .then((results: any[]) => {
            resolve(results);
          }).catch((e: any) => {
            const stack = ((<any>window).DG.Logger.translateStackTrace(e.stack)).then(() => {
              resolve({
                failReport: `${e.message}\n${stack}`,
              });
            });
          });
      });
    }, packages, core);

    await recorder?.stop();
    await browser?.close();

    return r;
  }, testCollectionTimeout);
  const testsList: Test[] = [];

  for (const testPackage of packageTestsData) {
    for (const key in testPackage.tests) {
      if (testPackage.tests.hasOwnProperty(key)) {
        for (const testValue of testPackage.tests[key]?.tests ?? []) {
          testValue.packageName = testPackage.packageName;
          testsList.push(testValue);
        }
      }
    }
  }
  return testsList;
}

export function addLogsToFile(filePath: string, stringToSave: string) {
  fs.appendFileSync(filePath, `${stringToSave}`);
}

export function printBrowsersResult(browserResult: ResultObject, verbose: boolean = false) {
  // Skip detailed summary if modernOutput was used (already printed per-category)
  if (!browserResult.modernOutput) {
    if (verbose) {
      if ((browserResult.passedAmount ?? 0) > 0 && (browserResult.verbosePassed ?? []).length > 0) {
        console.log('Passed: ');
        console.log(browserResult.verbosePassed);
      }
      if ((browserResult.skippedAmount ?? 0) > 0 && (browserResult.verboseSkipped ?? []).length > 0) {
        console.log('Skipped: ');
        console.log(browserResult.verboseSkipped);
      }
    }

    if ((browserResult.failedAmount ?? 0) > 0 && (browserResult.verboseFailed ?? []).length > 0) {
      console.log('Failed: ');
      console.log(browserResult.verboseFailed);
    }
    console.log('Passed amount:  ' + browserResult?.passedAmount);
    console.log('Skipped amount: ' + browserResult?.skippedAmount);
    console.log('Failed amount:  ' + browserResult?.failedAmount);
  }

  if (browserResult.failed) {
    if (browserResult.verboseFailed === 'Package not found')
      color.fail('Tests not found');
    else
      color.fail('Tests failed.');
  } else
    color.success('Tests passed.');
  
}

export function saveCsvResults(stringToSave: string[], csvReportDir: string) {
  const modifiedStrings = stringToSave.map((str, index) => {
    if (index === 0) return str;
    return str.split('\n').slice(1).join('\n');
  });

  fs.writeFileSync(csvReportDir, modifiedStrings.join('\n'), 'utf8');
  color.info('Saved `test-report.csv`\n');
}

async function runTests(testParams: { package: any, params: any }): Promise<any> {
    let failed = false;
    let verbosePassed = '';
    let verboseSkipped = '';
    let verboseFailed = '';
    let countPassed = 0;
    let countSkipped = 0;
    let countFailed = 0;

    try {
        // Check if retry is supported by looking for skipToCategory parameter
        let retrySupported = false;
        try {
            const funcs = (<any>window).DG.Func.find({package: testParams.package, name: 'test'});
            if (funcs && funcs.length > 0) {
                const testFunc = funcs[0];
                retrySupported = testFunc.inputs?.some((input: any) => input.name === 'skipToCategory') === true;
            }
        } catch (e) {
            retrySupported = false;
        }

        const testCallParams: any = {};
        if (testParams.params?.category)
            testCallParams.category = testParams.params.category;
        if (testParams.params?.test)
            testCallParams.test = testParams.params.test;
        // Only pass retry-related params if supported
        if (retrySupported) {
            if (testParams.params?.skipToCategory)
                testCallParams.skipToCategory = testParams.params.skipToCategory;
            if (testParams.params?.skipToTest)
                testCallParams.skipToTest = testParams.params.skipToTest;
            if (testParams.params?.returnOnFail)
                testCallParams.returnOnFail = testParams.params.returnOnFail;
        }

        const df: DG.DataFrame = await (<any>window).grok.functions.call(
            testParams.package + ':test',
            Object.keys(testCallParams).length > 0 ? testCallParams : undefined
        );
        if (!df.getCol('flaking')) {
            const flakingCol = (<any>window).DG.Column.fromType((<any>window).DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount);
            df.columns.add(flakingCol);
        }

        if (!df.getCol('package')) {
            const packageNameCol =
                (<any>window).DG.Column.fromList((<any>window).DG.COLUMN_TYPE.STRING, 'package', Array(df.rowCount).fill(testParams.package));
            df.columns.add(packageNameCol);
        }

        df.columns
            .setOrder([
                'date', 'category', 'name', 'success', 'result', 'ms', 'skipped', 'logs', 'owner', 'package', 'widgetsDifference', 'flaking']);

        df.changeColumnType('result', (<any>window).DG.COLUMN_TYPE.STRING);
        df.changeColumnType('logs', (<any>window).DG.COLUMN_TYPE.STRING);
        if (df.col('widgetsDifference'))
            df.changeColumnType('widgetsDifference', (<any>window).DG.COLUMN_TYPE.STRING);
        // df.changeColumnType('memoryDelta', (<any>window).DG.COLUMN_TYPE.BIG_INT);

        let lastFailedTest: {category: string, test: string} | undefined = undefined;
        for (let row: DG.Row of df.rows) {
            if (row['skipped']) {
                verboseSkipped += `${row['category']}: ${row['name']} (${row['ms']} ms) :  ${row['result']}\n`;
                countSkipped += 1;
            } else if (row['success']) {
                verbosePassed += `${row['category']}: ${row['name']} (${row['ms']} ms) :  ${row['result']}\n`;
                countPassed += 1;
            } else {
                verboseFailed += `${row['category']}: ${row['name']} (${row['ms']} ms) :  ${row['result']}\n`;
                countFailed += 1;
                failed = true;
                if (row['category'] != 'Unhandled exceptions' && row['name'] != 'before' && row['name'] != 'after')
                  lastFailedTest = {category: row['category'], test: row['name']};
            }
        }

        return {
            verbosePassed: verbosePassed,
            verboseSkipped: verboseSkipped,
            verboseFailed: verboseFailed,
            passedAmount: countPassed,
            skippedAmount: countSkipped,
            failedAmount: countFailed,
            csv: df.toCsv(),
            failed: failed,
            lastFailedTest: lastFailedTest,
            retrySupported: retrySupported,
            // df: resultDF?.toJson()
        };
    } catch (e) {
        console.log(`DEBUG: runTests: IN CATCH: ERROR: ${e}`);
        return {failed: true, retrySupported: false, error: `${e}, ${await (<any>window).DG.Logger.translateStackTrace((e as any).stack)}`}
    }

    /*
        for (const testParam of testsParams) {
          lastTest = testParam;
          console.log(testParam);
          const df: DG.DataFrame = await (<any>window).grok.functions.call(testParam.package + ':test', testParam.params);

          if (!df.getCol('flaking')) {
              const flakingCol = (<any>window).DG.Column.fromType((<any>window).DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount);
              df.columns.add(flakingCol);
          }

          if (!df.getCol('package')) {
            const packageNameCol =
                    (<any>window).DG.Column.fromList((<any>window).DG.COLUMN_TYPE.STRING, 'package', Array(df.rowCount).fill(testParam.package));
            df.columns.add(packageNameCol);
          }

          df.columns
            .setOrder([
              'date', 'category', 'name', 'success', 'result', 'ms', 'skipped', 'logs', 'owner', 'package', 'widgetsDifference', 'flaking']);
          // addColumn('flaking', (<any>window).DG.Column.fromType((<any>window).DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount), df);
          // addColumn('package', (<any>window).DG.Column.fromType((<any>window).DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount), df);
          // if (!df.getCol('flaking')) {
          //   const flakingCol = (<any>window).DG.Column.fromType((<any>window).DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount);
          //   df.columns.add(flakingCol);
          // }
          // if (!df.getCol('package')) {
          //   const packageNameCol =
          //     (<any>window).DG.Column.fromList((<any>window).DG.COLUMN_TYPE.STRING, 'package', Array(df.rowCount).fill(testParam.package));
          //   df.columns.add(packageNameCol);
          // }
          if (df.rowCount === 0) {
            verboseFailed += `Test result : Invocation Fail : ${testParam.params.category}: ${testParam.params.test}\n`;
            countFailed += 1;
            failed = true;
            continue;
          }

          let row = df.rows.get(0);
          console.log(`DEBUG: runTests: IN A LOOP: rows in df ${df.rowCount}`);
          if (df.rowCount > 1) {
            const unhandledErrorRow = df.rows.get(1);
            if (!unhandledErrorRow.get('success')) {
              unhandledErrorRow['category'] = row.get('category');
              unhandledErrorRow['name'] = row.get('name');
              row = unhandledErrorRow;
            }
          }
          const category = row.get('category');
          const testName = row.get('name');
          const time = row.get('ms');
          const result = row.get('result');
          const success = row.get('success');
          const skipped = row.get('skipped');
          row['flaking'] = success && (<any>window).DG.Test.isReproducing;

          df.changeColumnType('result', (<any>window).DG.COLUMN_TYPE.STRING);
          df.changeColumnType('logs', (<any>window).DG.COLUMN_TYPE.STRING);
          if (df.col('widgetsDifference'))
            df.changeColumnType('widgetsDifference', (<any>window).DG.COLUMN_TYPE.STRING);
          // df.changeColumnType('memoryDelta', (<any>window).DG.COLUMN_TYPE.BIG_INT);

          if (resultDF === undefined)
            resultDF = df;
          else {
            console.log(`DEBUG: COLUMN NAMES IN RESULT_DF: ${resultDF.columns.names()}`);
            console.log(`DEBUG: COLUMN TYPES IN RESULT_DF: ${resultDF.columns.names().map((c) => `${c}: ${resultDF.col(c)?.type}`)}`);
            console.log(`DEBUG: COLUMN NAMES IN DF: ${df.columns.names()}`);
            console.log(`DEBUG: COLUMN TYPES IN DF: ${df.columns.names().map((c) => `${c}: ${df.col(c)?.type}`)}`);
            resultDF = resultDF.append(df);
          }

          if (row['skipped']) {
            verboseSkipped += `Test result : Skipped : ${time} : ${category}: ${testName} :  ${result}\n`;
            countSkipped += 1;
          } else if (row['success']) {
            verbosePassed += `Test result : Success : ${time} : ${category}: ${testName} :  ${result}\n`;
            countPassed += 1;
          } else {
            verboseFailed += `Test result : Failed : ${time} : ${category}: ${testName} :  ${result}\n`;
            countFailed += 1;
            failed = true;
          }
          if ((success !== true && skipped !== true) && stopOnFail) {
            break;
          }
        }

        if ((<any>window).DG.Test.isInDebug) {
          console.log('on browser closing debug point');
          debugger;
        }

        res = '';
        console.log(`DEBUG: runTests: AFTER THE LOOP: rows in resultDF ${resultDF?.rowCount}`);
        if (resultDF) {
          const bs = (<any>window).DG.BitSet.create(resultDF.rowCount);
          bs.setAll(true);
          for (let i = 0; i < resultDF.rowCount; i++) {
            if (resultDF.rows.get(i).get('category') === 'Unhandled exceptions')
              bs.set(i, false);

          }
          resultDF = resultDF.clone(bs);
          console.log(`DEBUG: runTests: IN IF CONDITION: rows in resultDF ${resultDF.rowCount}`);
          res = resultDF.toCsv();
          console.log(`DEBUG: runTests: IN IF CONDITION: csv length ${res?.length}`);
        }
      } catch (e) {
        failed = true;
        console.log(`DEBUG: runTests: IN CATCH: ERROR: ${e}`);
        error = lastTest ?
            `category: ${
                lastTest.params.category}, name: ${
                lastTest.params.test}, error: ${e}, ${await (<any>window).DG.Logger.translateStackTrace((e as any).stack)}` :
            `test: null, error: ${e}, ${await (<any>window).DG.Logger.translateStackTrace((e as any).stack)}`;
      }
      console.log(`DEBUG: runTests: BEFORE RETURN: RES: ${res.length}, FAILED: ${failed}, ERROR: ${error}`);
      return {
        failed: failed,
        verbosePassed: verbosePassed,
        verboseSkipped: verboseSkipped,
        verboseFailed: verboseFailed,
        passedAmount: countPassed,
        skippedAmount: countSkipped,
        failedAmount: countFailed,
        error: error,
        csv: res,
        // df: resultDF?.toJson()
      };*/
  }

export async function runBrowser(
  testExecutionData: OrganizedTests,
  browserOptions: BrowserOptions,
  browsersId: number,
  testInvocationTimeout: number = 3600000,
  existingBrowserSession?: { browser: Browser, page: Page, webUrl: string }): Promise<ResultObject & { browserSession?: { browser: Browser, page: Page, webUrl: string } }> {
  const testsToRun = {
    func: runTests.toString(),
    testParams: testExecutionData,
  };
  return await timeout(async () => {
    let browser: Browser;
    let page: Page;
    let webUrl: string;

    if (existingBrowserSession) {
      // Reuse existing browser - just refresh the page
      browser = existingBrowserSession.browser;
      page = existingBrowserSession.page;
      webUrl = existingBrowserSession.webUrl;

      // Remove old console listeners before refresh to avoid stale state references
      page.removeAllListeners('console');

      await page.goto(webUrl);
      try {
        await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, {timeout: 3600000});
      } catch (error) {
        throw error;
      }
    } else {
      // Create new browser
      const params = Object.assign({
        devtools: browserOptions.debug,
      }, defaultLaunchParameters);
      if (browserOptions.gui)
        params['headless'] = false;
      const out = await getBrowserPage(puppeteer, params);
      browser = out.browser;
      page = out.page;
      webUrl = await getWebUrlFromPage(page);
    }
    const recorder = new PuppeteerScreenRecorder(page, recorderConfig);
    const currentBrowserNum = browsersId;
    const logsDir = `./test-console-output-${currentBrowserNum}.log`;
    const recordDir = `./test-record-${currentBrowserNum}.mp4`;

    if (browserOptions.record && !existingBrowserSession) {
      // Only set up recording on initial browser creation, not on retry
      await recorder.start(recordDir);
      await page.exposeFunction('addLogsToFile', addLogsToFile);

      fs.writeFileSync(logsDir, ``);
      page.on('console', (msg) => {addLogsToFile(logsDir, `CONSOLE LOG ENTRY: ${msg.text()}\n`);});
      page.on('pageerror', (error) => {addLogsToFile(logsDir, `CONSOLE LOG ERROR: ${error.message}\n`);});
      page.on('response', (response) => {
        addLogsToFile(logsDir, `CONSOLE LOG REQUEST: ${response.status()}, ${response.url()}\n`);
      });
    }

    // State tracking for test output formatting
    const categoryResults: Map<string, {passed: number, failed: number}> = new Map();
    let currentCategory: string | null = null;
    let currentTestName: string | null = null;
    // Store failed tests with their errors to print after category header
    let pendingFailures: {testName: string, error: string | null}[] = [];
    let pendingBeforeFailure: {error: string | null} | null = null;
    let pendingAfterFailure: {error: string | null} | null = null;
    let modernOutput = false;

    const printCategorySummary = (category: string) => {
      const results = categoryResults.get(category);
      if (!results) return;
      const formattedCategory = category.replace(/: /g, ', ');
      const passedCount = results.passed;
      if (results.failed > 0 || pendingBeforeFailure || pendingAfterFailure) {
        console.log(`\x1b[31m❌ ${formattedCategory} (${passedCount} passed)\x1b[0m`);
        // Print before() failure first
        if (pendingBeforeFailure) {
          console.log(`  \x1b[31m❌ before\x1b[0m`);
          if (pendingBeforeFailure.error)
            console.log(`    \x1b[31m${pendingBeforeFailure.error}\x1b[0m`);
        }
        // Print after() failure
        if (pendingAfterFailure) {
          console.log(`  \x1b[31m❌ after\x1b[0m`);
          if (pendingAfterFailure.error)
            console.log(`    \x1b[31m${pendingAfterFailure.error}\x1b[0m`);
        }
        // Print test failures
        for (const failure of pendingFailures) {
          console.log(`  \x1b[31m❌ ${failure.testName}\x1b[0m`);
          if (failure.error)
            console.log(`    \x1b[31m${failure.error}\x1b[0m`);
        }
      } else {
        console.log(`\x1b[32m✅ ${formattedCategory} (${passedCount} passed)\x1b[0m`);
      }
      pendingFailures = [];
      pendingBeforeFailure = null;
      pendingAfterFailure = null;
    };

    // Subscribe to page console events for modern output formatting
    // On retry, old listeners were removed so we need to re-attach
    page.on('console', (msg) => {
      const text = msg.text();
      if (!text.startsWith('Package testing: '))
        return;
      modernOutput = true;
      // Extract tokens in {{...}} format
      const tokens: string[] = [];
      const tokenRegex = /\{\{([^}]+)\}\}/g;
      let match;
      while ((match = tokenRegex.exec(text)) !== null)
        tokens.push(match[1]);

      // Category start: "Package testing: Started {{Category}}"
      if (text.includes('Started') && tokens.length === 1) {
        // Print summary of previous category if exists
        if (currentCategory && categoryResults.has(currentCategory))
          printCategorySummary(currentCategory);
        currentCategory = tokens[0];
        categoryResults.set(currentCategory, {passed: 0, failed: 0});
        pendingFailures = [];
        pendingBeforeFailure = null;
        pendingAfterFailure = null;
      } else if (text.includes('Started') && tokens.length === 2) {
        // Test start: "Package testing: Started {{Category}} {{TestName}}"
        const category = tokens[0].replace(/: /g, ', ');
        currentTestName = tokens[1];
        process.stdout.write(`${category}: ${currentTestName}...`);
      } else if (text.includes('Finished') && tokens.length === 3) {
        // Test finish: "Package testing: Finished {{Category}} {{TestName}} with {{success/error}} for X ms"
        const category = tokens[0];
        const testName = tokens[1];
        const status = tokens[2];
        const results = categoryResults.get(category);

        // Clear the current test line
        process.stdout.write('\r\x1b[K');

        if (status === 'success') {
          if (results) results.passed++;
        } else {
          if (results) results.failed++;
          pendingFailures.push({testName, error: null});
        }
        currentTestName = null;
      } else if (text.includes('Result for') && tokens.length === 2) {
        // Error result: "Package testing: Result for {{Category}} {{TestName}}: error message"
        const errorMsg = text.split(': ').slice(-1)[0];
        // Attach error to the last pending failure
        if (pendingFailures.length > 0)
          pendingFailures[pendingFailures.length - 1].error = errorMsg;
      } else if (text.includes('Category before()') && text.includes('failed') && tokens.length >= 1) {
        // Category before() failed: "Package testing: Category before() {{Category}} failed"
        process.stdout.write('\r\x1b[K');
        pendingBeforeFailure = {error: null};
      } else if (text.includes('Result for') && text.includes('before:') && tokens.length >= 1) {
        // Before error result: "Package testing: Result for {{Category}} before: error message"
        const errorMsg = text.split('before: ').slice(-1)[0];
        if (pendingBeforeFailure)
          pendingBeforeFailure.error = errorMsg;
      } else if (text.includes('Category after()') && text.includes('failed') && tokens.length >= 1) {
        // Category after() failed: "Package testing: Category after() {{Category}} failed"
        process.stdout.write('\r\x1b[K');
        pendingAfterFailure = {error: null};
      } else if (text.includes('Result for') && text.includes('after:') && tokens.length >= 1) {
        // After error result: "Package testing: Result for {{Category}} after: error message"
        const errorMsg = text.split('after: ').slice(-1)[0];
        if (pendingAfterFailure)
          pendingAfterFailure.error = errorMsg;
      } else if (text.includes('Unhandled Exception')) {
        // Unhandled exception: "Package testing: Unhandled Exception: ..."
        // Clear any pending test line and move to new line
        process.stdout.write('\r\x1b[K\n');
        const errorMsg = text.replace('Package testing: Unhandled Exception: ', '');
        if (errorMsg && errorMsg !== 'null')
          console.log(`\x1b[31m❌ Unhandled Exception: ${errorMsg}\x1b[0m`);
      }
    });

    const printFinalCategorySummary = () => {
      if (currentCategory && categoryResults.has(currentCategory))
        printCategorySummary(currentCategory);
    };

    const testingResults = await page.evaluate((testData, options): Promise<ResultObject> => {

      if (options.benchmark)
        (<any>window).DG.Test.isInBenchmark = true;
      if (options.reproduce)
        (<any>window).DG.Test.isReproducing = true;
      if (options.ciCd)
        (<any>window).DG.Test.isCiCd = true;
      if (options.debug)
        (<any>window).DG.Test.isInDebug = true;

      return new Promise<any>((resolve, reject) => {

        (<any>window).runTests = eval('(' + testData.func + ')');

        (<any>window).runTests(testData.testParams, options.stopOnTimeout).then((results: any) => {
          resolve(results);
        })
          .catch((e: any) => {
            resolve({
              failed: true,
              verbosePassed: '',
              verboseSkipped: '',
              verboseFailed: '',
              error: JSON.stringify(e),
              passedAmount: 0,
              skippedAmount: 0,
              failedAmount: 1,
              csv: '',
              df: undefined,
            });
          });

        // (<any>window).DG.Utils.executeTests(testData.tests, options.stopOnTimeout)
        //   .then((results: any) => {
        //     resolve(results);
        //   })
        //   .catch((e: any) => {
        //     resolve({
        //       failed: true,
        //       verbosePassed: "",
        //       verboseSkipped: "",
        //       verboseFailed: "Tests execution failed",
        //       passedAmount: 0,
        //       skippedAmount: 0,
        //       failedAmount: 1,
        //       csv: "",
        //       df: undefined
        //     })
        //   });
      });
    }, testsToRun, browserOptions);

    // Print the final category summary
    printFinalCategorySummary();

    if (browserOptions.record && !existingBrowserSession)
      await recorder.stop();

    if (modernOutput) {
      testingResults.verbosePassed = '';
      testingResults.verboseSkipped = '';
      testingResults.verboseFailed = '';
    }

    // Only close browser if not keeping it open for retry
    if (!browserOptions.keepBrowserOpen && browser != null)
      await browser.close();

    // Return browser session for potential reuse
    return {
      ...testingResults,
      browserSession: browserOptions.keepBrowserOpen ? {browser, page, webUrl} : undefined,
      modernOutput,
    };
  }, testInvocationTimeout);
}

export async function mergeBrowsersResults(browsersResults: ResultObject[]): Promise<ResultObject> {
  const mergedResult: ResultObject = {
    failed: browsersResults[0].failed,
    verbosePassed: browsersResults[0].verbosePassed,
    verboseSkipped: browsersResults[0].verboseSkipped,
    verboseFailed: browsersResults[0].verboseFailed,
    passedAmount: browsersResults[0].passedAmount,
    skippedAmount: browsersResults[0].skippedAmount,
    failedAmount: browsersResults[0].failedAmount,
    csv: browsersResults[0].csv,
    error: '',
    modernOutput: browsersResults.some((r) => r.modernOutput),
  };

  for (const browsersResult of browsersResults) {
    if (mergedResult.csv === browsersResult.csv)
      continue;

    mergedResult.failed = mergedResult.failed || browsersResult.failed;
    mergedResult.verbosePassed = `${mergedResult.verbosePassed.trim()}\n${browsersResult.verbosePassed.trim()}`;
    mergedResult.verboseFailed = `${mergedResult.verboseFailed.trim()}\n${browsersResult.verboseFailed.trim()}`;
    mergedResult.verboseSkipped = `${mergedResult.verboseSkipped.trim()}\n${browsersResult.verboseSkipped.trim()}`;

    mergedResult.passedAmount += browsersResult.passedAmount;
    mergedResult.failedAmount += browsersResult.failedAmount;
    mergedResult.skippedAmount += browsersResult.skippedAmount;

    const resultToMerdge1 = mergedResult.csv.trim().split('\n');
    const resultToMerdge2 = browsersResult.csv.trim().split('\n');

    const header = resultToMerdge1[0];
    mergedResult.csv = [
      header,
      ...resultToMerdge1.slice(1),
      ...resultToMerdge2.slice(1),
    ].join('\n');
  }
  return mergedResult;
}

export interface BrowserOptions {
  path?: string, catchUnhandled?: boolean, core?: boolean,
  report?: boolean, record?: boolean, verbose?: boolean, benchmark?: boolean, platform?: boolean, category?: string, test?: string,
  stressTest?: boolean, gui?: boolean, stopOnTimeout?: boolean, reproduce?: boolean, ciCd?: boolean, debug?: boolean,
  skipToCategory?: string, skipToTest?: string, keepBrowserOpen?: boolean
}

export type ResultObject = {
  failed: boolean,
  verbosePassed: string,
  verboseSkipped: string,
  verboseFailed: string,
  passedAmount: number,
  skippedAmount: number,
  failedAmount: number,
  error?: string,
  csv: string,
  lastFailedTest?: {category: string, test: string},
  retrySupported?: boolean,
  modernOutput?: boolean
};

export class TestContext {
  catchUnhandled = true;
  report = false;

  constructor(catchUnhandled?: boolean, report?: boolean) {
    if (catchUnhandled !== undefined) this.catchUnhandled = catchUnhandled;
    if (report !== undefined) this.report = report;
  };
}

export type LoadedTest = {
  category: string,
  name: string,
}
export type LoadedPackageData = {
  packageName: string,
  tests: { [key: string]: LoadedTest }
}

export type Test = {
  category: string,
  name: string,
  packageName: string,
  options: {
    tags: string[],
    stressTest: boolean,
    benchmark: boolean
  }
}

export type OrganizedTests = {
  package: string,
  params: {
    category: string,
    test: string,
    skipToCategory?: string,
    skipToTest?: string,
    returnOnFail?: boolean,
    options: {
      catchUnhandled: boolean | undefined,
      report: boolean | undefined
    }
  }
}

export function addColumnToCsv(
    csv: string,
    columnName: string,
    defaultValue: string | number | boolean
): string {
  const results = Papa.parse(csv, {
    header: true,
    skipEmptyLines: true,
  });

  const dataWithDefaultColumn = (results.data as any[]).map((row) => ({
    ...row,
    [columnName]: defaultValue,
  }));

  return Papa.unparse(dataWithDefaultColumn, { header: true });
}

export function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
