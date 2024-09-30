/* eslint-disable max-len */
import { exec } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import puppeteer from 'puppeteer';
import { Browser, Page } from 'puppeteer';
import { PuppeteerScreenRecorder } from 'puppeteer-screen-recorder';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import * as testUtils from '../utils/test-utils';
import { WorkerOptions } from '../utils/test-utils';


const testInvocationTimeout = 3600000;

const availableCommandOptions = ['host', 'package', 'csv', 'gui', 'catchUnhandled', 'platform', 'core',
  'report', 'skip-build', 'skip-publish', 'path', 'record', 'verbose', 'benchmark', 'category', 'test', 'stress-test', 'link', 'tag'];

const curDir = process.cwd();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const consoleLogOutputDir = path.join(curDir, 'test-console-output.log');

let browser: Browser;
let page: Page;
let recorder: PuppeteerScreenRecorder;

export async function test(args: TestArgs): Promise<boolean> {
  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

  isArgsValid(args);
  utils.setHost(args.host, config);

  const packageJsonData = JSON.parse(fs.readFileSync(path.join(curDir, 'package.json'), { encoding: 'utf-8' }));
  let packageName = args.package ? utils.kebabToCamelCase(args.package) : utils.kebabToCamelCase(utils.removeScope(packageJsonData.name));
  let packagesDir = path.basename(curDir) === "packages" ? path.dirname(curDir) : curDir;
  let categoryToCheck: string | undefined = undefined;
  let testToCheck: string | undefined = undefined;

  if (args.category) {
    categoryToCheck = args.category.toString();
    if (args.test)
      testToCheck = args.test.toString();
  }

  console.log('Environment variable `TARGET_PACKAGE` is set to', packageName);

  if (args.platform && packageName !== 'ApiTests')
    color.warn('--platform flag can only be used in the ApiTests package');
  if (args.core && packageName !== 'DevTools')
    color.warn('--core flag can only be used in the DevTools package');

  if (!args.package) {
    await testUtils.loadPackages(packagesDir,
      packageName,
      args.host,
      args['skip-publish'],
      args['skip-build'], args.link);
  }

  process.env.TARGET_PACKAGE = packageName;
  let failed = await runTesting(args, categoryToCheck, testToCheck);

  if (failed)
    testUtils.exitWithCode(1);
  else
    testUtils.exitWithCode(0);
  return true;
}

function isArgsValid(args: TestArgs): boolean {
  const options = Object.keys(args).slice(1);
  if (args['_'].length > 1 || options.length > availableCommandOptions.length || (options.length > 0 &&
    !options.every((op) => availableCommandOptions.includes(op))))
    return false;
  return true;
}

async function runTesting(args: TestArgs, categoryToCheck: string | undefined, testToCheck: string | undefined): Promise<boolean> {
  color.info('Starting tests...');

  const r = await runTests(testInvocationTimeout, {
    path: args.path, verbose: args.verbose, platform: args.platform,
    catchUnhandled: args.catchUnhandled, report: args.report, record: args.record, benchmark: args.benchmark,
    core: args.core, category: categoryToCheck, test: testToCheck, stressTest: args['stress-test'], gui: args.gui
  });

  printResults(r, args.csv, args.verbose);
  //@ts-ignore
  if (browser != null)
    await browser.close();
  return r.failed;
}


function runTests(timeout: number, options: WorkerOptions = {}): Promise<ResultObject> {
  return testUtils.timeout(async () => {
    const params = Object.assign({}, testUtils.defaultLaunchParameters);
    if (options.gui)
      params['headless'] = false;

    const out = await testUtils.getBrowserPage(puppeteer, params);
    browser = out.browser;
    page = out.page;
    recorder = new PuppeteerScreenRecorder(page, testUtils.recorderConfig);

    function addLogsToFile(msg: any) {
      fs.appendFileSync(consoleLogOutputDir, `${msg}`);
    }
    await page.exposeFunction("addLogsToFile", addLogsToFile);
    if (options.record) {
      fs.writeFileSync(consoleLogOutputDir, ``);
      await recorder.start('./test-record.mp4');
      page.on('console', (msg) => { addLogsToFile(`CONSOLE LOG ENTRY: ${msg.text()}\n`); });
      page.on('pageerror', (error) => {
        addLogsToFile(`CONSOLE LOG ERROR: ${error.message}\n`);
      });
      page.on('response', (response) => {
        addLogsToFile(`CONSOLE LOG REQUEST: ${response.status()}, ${response.url()}\n`);
      });
    }

    return await runWorker(options);
  }, timeout);
}

async function runWorker(options: WorkerOptions = {}) {
  const targetPackage: string = process.env.TARGET_PACKAGE ?? '#{PACKAGE_NAMESPACE}';
  console.log(`Testing ${targetPackage} package...\n`);
  const r: ResultObject = await page.evaluate((targetPackage, options, testContext): Promise<ResultObject> => {
    if (options.benchmark)
      (<any>window).DG.Test.isInBenchmark = true;
    return new Promise<ResultObject>((resolve, reject) => {
      const params: {
        category?: string,
        test?: string,
        testContext: testUtils.TestContext,
        skipCore?: boolean,
        verbose?: boolean
        stressTest?: boolean
      } = {
        testContext: testContext,
        category: options.category,
        test: options.test
      };

      if (options.stressTest)
        params.stressTest = options.stressTest;

      if (options.path) {
        const split = options.path.split(' -- ');
        params.category = split[0];
        params.test = split[1];
      }

      if (targetPackage === 'DevTools')
        params.skipCore = options.core ? false : true;

      params.verbose = options.verbose === true;

      (<any>window).grok.functions.call(`${targetPackage}:${options.platform ? 'testPlatform' : 'test'}`, params).then((df: any) => {
        let failed = false;
        let skipReport = '';
        let passReport = '';
        let failReport = '';
        const countReport = { skip: 0, pass: 0 };

        if (df == null) {
          failed = true;
          failReport = `Fail reason: No package tests found${options.path ? ' for path "' + options.path + '"' : ''}`;
          resolve({ failReport, skipReport, passReport, failed, countReport });
          return;
        }

        const csv = df.toCsv();
        const cStatus = df.columns.byName('success');
        const cSkipped = df.columns.byName('skipped');
        const cMessage = df.columns.byName('result');
        const cCat = df.columns.byName('category');
        const cName = df.columns.byName('name');
        const cTime = df.columns.byName('ms');

        for (let i = 0; i < df.rowCount; i++) {
          if (cStatus.get(i)) {
            if (cSkipped.get(i)) {
              skipReport += `Test result : Skipped : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
              countReport.skip += 1;
            } else {
              passReport += `Test result : Success : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
              countReport.pass += 1;
            }
          } else {
            failed = true;
            failReport += `Test result : Failed : ${cTime.get(i)} : ${targetPackage}.${cCat.get(i)}.${cName.get(i)} : ${cMessage.get(i)}\n`;
          }
        }

        if (!options.verbose)
          df.rows.removeWhere((r: any) => r.get('success'));

        resolve({ failReport, skipReport, passReport, failed, csv, countReport });
      }).catch((e: any) => {
        const stack = ((<any>window).DG.Logger.translateStackTrace(e.stack)).then(() => {
          resolve({
            failReport: `${e.message}\n${stack}`,
            skipReport: '',
            passReport: '',
            failed: true,
            csv: '',
            countReport: { skip: 0, pass: 0 }
          });
        });
      });
    });
  }, targetPackage, options, new testUtils.TestContext(options.catchUnhandled, options.report));

  if (options.record) {
    await recorder.stop();
  }
  return r;
}

function printResults(results: ResultObject, csv: boolean = false, verbose: boolean = false) {
  if (results.csv && csv) {
    fs.writeFileSync(path.join(curDir, 'test-report.csv'), results.csv, 'utf8');
    color.info('Saved `test-report.csv`\n');
  }

  if (results.passReport && verbose)
    console.log(results.passReport);
  else
    console.log('Passed tests: ' + results.countReport.pass);

  if (results.skipReport && verbose)
    console.log(results.skipReport);
  else
    console.log('Skipped tests: ' + results.countReport.skip);

  if (results.failed) {
    console.log(results.failReport);
    color.fail('Tests failed.');
  } else {
    color.success('Tests passed.');
  }
}

interface TestArgs {
  _: string[],
  category?: any,
  test?: any,
  path?: string,
  host?: string,
  package?: string,
  csv?: boolean,
  gui?: boolean,
  catchUnhandled?: boolean,
  report?: boolean,
  record?: boolean,
  'skip-build'?: boolean,
  'skip-publish'?: boolean,
  link?: boolean,
  verbose?: boolean,
  benchmark?: boolean,
  platform?: boolean,
  core?: boolean,
  'stress-test'?: boolean,
}

type ResultObject = {
  failReport: string,
  skipReport: string,
  passReport: string,
  failed: boolean,
  csv?: string,
  countReport: { skip: number, pass: number }
};