/* eslint-disable max-len */
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
import { setRandomOrder, setAlphabeticalOrder, setPackageRandomOrder, setPackageAlphabeticalOrder } from '../utils/order-functions';
import { WorkerOptions } from '../utils/test-utils';

enum order {
  random = 0,
  alphabetical = 1,
  packageRandom = 2,
  packageAlphabetical = 3,
}

function getEnumOrder(orderStr: string): order {
  switch (orderStr.toLowerCase().replaceAll("-", "")) {
    case "random":
    case "shuffle":
      return order.random;
      break;
    case "alphabetical":
    case "atoz":
      return order.alphabetical;
      break;
    case "packagerandom":
    case "packageshuffle":
      return order.packageRandom;
      break;
    case "packagealphabetical":
    case "packageatoz":
      return order.packageAlphabetical;
      break;
  }
  return order.random;
}

const curDir = process.cwd();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');

const csvReportDir = path.join(curDir, 'test-report.csv');

const testCollectionTimeout = 100000;
const testInvocationTimeout = 7200000;

const orderingFunctions: Map<order, (tests: any, workersAmount: number, testRepeats: number) => any[][]> = new Map<order, (tests: any, workersAmount: number, testRepeats: number) => any[][]>([
  [order.random, setRandomOrder],
  [order.alphabetical, setAlphabeticalOrder],
  [order.packageRandom, setPackageRandomOrder],
  [order.packageAlphabetical, setPackageAlphabeticalOrder]
]);
let workersStarted: number = 0;

export async function testAll(args: TestArgs): Promise<boolean> {
  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

  utils.setHost(args.host, config);
  let packagesToRun = await testUtils.loadPackages(curDir, args.packages, args.host, args['skip-publish'], args['skip-build']);

  let testsObj = await loadTestsList(packagesToRun, args.core);
  let filteredTests = await filterTests(testsObj, (args.tags ?? "").split(" "), args['stress-test'], args.benchmark);
  let workersOrder = await setWorkersOrder(filteredTests, getEnumOrder(args.order ?? ''), args.workersCount, args.testRepeat);

  let testsResults = await runTests(workersOrder, {
    benchmark: args.benchmark ?? false,
    catchUnhandled: args.catchUnhandled ?? false,
    gui: args.gui ?? false,
    record: args.record ?? false,
    report: args.report ?? false,
    verbose: args.verbose ?? false
  });

  let i = 0;
  for (let result of testsResults) {
    console.log(`\nWorker #${i++} `)
    printWorkersResult(result, args.verbose);
  }

  if (args.csv) {
    saveCsvResults(testsResults.map(result => result.csv));
  }

  return !(testsResults.map((test) => test.failed)).some(failStatus => failStatus === true);
}

async function loadTestsList(packages: string[], core: boolean = false): Promise<Object[]> {
  var packageTestsData = await testUtils.timeout(async () => {
    const params = Object.assign({}, testUtils.defaultLaunchParameters);
    // params['headless'] = false;
    const out = await testUtils.getBrowserPage(puppeteer, params);
    const browser: Browser = out.browser;
    const page: Page = out.page;

    const r = await page.evaluate((packages, coreTests): Promise<any> => {
      return new Promise<any>((resolve, reject) => {
        const promises: any[] = [];
        try {
          packages.map((packageName: string) => {
            const p = (<any>window).DG.Func.find({ package: packageName, name: 'test' })[0]?.package;
            if (p) {
              try {
                promises.push(p.getTests(coreTests).catch((e: any) => {
                  console.error('something else went wrong with collecting package tests')
                  console.error(e?.message)
                })
                  .then((ts: any) => ({ packageName: packageName, tests: ts })))
              } catch (e: any) {
                console.error('something went wrong while adding test collection promise');
                console.error(e.message);
              }
            }
          });

        } catch (err) {
          console.error("Error during evaluation in browser context:", err);
          reject();
        }
        Promise.all(promises)
          .then((results: any[]) => {
            resolve(results);
          }).catch((e: any) => {
            const stack = ((<any>window).DG.Logger.translateStackTrace(e.stack)).then(() => {
              resolve({
                failReport: `${e.message}\n${stack}`
              });
            });
          });
      });
    }, packages, core);

    if (browser != null) {
      await browser.close();
    }

    return r;
  }, testCollectionTimeout);

  let testsList: Object[] = [];

  for (let testPackage of packageTestsData) {
    for (const key in testPackage.tests) {
      if (testPackage.tests.hasOwnProperty(key)) {
        for (let testValue of testPackage.tests[key].tests) {
          testValue.packageName = testPackage.packageName;
          testsList.push(testValue);
        }
      }
    }
  }
  return testsList;
}

async function filterTests(tests: any[], tags: string[], stressTest: boolean = false, benchmark: boolean = false): Promise<any[]> {
  let filteredTests: any[] = [];
  let stressTestValue: boolean = tags.includes("stress-test") || stressTest;
  let benchmarkValue: boolean = benchmark;

  tags = tags.filter(tag => tag !== "stress-test" && tag != "")

  for (let test of tests) {
    if (benchmarkValue && test.options?.benchmark !== true)
      continue;
    if (stressTestValue && test.options?.stressTest !== true && !test.options?.tags?.includes("stressTest"))
      continue;

    const hasMatchingTags = tags.some(tag => test.options?.tags?.includes(tag));

    if (hasMatchingTags || tags?.length === 0)
      filteredTests.push(test);
  }

  return filteredTests;
}

async function setWorkersOrder(tests: any[], invocationOrder: order = 0, countOfWorkers: number = 1, testRepeats: number = 1): Promise<any[][]> {
  let resultOrder: any[][] = [];

  let orderingFunction = orderingFunctions.get(invocationOrder);
  if (orderingFunction !== undefined)
    resultOrder = orderingFunction(tests, countOfWorkers, testRepeats);
  else
    throw new Error("Cannot find ordering function");

  return resultOrder;
}

async function runTests(workersOrder: any[][], workerOptions: WorkerOptions): Promise<ResultObject[]> {
  let workersCommands: any[][] = [];

  for (let workerOrder of workersOrder)
    workersCommands.push(workerOrder.map(testObj => ({
      package: testObj.packageName,
      params: {
        category: testObj.category,
        test: testObj.name,
        options: {
          catchUnhandled: workerOptions.catchUnhandled,
          report: workerOptions.report
        }
      }
    })))

  let workersPromises: Promise<ResultObject>[] = [];

  for (let workerCommands of workersCommands) {
    workersPromises.push(runWorker(workerCommands, workerOptions));
    await workersPromises[workersPromises.length];
  }
  let resultObjects = await Promise.all(workersPromises);
  return resultObjects;
}

async function runWorker(testExecutionData: any[], workerOptions: WorkerOptions): Promise<ResultObject> {
  return await testUtils.timeout(async () => {
    const params = Object.assign({}, testUtils.defaultLaunchParameters);
    if (workerOptions.gui)
      params['headless'] = false;
    const out = await testUtils.getBrowserPage(puppeteer, params);
    const browser: Browser = out.browser;
    const page: Page = out.page;
    const recorder = new PuppeteerScreenRecorder(page, testUtils.recorderConfig);

    const currentWorkerNum = workersStarted++;
    const logsDir = `./test-console-output-${currentWorkerNum}.log`;
    const recordDir = `./test-record-${currentWorkerNum}.mp4`;

    if (workerOptions.record) {
      await recorder.start(recordDir);
      await page.exposeFunction("addLogsToFile", addLogsToFile);

      fs.writeFileSync(logsDir, ``);
      page.on('console', (msg) => { addLogsToFile(logsDir, `CONSOLE LOG ENTRY: ${msg.text()}\n`); });
      page.on('pageerror', (error) => { addLogsToFile(logsDir, `CONSOLE LOG ERROR: ${error.message}\n`); });
      page.on('response', (response) => {
        addLogsToFile(logsDir, `CONSOLE LOG REQUEST: ${response.status()}, ${response.url()}\n`);
      });
    }

    let testingResults = await page.evaluate((testData, options): Promise<any> => {
      if (options.benchmark)
        (<any>window).DG.Test.isInBenchmark = true;

      return new Promise<any>((resolve, reject) => {
        (<any>window).DG.Utils.executeTests(testData)
          .then((results: any) => {
            resolve(results);
          })
          .catch((e: any) => {
            resolve({
              failed: true,
              verbosePassed: "",
              verboseSkipped: "",
              verboseFailed: "Tests execution failed",
              passedAmount: 0,
              skippedAmount: 0,
              failedAmount: 1,
              csv: "",
              df: undefined
            })
          });
      })
    }, testExecutionData, workerOptions);

    if (workerOptions.record) {
      await recorder.stop();
    }

    if (browser != null) {
      await browser.close();
    }
    return testingResults;
  }, testInvocationTimeout);
}

function addLogsToFile(filePath: string, stringToSave: any) {
  fs.appendFileSync(filePath, `${stringToSave}`);
}

function printWorkersResult(workerResult: ResultObject, verbose: boolean = false) {
  if (verbose) {
    if ((workerResult.passedAmount ?? 0) > 0 && (workerResult.verbosePassed ?? []).length > 0) {
      console.log("Passed: ");
      console.log(workerResult.verbosePassed);
    }
    if ((workerResult.skippedAmount ?? 0) > 0 && (workerResult.verboseSkipped ?? []).length > 0) {
      console.log("Skipped: ");
      console.log(workerResult.verboseSkipped);
    }
  }

  if ((workerResult.failedAmount ?? 0) > 0 && (workerResult.verboseFailed ?? []).length > 0) {
    console.log("Failed: ");
    console.log(workerResult.verboseFailed);
  }
  console.log("Passed amount:  " + workerResult?.passedAmount);
  console.log("Skipped amount: " + workerResult?.skippedAmount);
  console.log("Failed amount:  " + workerResult?.failedAmount);

  if (workerResult.failed) {
    color.fail('Tests failed.');
  } else {
    color.success('Tests passed.');
  }
}

function saveCsvResults(stringToSave: string[]) { 
  const modifiedStrings = stringToSave.map((str, index) => {
    if (index === 0) return str;
    return str.split('\n').slice(1).join('\n');
  }); 
  
  fs.writeFileSync(csvReportDir, modifiedStrings.join('\n'), 'utf8'); 
  color.info('Saved `test-report.csv`\n');
} 

type ResultObject = {
  failed: boolean,
  verbosePassed: string,
  verboseSkipped: string,
  verboseFailed: string,
  passedAmount: number,
  skippedAmount: number,
  failedAmount: number,
  csv: string,
  df: any
};

interface TestArgs {
  _: string[];

  host?: string;
  'link-package'?: boolean;
  'skip-build'?: boolean;
  'skip-publish'?: boolean;

  catchUnhandled?: boolean,
  core?: boolean,
  csv?: boolean;
  gui?: boolean;
  record?: boolean;
  report?: boolean;
  verbose?: boolean;

  benchmark?: boolean;
  order?: string;
  packages?: string;
  'stress-test'?: boolean;
  tags?: string;
  testRepeat?: number;
  workersCount?: number;
}
