/* eslint-disable max-len */
import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as testUtils from '../utils/test-utils';
import { setRandomOrder, setAlphabeticalOrder, setPackageRandomOrder, setPackageAlphabeticalOrder } from '../utils/order-functions';
import { WorkerOptions, loadTestsList, saveCsvResults, printWorkersResult, runWorker, ResultObject, Test, OrganizedTests } from '../utils/test-utils';

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

const testInvocationTimeout = 7200000;

const orderingFunctions: Map<order, (tests: Test[], workersAmount: number, testRepeats: number) => Test[][]> = new Map<order, (tests: Test[], workersAmount: number, testRepeats: number) => Test[][]>([
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
  let filteredTests: Test[] = await filterTests(testsObj, (args.tags ?? "").split(" "), args['stress-test'], args.benchmark);
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
    saveCsvResults(testsResults.map(result => result.csv), csvReportDir);
  }

  return !(testsResults.map((test) => test.failed)).some(failStatus => failStatus === true);
}

async function filterTests(tests: Test[], tags: string[], stressTest: boolean = false, benchmark: boolean = false): Promise<Test[]> {
  let filteredTests: Test [] = [];
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

async function setWorkersOrder(tests: Test[], invocationOrder: order = 0, countOfWorkers: number = 1, testRepeats: number = 1): Promise<Test[][]> {
  let resultOrder: Test[][] = [];

  let orderingFunction = orderingFunctions.get(invocationOrder);
  if (orderingFunction !== undefined)
    resultOrder = orderingFunction(tests, countOfWorkers, testRepeats);
  else
    throw new Error("Cannot find ordering function");
  return resultOrder;
}

async function runTests(workersOrder: Test[][], workerOptions: WorkerOptions): Promise<ResultObject[]> {
  let workersCommands: OrganizedTests[][] = [];

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
    workersPromises.push(runWorker(workerCommands, workerOptions, workersStarted++, testInvocationTimeout));
    await workersPromises[workersPromises.length];
  }
  let resultObjects = await Promise.all(workersPromises);
  return resultObjects;
}

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
