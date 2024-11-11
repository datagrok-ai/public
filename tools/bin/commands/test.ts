/* eslint-disable max-len */
import { exec } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import puppeteer from 'puppeteer';
import { Browser, Page } from 'puppeteer';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import * as testUtils from '../utils/test-utils';
import { WorkerOptions, loadTestsList, runWorker, ResultObject, saveCsvResults, printWorkersResult, mergeWorkersResults, Test, OrganizedTests} from '../utils/test-utils';
import { setAlphabeticalOrder } from '../utils/order-functions';

const testInvocationTimeout = 3600000;

const availableCommandOptions = ['host', 'package', 'csv', 'gui', 'catchUnhandled', 'platform', 'core',
  'report', 'skip-build', 'skip-publish', 'path', 'record', 'verbose', 'benchmark', 'category', 'test', 'stress-test', 'link', 'tag'];

const curDir = process.cwd();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const consoleLogOutputDir = path.join(curDir, 'test-console-output.log');
const csvReportDir = path.join(curDir, 'test-report.csv');

export async function test(args: TestArgs): Promise<boolean> {
  const config = yaml.load(fs.readFileSync(confPath, { encoding: 'utf-8' })) as utils.Config;

  isArgsValid(args);
  utils.setHost(args.host, config);

  let packageJsonData = undefined;
  if (!args.package)
    packageJsonData = JSON.parse(fs.readFileSync(path.join(curDir, 'package.json'), { encoding: 'utf-8' }))
  let packageName = args.package ? utils.kebabToCamelCase(args.package) : utils.kebabToCamelCase(utils.removeScope(packageJsonData.name));
  let packagesDir = path.basename(curDir) === "packages" ? curDir : path.dirname(curDir);
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
  let res = await runTesting(args, categoryToCheck, testToCheck);
  if (args.csv)
    saveCsvResults([res.csv], csvReportDir);
  printWorkersResult(res, args.verbose)
  if (!res)
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

async function runTesting(args: TestArgs, categoryToCheck: string | undefined, testToCheck: string | undefined): Promise<ResultObject> {
  color.info('Loading tests...');
  const testsObj = await loadTestsList([process.env.TARGET_PACKAGE ?? ''], args.core);
  const parsed: Test[][] = (setAlphabeticalOrder(testsObj, 1, 1));
  let organized : OrganizedTests[]= parsed[0].map(testObj => ({
    package: testObj.packageName,
    params: {
      category: testObj.category,
      test: testObj.name,
      options: {
        catchUnhandled: args.catchUnhandled,
        report: args.report
      }
    }
  }));
  let filtered: OrganizedTests[] = []
  if (args.category) {
    for (let element of organized) {
      if (element.params.category === args.category) {
        if (element.params.test === args.test || !args.test)
          filtered.push(element);
      }
    }
    organized = filtered;
  }

  color.info('Starting tests...');
  let testsResults: ResultObject[] = [];
  let r :ResultObject;
  do {
    r = await runWorker(organized, {
      benchmark: args.benchmark ?? false,
      catchUnhandled: args.catchUnhandled ?? false,
      gui: args.gui ?? false,
      record: args.record ?? false,
      report: args.report ?? false,
      verbose: args.verbose ?? false,
      stopOnTimeout: true
    }, 1, testInvocationTimeout);
    testsResults.push(r);
    let testsLeft: OrganizedTests[] = [];
    for(let testData of organized){
      if(!r.csv.includes(`${testData.params.category},${testData.params.test}`))
        testsLeft.push(testData);
    }
    organized = testsLeft;
  }
  while (r.verboseFailed.includes('EXECUTION TIMEOUT'));
  return await mergeWorkersResults(testsResults);
}

interface TestArgs {
  _: string[],
  category?: string,
  test?: string,
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