/* eslint-disable max-len */
import {exec, execFile, spawnSync} from 'child_process';
import {promisify} from 'util';
import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import {discoverPackages, applyFilter, PackageInfo, confirm} from './build';

const execAsync = promisify(exec);
const execFileAsync = promisify(execFile);
import * as Papa from 'papaparse';
import * as testUtils from '../utils/test-utils';
import {BrowserOptions, loadTestsList, runBrowser, ResultObject, saveCsvResults, printBrowsersResult, mergeBrowsersResults, Test, OrganizedTests as OrganizedTest, timeout, addColumnToCsv} from '../utils/test-utils';
import {setAlphabeticalOrder} from '../utils/order-functions';

const testInvocationTimeout = 3600000;

const availableCommandOptions = ['host', 'package', 'csv', 'gui', 'catchUnhandled', 'platform', 'core',
  'report', 'skip-build', 'skip-publish', 'path', 'record', 'verbose', 'benchmark', 'category', 'test', 'stress-test', 'link', 'tag', 'ci-cd', 'debug', 'no-retry'];

const curDir = process.cwd();
const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const consoleLogOutputDir = path.join(curDir, 'test-console-output.log');
const csvReportDir = path.join(curDir, 'test-report.csv');

/**
 * Detects if the current directory is within a Dart library folder (d4, xamgle, ddt, dml)
 * and returns the appropriate test category to use.
 * @returns The category string (e.g., "Core: d4") or undefined if not in a recognized Dart folder
 */
function detectDartLibraryCategory(): string | undefined {
  const normalizedPath = curDir.replace(/\\/g, '/');

  if (normalizedPath.includes('/d4/') || normalizedPath.endsWith('/d4'))
    return 'Core: d4';
  if (normalizedPath.includes('/xamgle/') || normalizedPath.endsWith('/xamgle'))
    return 'Core: xamgle';
  if (normalizedPath.includes('/ddt/') || normalizedPath.endsWith('/ddt'))
    return 'Core: ddt';
  if (normalizedPath.includes('/dml/') || normalizedPath.endsWith('/dml'))
    return 'Core: dml';

  return undefined;
}

/**
 * Traverses up from startDir to find the git repository root (.git directory).
 */
function findGitRoot(startDir: string): string | undefined {
  let current = startDir;
  while (true) {
    if (fs.existsSync(path.join(current, '.git')))
      return current;
    const parent = path.dirname(current);
    if (parent === current)
      return undefined;
    current = parent;
  }
}

export async function test(args: TestArgs): Promise<boolean> {
  if (args.recursive)
    return await testRecursive(process.cwd(), args);

  const config = yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as utils.Config;

  isArgsValid(args);

  // If running from a core Dart library directory, delegate to DevTools package
  if (!args.package) {
    const detectedCategory = detectDartLibraryCategory();
    if (detectedCategory) {
      const category = args.category ?? detectedCategory;
      const gitRoot = findGitRoot(curDir);
      const devToolsDir = gitRoot ? path.join(gitRoot, 'public', 'packages', 'DevTools') : undefined;
      if (!devToolsDir || !fs.existsSync(devToolsDir)) {
        color.error(`Cannot run core tests from this directory: DevTools package not found.`);
        color.error(`Run 'grok test --category="${category}"' from 'public/packages/DevTools' instead.`);
        process.exit(1);
      }
      const hostAlias = args.host ?? config.default;
      color.info(`Detected core library directory. Delegating to DevTools with category: "${category}"${hostAlias ? ` (host: ${hostAlias})` : ''}`);
      const cmdArgs = ['test', `--category=${category}`];
      if (args.host) cmdArgs.push(`--host=${args.host}`);
      if (args.test) cmdArgs.push(`--test=${args.test}`);
      if (args.csv) cmdArgs.push('--csv');
      if (args.gui) cmdArgs.push('--gui');
      if (args.verbose) cmdArgs.push('--verbose');
      if (args.benchmark) cmdArgs.push('--benchmark');
      if (args['stress-test']) cmdArgs.push('--stress-test');
      if (args['skip-build']) cmdArgs.push('--skip-build');
      if (args['skip-publish']) cmdArgs.push('--skip-publish');
      if (args.link) cmdArgs.push('--link');
      if (args.record) cmdArgs.push('--record');
      if (args.catchUnhandled) cmdArgs.push('--catchUnhandled');
      if (args.report) cmdArgs.push('--report');
      if (args.debug) cmdArgs.push('--debug');
      if (args['ci-cd']) cmdArgs.push('--ci-cd');
      if (args['no-retry']) cmdArgs.push('--no-retry');
      if (!args['skip-publish']) {
        const isDevToolsOnServer = await testUtils.isPackageOnServer(args.host ?? '', 'DevTools');
        if (isDevToolsOnServer) {
          cmdArgs.push('--skip-publish');
          if (!args['skip-build']) cmdArgs.push('--skip-build');
        }
      }
      const grokScript = path.resolve(__dirname, '..', 'grok.js');
      const result = spawnSync(process.execPath, [grokScript, ...cmdArgs], {cwd: devToolsDir, stdio: 'inherit'});
      process.exit(result.status ?? 1);
    }
  }

  utils.setHost(args.host, config, true);

  let packageJsonData = undefined;
  if (!args.package)
    packageJsonData = JSON.parse(fs.readFileSync(path.join(curDir, 'package.json'), {encoding: 'utf-8'}));
  const packageName = args.package ? utils.kebabToCamelCase(args.package) : utils.kebabToCamelCase(utils.removeScope(packageJsonData.name));
  const packagesDir = path.basename(curDir) === 'packages' ? curDir : path.dirname(curDir);

  console.log(`HOST: ${process.env.HOST}, TARGET_PACKAGE: ${packageName}`);

  if (args.platform && packageName !== 'ApiTests')
    color.warn('--platform flag can only be used in the ApiTests package');
  // if (args.core && packageName !== 'DevTools')
  //   color.warn('--core flag can only be used in the DevTools package');


  if (!args.package) {
      try {
          await testUtils.loadPackages(packagesDir,
              packageName,
              args.host,
              args['skip-publish'],
              args['skip-build'], args.link);
      } catch (e) {
          console.error('\n');
          // @ts-ignore
          console.error(e.message);
          process.exit(1);
      }
  }
  process.env.TARGET_PACKAGE = packageName;
  const res = await runTesting(args);
  if (args.csv) {
      res.csv = addColumnToCsv(res.csv, 'stress_test', args['stress-test'] ?? false);
      res.csv = addColumnToCsv(res.csv, 'benchmark', args.benchmark ?? false);
      saveCsvResults([res.csv], csvReportDir);
  }
  printBrowsersResult(res, args.verbose);
  if (res.failed) {
    if (res.verboseFailed === 'Package not found')
      testUtils.exitWithCode(0);
    testUtils.exitWithCode(1);
  } else
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

// Retry state - persists across test runs in the session
const retriedTests: Set<string> = new Set();
let totalRetries = 0;
const MAX_RETRIES_PER_SESSION = 10;
let retryEnabled = true;

async function runTesting(args: TestArgs): Promise<ResultObject> {

  retryEnabled = args['retry'] ?? true;
  if (args.test || args.category)
      retryEnabled = false;
  let organized: OrganizedTest = {
    package: process.env.TARGET_PACKAGE ?? '',
    params: {
      category: args.category ?? '',
      test: args.test ?? '',
      options: {
        catchUnhandled: args.catchUnhandled,
        report: args.report,
      },
    },
  };

  color.info('Starting tests...');
  const testsResults: ResultObject[] = [];
  let r: ResultObject & { browserSession?: { browser: any, page: any, webUrl: string } };
  let browserId = 1;
  let retrySupported: boolean | undefined = undefined; // Will be set after first run
  let browserSession: { browser: any, page: any, webUrl: string } | undefined = undefined;

  await timeout(async () => {
    let shouldRetry = true;
    let currentSkipToCategory: string | undefined = undefined;
    let currentSkipToTest: string | undefined = undefined;

    let oneLastTry = false;
    while (shouldRetry || oneLastTry) {
      shouldRetry = false;
      oneLastTry = false;
      // On first run, assume retry is supported; after first run, use actual value
      const useRetry = retryEnabled && (retrySupported === undefined || retrySupported);

      const testParams: OrganizedTest = {
        ...organized,
        params: {
          ...organized.params,
          skipToCategory: currentSkipToCategory,
          skipToTest: currentSkipToTest,
          returnOnFail: useRetry,
        },
      };

      r = await runBrowser(testParams, {
        benchmark: args.benchmark ?? false,
        stressTest: args['stress-test'] ?? false,
        catchUnhandled: args.catchUnhandled ?? false,
        gui: args.gui ?? false,
        record: args.record ?? false,
        report: args.report ?? false,
        verbose: args.verbose ?? false,
        ciCd: args['ci-cd'] ?? false,
        stopOnTimeout: false,
        debug: args['debug'] ?? false,
        skipToCategory: currentSkipToCategory,
        skipToTest: currentSkipToTest,
        keepBrowserOpen: useRetry, // Keep browser open if retry is enabled
      }, browserId, testInvocationTimeout, browserSession);

      // Store browser session for potential reuse
      if (r.browserSession)
        browserSession = r.browserSession;

      if (r.error) {
        console.log(`\nexecution error:`);
        console.log(r.error);
        // Close browser on error
        if (browserSession?.browser)
          await browserSession.browser.close();
        // Add the failed result before returning
        testsResults.push(r);
        return;
      }
      // Check retry support from first run result
      if (useRetry && r.failed && r.lastFailedTest && retrySupported === undefined) {
        retrySupported = r.retrySupported === true;
        if (!retrySupported)
          color.warn('Retry not supported: test() function does not have skipToCategory parameter');
      }

      // Check if we should retry on failure
      if (useRetry && r.failed && r.lastFailedTest && retryEnabled && retrySupported) {
        const testKey = `${r.lastFailedTest.category}::${r.lastFailedTest.test}`;

        if (totalRetries >= MAX_RETRIES_PER_SESSION) {
          color.warn(`Maximum retries (${MAX_RETRIES_PER_SESSION}) reached. Disabling retries.`);
          retryEnabled = false;
          oneLastTry = true;
        } else {
          // Retry this test
          retriedTests.add(testKey);
          totalRetries++;
          console.log('Refreshing page for retry...');
          color.info(`Retrying from "${r.lastFailedTest.category}: ${r.lastFailedTest.test}" (retry ${totalRetries}/${MAX_RETRIES_PER_SESSION})...`);

          // Store results from previous run
          testsResults.push(r);

          // Set skip params for next run
          currentSkipToCategory = r.lastFailedTest.category;
          currentSkipToTest = r.lastFailedTest.test;
          shouldRetry = true;
          browserId++;
        }
      }

      if (!shouldRetry) {
        testsResults.push(r);
      }

    }

    // Close browser after all retries are done
    if (browserSession?.browser)
      await browserSession.browser.close();
  }, testInvocationTimeout);

  // Handle empty results (shouldn't happen but safety check)
  if (testsResults.length === 0) {
    return {
      failed: true,
      verbosePassed: '',
      verboseSkipped: '',
      verboseFailed: 'No test results collected',
      passedAmount: 0,
      skippedAmount: 0,
      failedAmount: 0,
      csv: '',
      error: 'No test results collected',
    };
  }

  return await mergeBrowsersResults(testsResults);
}

async function reproducedTest(args: TestArgs, testsToReproduce: OrganizedTest[]): Promise<Map<OrganizedTest, ResultObject>> {
  const res: Map<OrganizedTest, ResultObject> = new Map<OrganizedTest, ResultObject>();
  for (const test of testsToReproduce) {
    const r = await runBrowser([test], {
      benchmark: args.benchmark ?? false,
      catchUnhandled: false,
      gui: false,
      record: false,
      report: false,
      verbose: false,
      stopOnTimeout: true,
      reproduce: true,
      ciCd: args['ci-cd'] ?? false,
    }, 0, testInvocationTimeout);
    if (test.params.category && test.params.test)
      res.set(test, r);
  }
  return res;
}

async function updateResultsByReproduced(curentResult: ResultObject, reproducedResult: ResultObject, testsParams: OrganizedTest): Promise<ResultObject> {
  const table2Dict: Record<string, Record<string, string>> = {};
  const table1 = readCSVResultData(curentResult.csv);
  const table2 = readCSVResultData(reproducedResult.csv);
  const flakingMap: Record<string, string> = {};
  table2.rows.forEach((row) => {
    const key = `${row['category']},${row['name']}`;
    flakingMap[key] = row['flaking'];
  });

  table1.rows.forEach((row) => {
    const key = `${row['category']},${row['name']}`;
    if (key in flakingMap)
      row['flaking'] = flakingMap[key];

  });

  curentResult.csv = Papa.unparse(table1.rows, {columns: table1.headers}); ;
  curentResult.verboseFailed = curentResult.verboseFailed.replaceAll(`${testsParams.params.category}: ${testsParams.params.test} :  Error:`, `${testsParams.params.category}: ${testsParams.params.test} : Flaking Error:`);
  return curentResult;
}

function readCSVResultData(data: string): { headers: string[], rows: Record<string, string>[] } {
  const parsed = Papa.parse(data, {
    header: true,
    skipEmptyLines: true,
  });
  if (parsed.errors.length > 0)
    throw new Error(`Error parsing CSV file: ${parsed.errors[0].message}`);

  return {headers: parsed.meta.fields || [], rows: parsed.data as Record<string, string>[]};
}

interface TestArgs {
  _: string[],
  category?: string,
  test?: string,
  host?: string,
  package?: string,
  csv?: boolean,
  gui?: boolean,
  'debug'?: boolean,
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
  'ci-cd'?: boolean,
  'no-retry'?: boolean,
  recursive?: boolean,
  filter?: string,
  parallel?: number,
}

interface TestResult {
  name: string;
  version: string;
  tests: string;
  time: string;
  success: boolean;
}

async function testRecursive(baseDir: string, args: TestArgs): Promise<boolean> {
  const packages = discoverPackages(baseDir);
  if (packages.length === 0) {
    color.warn('No packages found in the current directory');
    return false;
  }

  const filtered = args.filter ? applyFilter(packages, args.filter) : packages;
  if (filtered.length === 0) {
    color.warn('No packages match the filter');
    return false;
  }

  console.log(`Found ${filtered.length} package(s): ${filtered.map((p) => p.friendlyName).join(', ')}`);

  const confirmed = await confirm(`\nTest ${filtered.length} package(s)?`);
  if (!confirmed) {
    console.log('Aborted.');
    return false;
  }

  const maxParallel = args.parallel || 4;
  const results = await testParallel(filtered, args, maxParallel);
  return results.every((r) => r.success);
}

function buildTestArgs(args: TestArgs): string[] {
  const grokScript = path.resolve(__dirname, '..', 'grok.js');
  const parts = [grokScript, 'test'];
  if (args.host) parts.push(`--host=${args.host}`);
  if (args['skip-build']) parts.push('--skip-build');
  if (args['skip-publish']) parts.push('--skip-publish');
  if (args.link) parts.push('--link');
  if (args.verbose) parts.push('--verbose');
  if (args.csv) parts.push('--csv');
  if (args.catchUnhandled) parts.push('--catchUnhandled');
  if (args.report) parts.push('--report');
  if (args.benchmark) parts.push('--benchmark');
  if (args['stress-test']) parts.push('--stress-test');
  if (args['ci-cd']) parts.push('--ci-cd');
  if (args['no-retry']) parts.push('--no-retry');
  return parts;
}

function parseTestOutput(stdout: string): {passed: number, failed: number, skipped: number} | null {
  const passedMatch = stdout.match(/Passed amount:\s*(\d+)/);
  const failedMatch = stdout.match(/Failed amount:\s*(\d+)/);
  const skippedMatch = stdout.match(/Skipped amount:\s*(\d+)/);
  if (!passedMatch && !failedMatch)
    return null;
  return {
    passed: passedMatch ? parseInt(passedMatch[1]) : 0,
    failed: failedMatch ? parseInt(failedMatch[1]) : 0,
    skipped: skippedMatch ? parseInt(skippedMatch[1]) : 0,
  };
}

async function testParallel(packages: PackageInfo[], args: TestArgs, maxParallel: number): Promise<TestResult[]> {
  const results: TestResult[] = [];
  const cmdArgs = buildTestArgs(args);

  const headers = ['Plugin', 'Version', 'Tests', 'Time'];
  const widths = [
    Math.max(headers[0].length, ...packages.map((p) => p.friendlyName.length)),
    Math.max(headers[1].length, ...packages.map((p) => p.version.length)),
    Math.max(headers[2].length, 18),
    Math.max(headers[3].length, 8),
  ];
  const pad = (s: string, w: number) => s + ' '.repeat(Math.max(0, w - s.length));

  console.log(`\nTesting with ${maxParallel} parallel job(s)...`);
  console.log(headers.map((h, i) => pad(h, widths[i])).join(' | '));
  console.log(widths.map((w) => '-'.repeat(w)).join('-+-'));

  const testOne = async (pkg: PackageInfo): Promise<TestResult> => {
    const start = Date.now();
    let success = true;
    let tests = '';
    let time = '';

    try {
      const {stdout} = await execFileAsync(process.execPath, cmdArgs, {
        cwd: pkg.dir,
        maxBuffer: 50 * 1024 * 1024,
        timeout: testInvocationTimeout,
      });
      const elapsed = (Date.now() - start) / 1000;
      time = `${elapsed.toFixed(1)}s`;
      const counts = parseTestOutput(stdout);
      if (counts) {
        const total = counts.passed + counts.failed + counts.skipped;
        tests = `${counts.passed} of ${total} passed`;
        success = counts.failed === 0;
      }
      else {
        tests = 'Done (no counts)';
      }
    }
    catch (error: any) {
      success = false;
      const elapsed = (Date.now() - start) / 1000;
      time = `${elapsed.toFixed(1)}s`;
      const stdout = error.stdout || '';
      const counts = parseTestOutput(stdout);
      if (counts) {
        const total = counts.passed + counts.failed + counts.skipped;
        tests = `${counts.passed} of ${total} passed`;
      }
      else {
        const raw = (error.stderr || error.stdout || error.message || 'Unknown error').trim();
        tests = raw.replace(/\r?\n/g, ' ').replace(/\s+/g, ' ').substring(0, 30);
      }
    }

    const result: TestResult = {name: pkg.friendlyName, version: pkg.version, tests, time, success};
    results.push(result);

    const cells = [result.name, result.version, result.tests, result.time];
    const line = cells.map((cell, j) => pad(cell, widths[j])).join(' | ');
    if (success)
      color.info(line);
    else
      color.error(line);

    return result;
  };

  let idx = 0;
  const next = async (): Promise<void> => {
    while (idx < packages.length) {
      const pkg = packages[idx++];
      await testOne(pkg);
    }
  };
  const workers = Array.from({length: Math.min(maxParallel, packages.length)}, () => next());
  await Promise.all(workers);

  const succeeded = results.filter((r) => r.success).length;
  const failed = results.length - succeeded;
  console.log('');
  if (failed === 0)
    color.success(`All ${results.length} package(s) passed`);
  else
    color.warn(`${succeeded} passed, ${failed} failed`);

  return results;
}
