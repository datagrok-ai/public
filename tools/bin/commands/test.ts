/* eslint-disable max-len */
import {exec, execFile, spawn, spawnSync} from 'child_process';
import {promisify} from 'util';
import fs from 'fs';
import os from 'os';
import path from 'path';
import readline from 'readline';
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
  'report', 'skip-build', 'skip-publish', 'path', 'record', 'verbose', 'benchmark', 'category', 'test', 'stress-test', 'link', 'tag', 'ci-cd', 'debug', 'no-retry', 'dartium', 'f', 'params'];

const curDir = process.cwd();

/** Expands camelCase to space-separated lowercase: "dataManipulation" → "data manipulation" */
function expandCamelCase(s: string): string {
  return s.replace(/([a-z])([A-Z])/g, '$1 $2').toLowerCase();
}

/** Checks if a segment matches a filter segment (case-insensitive substring, with camelCase expansion) */
function segmentMatches(testSegment: string, filterSegment: string): boolean {
  if (filterSegment === '*')
    return true;
  const t = expandCamelCase(testSegment.trim());
  const f = expandCamelCase(filterSegment.trim());
  return t.includes(f);
}

/**
 * Fluent test name filter. Matches a full test name like "Core: d4 | Viewers | Data Manipulation | Scatter plot".
 *
 * Supports: substring search, `/`-anchored paths, `|`/`/` delimited segments, camelCase expansion, `*` wildcards.
 */
export function matchesFilter(fullName: string, filter: string): boolean {
  // Normalize test name: strip "Core: " prefix, split by " | "
  const normalized = fullName.replace(/^Core:\s*/, '');
  const testSegments = normalized.split(/\s*\|\s*/);

  // Detect anchored mode: ^ prefix
  const anchored = filter.startsWith('^');
  const rawFilter = anchored ? filter.slice(1) : filter;

  // Split filter by / or | (with optional spaces)
  const filterParts = rawFilter.split(/\s*[/|]\s*/).map((s) => s.trim()).filter((s) => s.length > 0);

  // Substring mode: single segment, no delimiters, no anchor
  if (filterParts.length <= 1 && !anchored) {
    const f = expandCamelCase(filterParts[0] || '');
    const full = expandCamelCase(normalized);
    return full.includes(f);
  }

  // Segment mode: find filterParts in order within testSegments
  if (anchored) {
    // Anchored: filter[0] must match testSegments[0], filter[1] must match testSegments[1], etc.
    if (filterParts.length > testSegments.length)
      return false;
    for (let i = 0; i < filterParts.length; i++) {
      if (!segmentMatches(testSegments[i], filterParts[i]))
        return false;
    }
    return true;
  }

  // Unanchored: find filter segments in order (with possible gaps, unless wildcard forces position)
  let ti = 0;
  for (let fi = 0; fi < filterParts.length; fi++) {
    let found = false;
    while (ti < testSegments.length) {
      if (segmentMatches(testSegments[ti], filterParts[fi])) {
        ti++;
        found = true;
        break;
      }
      // If previous filter was wildcard (exact position), don't skip
      if (fi > 0 && filterParts[fi - 1] === '*') {
        // Wildcard consumed exactly one segment, so this filter must match next
        return false;
      }
      ti++;
    }
    if (!found)
      return false;
  }
  return true;
}
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
  if (args.dartium)
    return await testDartium(args);

  if (args.recursive)
    return await testRecursive(process.cwd(), args);

  const config = yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as utils.Config;

  isArgsValid(args);

  // Resolve fluent filter into category/test for normal mode (best-effort)
  const filter = resolveFilter(args);
  if (filter && !args.category && !args.test) {
    const anchored = filter.startsWith('/');
    const raw = anchored ? filter.slice(1) : filter;
    const parts = raw.split(/\s*[/|]\s*/).map((s) => s.trim()).filter((s) => s.length > 0);
    if (parts.length > 1) {
      // Multi-segment: use all but last as category, last as test
      args.category = 'Core: ' + parts.slice(0, -1).join(': ');
      args.test = parts[parts.length - 1];
    }
    else if (parts.length === 1) {
      // Single segment: use as test name filter
      args.test = parts[0];
    }
    color.info(`Filter "${filter}" → category: "${args.category || ''}", test: "${args.test || ''}"`);
  }

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
      } catch (e: any) {
          console.error(e.message || 'Package build/publish failed with no output. Run with --verbose for details.');
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

function resolveFilter(args: TestArgs): string | undefined {
  return args.f || (args['_'].length > 1 ? String(args['_'][1]) : undefined);
}

function isArgsValid(args: TestArgs): boolean {
  const options = Object.keys(args).slice(1);
  if (args['_'].length > 2 || options.length > availableCommandOptions.length || (options.length > 0 &&
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
        keepBrowserOpen: useRetry,
        urlParams: args.params,
      }, browserId, testInvocationTimeout, browserSession);

      // Store browser session for potential reuse
      if (r.browserSession)
        browserSession = r.browserSession;

      if (r.error) {
        color.error(`\nTest execution failed:`);
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
  dartium?: boolean | string,
  f?: string,
  params?: string,
}

interface TestResult {
  name: string;
  version: string;
  tests: string;
  time: string;
  success: boolean;
}

const DARTIUM_WELL_KNOWN_PATHS = [
  'C:\\programs\\dartium-win-ia32-stable-1.24.2.0\\chrome.exe',
  path.join(os.homedir(), 'dartium', 'chrome.exe'),
  path.join(os.homedir(), 'dartium', 'chrome'),
];

const DARTIUM_INACTIVITY_TIMEOUT = 180000; // 3 minutes

function resolveDartiumPath(arg: boolean | string): string {
  if (typeof arg === 'string' && arg !== 'true') {
    if (!fs.existsSync(arg))
      throw new Error(`Dartium not found at: ${arg}`);
    return arg;
  }
  const envPath = process.env.DARTIUM_PATH;
  if (envPath) {
    if (!fs.existsSync(envPath))
      throw new Error(`DARTIUM_PATH set but not found: ${envPath}`);
    return envPath;
  }
  for (const p of DARTIUM_WELL_KNOWN_PATHS) {
    if (fs.existsSync(p))
      return p;
  }
  throw new Error('Dartium not found. Use --dartium=/path/to/chrome.exe or set DARTIUM_PATH');
}

function extractConsoleMessage(line: string): string | null {
  const match = line.match(/INFO:CONSOLE\(\d+\)\] "(.*)", source:/);
  return match ? match[1] : null;
}

interface DartiumTestResult {
  name: string;
  category: string;
  testName: string;
  success: boolean;
  ms: number;
  skipped: boolean;
  error: string;
}

function parseAutotestLine(msg: string): DartiumTestResult | null {
  // AUTOTEST_PASS: d4 | Viewers | Inputs (42ms)
  // AUTOTEST_FAIL: d4 | Viewers | Inputs (42ms) - error message
  // AUTOTEST_SKIP: d4 | Viewers | Inputs - reason
  const passMatch = msg.match(/^AUTOTEST_PASS: (.+?) \((\d+)ms\)$/);
  if (passMatch) {
    const parts = passMatch[1].split(' | ');
    return {
      name: passMatch[1], category: parts.slice(0, -1).join(' | '),
      testName: parts[parts.length - 1], success: true, ms: parseInt(passMatch[2]),
      skipped: false, error: '',
    };
  }
  const failMatch = msg.match(/^AUTOTEST_FAIL: (.+?) \((\d+)ms\) - (.*)$/);
  if (failMatch) {
    const parts = failMatch[1].split(' | ');
    return {
      name: failMatch[1], category: parts.slice(0, -1).join(' | '),
      testName: parts[parts.length - 1], success: false, ms: parseInt(failMatch[2]),
      skipped: false, error: failMatch[3],
    };
  }
  const skipMatch = msg.match(/^AUTOTEST_SKIP: (.+?) - (.*)$/);
  if (skipMatch) {
    const parts = skipMatch[1].split(' | ');
    return {
      name: skipMatch[1], category: parts.slice(0, -1).join(' | '),
      testName: parts[parts.length - 1], success: true, ms: 0,
      skipped: true, error: skipMatch[2],
    };
  }
  return null;
}

async function testDartium(args: TestArgs): Promise<boolean> {
  const config = yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as utils.Config;
  const dartiumPath = resolveDartiumPath(args.dartium!);
  color.info(`Using Dartium: ${dartiumPath}`);

  // Get auth token
  const {url, key} = testUtils.getDevKey(args.host ?? '');
  const token = await testUtils.getToken(url, key);
  const webUrl = await testUtils.getWebUrl(url, token);

  // Build URL
  const filter = resolveFilter(args) || '';
  const urlParams = new URLSearchParams();
  urlParams.set('token', token);
  urlParams.set('tests', filter);
  if (args.params)
    for (const pair of args.params.split('&')) {
      const [k, ...v] = pair.split('=');
      if (k) urlParams.set(k.trim(), v.join('='));
    }
  const testUrl = `${webUrl}/?${urlParams.toString()}`;

  // User-data-dir in temp
  const userDataDir = path.join(os.tmpdir(), 'dartium-grok-test');

  color.info(`Opening: ${webUrl} (tests: ${filter || 'all'})`);

  // Spawn Dartium
  const dartium = spawn(dartiumPath, [
    '--enable-logging=stderr',
    '--no-first-run',
    `--user-data-dir=${userDataDir}`,
    testUrl,
  ], {stdio: ['ignore', 'ignore', 'pipe']});

  const results: DartiumTestResult[] = [];
  let currentCategory = '';
  let passed = 0;
  let failed = 0;
  let skipped = 0;
  let totalExpected = 0;
  let done = false;
  let lastActivityTime = Date.now();

  // Category-level tracking for summary lines
  const categoryResults: Map<string, {passed: number, failed: number, skipped: number}> = new Map();
  const categoryFailures: Map<string, {testName: string, error: string}[]> = new Map();

  const printCategorySummary = (cat: string) => {
    const r = categoryResults.get(cat);
    if (!r) return;
    const skippedSuffix = r.skipped > 0 ? `, \x1b[33m${r.skipped} skipped\x1b[0m` : '';
    if (r.failed > 0) {
      console.log(`\x1b[31m\u274C ${cat}\x1b[31m (\x1b[32m${r.passed} passed${skippedSuffix}\x1b[31m, ${r.failed} failed)\x1b[0m`);
      const failures = categoryFailures.get(cat) || [];
      for (const f of failures)
        console.log(`  \x1b[31m\u274C ${f.testName}\x1b[0m${f.error ? `: ${f.error}` : ''}`);
    }
    else
      console.log(`\x1b[32m\u2714 ${cat} (${r.passed} passed${skippedSuffix})\x1b[0m`);
  };

  return new Promise<boolean>((resolve) => {
    const rl = readline.createInterface({input: dartium.stderr!});

    // Inactivity timeout check
    const inactivityCheck = setInterval(() => {
      if (Date.now() - lastActivityTime > DARTIUM_INACTIVITY_TIMEOUT && !done) {
        process.stdout.write('\r\x1b[K');
        color.error(`\nTest appears stuck (no output for ${DARTIUM_INACTIVITY_TIMEOUT / 1000}s). Killing Dartium.`);
        done = true;
        dartium.kill();
      }
    }, 10000);

    rl.on('line', (line: string) => {
      const msg = extractConsoleMessage(line);
      if (!msg || !msg.startsWith('AUTOTEST_'))
        return;

      lastActivityTime = Date.now();

      if (msg.startsWith('AUTOTEST_START:')) {
        const countMatch = msg.match(/(\d+) test/);
        totalExpected = countMatch ? parseInt(countMatch[1]) : 0;
        color.info(`Running ${totalExpected} test(s)...\n`);
        return;
      }

      if (msg.startsWith('AUTOTEST_RUN:')) {
        const name = msg.replace('AUTOTEST_RUN: ', '');
        process.stdout.write(`\r\x1b[K  \x1b[90m\u25B6 ${name}\x1b[0m`);
        return;
      }

      if (msg.startsWith('AUTOTEST_DONE:')) {
        process.stdout.write('\r\x1b[K');
        // Print last category summary
        if (currentCategory)
          printCategorySummary(currentCategory);

        console.log('');
        const doneMatch = msg.match(/total=(\d+) passed=(\d+) failed=(\d+)/);
        if (doneMatch) {
          const total = parseInt(doneMatch[1]);
          const p = parseInt(doneMatch[2]);
          const f = parseInt(doneMatch[3]);
          const s = total - p - f;
          if (f > 0)
            color.error(`\nResults: ${p} passed, ${f} failed${s > 0 ? `, ${s} skipped` : ''} (${total} total)`);
          else
            color.success(`\nResults: ${p} passed${s > 0 ? `, ${s} skipped` : ''} (${total} total)`);
        }
        done = true;
        dartium.kill();
        return;
      }

      const result = parseAutotestLine(msg);
      if (!result) return;

      process.stdout.write('\r\x1b[K');

      // Category change - print summary of previous category
      if (result.category !== currentCategory) {
        if (currentCategory)
          printCategorySummary(currentCategory);
        currentCategory = result.category;
      }

      // Track category results
      if (!categoryResults.has(result.category))
        categoryResults.set(result.category, {passed: 0, failed: 0, skipped: 0});
      const catR = categoryResults.get(result.category)!;

      if (result.skipped) {
        skipped++;
        catR.skipped++;
      }
      else if (result.success) {
        passed++;
        catR.passed++;
      }
      else {
        failed++;
        catR.failed++;
        if (!categoryFailures.has(result.category))
          categoryFailures.set(result.category, []);
        categoryFailures.get(result.category)!.push({testName: result.testName, error: result.error});
      }

      results.push(result);

      // Verbose: print every test
      if (args.verbose) {
        if (result.skipped)
          console.log(`  \x1b[33m\u25CB ${result.testName} (skipped: ${result.error})\x1b[0m`);
        else if (result.success)
          console.log(`  \x1b[32m\u2714 ${result.testName} (${result.ms}ms)\x1b[0m`);
        else
          console.log(`  \x1b[31m\u274C ${result.testName} (${result.ms}ms) - ${result.error}\x1b[0m`);
      }
    });

    dartium.on('close', () => {
      clearInterval(inactivityCheck);
      rl.close();

      if (!done) {
        // Print last category
        if (currentCategory)
          printCategorySummary(currentCategory);
        console.log('');
        color.warn(`Dartium exited before tests completed (${passed + failed + skipped}/${totalExpected || '?'} tests ran)`);
      }

      // CSV output
      if (args.csv && results.length > 0) {
        const now = new Date().toISOString();
        const csvRows = results.map((r) => ({
          date: now, category: r.category, name: r.testName,
          success: r.success, result: r.skipped ? r.error : (r.success ? 'OK' : r.error),
          ms: r.ms, skipped: r.skipped, error: r.success ? '' : r.error,
        }));
        const csv = Papa.unparse(csvRows);
        fs.writeFileSync(csvReportDir, csv, 'utf8');
        color.info('Saved `test-report.csv`');
      }

      console.log(`\nPassed tests: ${passed}`);
      console.log(`Failed tests: ${failed}`);
      console.log(`Skipped tests: ${skipped}`);

      if (failed > 0)
        testUtils.exitWithCode(1);
      else
        testUtils.exitWithCode(0);
      resolve(failed === 0);
    });

    // Handle ctrl+C
    process.on('SIGINT', () => {
      color.warn('\nInterrupted. Killing Dartium...');
      done = true;
      dartium.kill();
    });
  });
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
  const passedMatch = stdout.match(/Passed (?:amount|tests):\s*(\d+)/);
  const failedMatch = stdout.match(/Failed (?:amount|tests):\s*(\d+)/);
  const skippedMatch = stdout.match(/Skipped (?:amount|tests):\s*(\d+)/);
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
