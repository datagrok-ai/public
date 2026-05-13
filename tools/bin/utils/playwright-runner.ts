import {spawn} from 'child_process';
import fs from 'fs';
import path from 'path';
import Papa from 'papaparse';
import * as color from './color-utils';
import * as testUtils from './test-utils';
import {ResultObject} from './test-utils';

export function hasPlaywrightTests(pkgDir: string): string | null {
  const pkgJsonPath = path.join(pkgDir, 'package.json');
  if (!fs.existsSync(pkgJsonPath))
    return null;
  let pkgJson: any;
  try {
    pkgJson = JSON.parse(fs.readFileSync(pkgJsonPath, 'utf-8'));
  } catch {
    return null;
  }
  const rel = pkgJson.playwrightTests;
  if (typeof rel !== 'string' || rel.length === 0)
    return null;
  const abs = path.resolve(pkgDir, rel);
  if (!fs.existsSync(abs))
    return null;
  return abs;
}

interface PlaywrightArgs {
  test?: string;
  category?: string;
  gui?: boolean;
  verbose?: boolean;
  'no-retry'?: boolean;
  package?: string;
}

interface PwTestResult {
  status: 'passed' | 'failed' | 'timedOut' | 'skipped' | 'interrupted';
  duration: number;
  errors?: Array<{message?: string; stack?: string}>;
  stdout?: Array<{text?: string}>;
  stderr?: Array<{text?: string}>;
}

interface PwTest {
  results?: PwTestResult[];
}

interface PwSpec {
  title: string;
  file: string;
  tests?: Array<PwTest & {projectName?: string}>;
}

interface PwSuite {
  title?: string;
  file?: string;
  specs?: PwSpec[];
  suites?: PwSuite[];
}

interface PwReport {
  config?: {rootDir?: string};
  suites?: PwSuite[];
}

interface FlatRow {
  date: string;
  category: string;
  name: string;
  success: boolean;
  result: string;
  ms: number;
  skipped: boolean;
  logs: string;
  owner: string;
  package: string;
  widgetsDifference: string;
  flaking: boolean;
}

function flattenSuites(
  suites: PwSuite[] | undefined,
  testDir: string,
  pkgName: string,
  owner: string,
  verbose: boolean,
  rows: FlatRow[],
): void {
  if (!suites)
    return;
  const isoDate = new Date().toISOString();
  for (var suite of suites) {
    if (suite.specs) {
      for (var spec of suite.specs) {
        const specFile = spec.file || suite.file || '';
        const absSpec = path.isAbsolute(specFile) ? specFile : path.resolve(testDir, specFile);
        const category = path.relative(testDir, path.dirname(absSpec)).replace(/\\/g, '/');
        for (var t of (spec.tests || [])) {
          const result = (t.results && t.results[t.results.length - 1]) || undefined;
          if (!result)
            continue;
          const status = result.status;
          const skipped = status === 'skipped';
          const success = status === 'passed';
          var errMsg = '';
          if (!success && !skipped && result.errors && result.errors.length > 0)
            errMsg = result.errors.map((e) => e.message || e.stack || '').filter((s) => s.length > 0).join('\n');
          var logs = '';
          if (verbose) {
            const out = (result.stdout || []).map((s) => s.text || '').join('');
            const err = (result.stderr || []).map((s) => s.text || '').join('');
            logs = [out, err].filter((s) => s.length > 0).join('\n');
          }
          rows.push({
            date: isoDate,
            category: category,
            name: spec.title,
            success: success,
            result: errMsg,
            ms: Math.round(result.duration || 0),
            skipped: skipped,
            logs: logs,
            owner: owner,
            package: pkgName,
            widgetsDifference: '',
            flaking: false,
          });
        }
      }
    }
    if (suite.suites)
      flattenSuites(suite.suites, testDir, pkgName, owner, verbose, rows);
  }
}

function rowsToCsv(rows: FlatRow[]): string {
  // Column order must match what runTesting() actually serializes (see
  // test-utils.ts:485 — `setOrder` is a no-op for the CSV writer, so the order
  // is the natural order Dart's package-test emits). mergeBrowsersResults uses
  // the first CSV's header, so any drift here misaligns Playwright rows.
  const header = ['date', 'category', 'name', 'success', 'ms', 'skipped',
    'owner', 'package', 'flaking', 'result', 'logs', 'widgetsDifference'];
  return Papa.unparse({fields: header, data: rows.map((r) => header.map((h) => (r as any)[h]))});
}

export async function runPlaywrightTests(
  pkgDir: string,
  testDir: string,
  args: PlaywrightArgs,
  hostKey: string,
): Promise<ResultObject> {
  const empty: ResultObject = {
    failed: false, verbosePassed: '', verboseSkipped: '', verboseFailed: '',
    passedAmount: 0, skippedAmount: 0, failedAmount: 0, csv: '',
  };

  let url: string; let key: string;
  try {
    ({url, key} = testUtils.getDevKey(hostKey));
  } catch (e: any) {
    color.error(`Playwright: cannot resolve host '${hostKey}': ${e.message || e}`);
    return {...empty, failed: true, failedAmount: 1, verboseFailed: `Playwright: ${e.message || e}\n`};
  }

  let token: string;
  try {
    token = await testUtils.getToken(url, key);
  } catch (e: any) {
    color.error(`Playwright: cannot exchange dev key for token: ${e.message || e}`);
    return {...empty, failed: true, failedAmount: 1, verboseFailed: `Playwright: ${e.message || e}\n`};
  }

  let webUrl: string;
  try {
    webUrl = await testUtils.getWebUrl(url, token);
    if (webUrl.endsWith('/'))
      webUrl = webUrl.slice(0, -1);
  } catch {
    webUrl = url.replace(/\/api\/?$/, '');
  }

  let token2 = '';
  if (process.env.DATAGROK_DEV_KEY_2 && process.env.DATAGROK_DEV_KEY_2.length > 0) {
    try {
      token2 = await testUtils.getToken(url, process.env.DATAGROK_DEV_KEY_2);
    } catch (e: any) {
      color.warn(`Playwright: DATAGROK_DEV_KEY_2 set but failed to exchange for token: ${e.message || e}`);
    }
  }

  const configPath = path.join(testDir, 'playwright.config.ts');
  if (!fs.existsSync(configPath)) {
    color.error(`Playwright: ${configPath} not found.`);
    return {...empty, failed: true, failedAmount: 1, verboseFailed: 'Playwright: missing playwright.config.ts\n'};
  }

  const reportFile = path.join(pkgDir, 'test-playwright-report.json');
  if (fs.existsSync(reportFile))
    fs.unlinkSync(reportFile);

  const cliArgs = ['--no-install', 'playwright', 'test', `--config=${configPath}`];
  if (!args.gui)
    cliArgs.push(`--reporter=json`);
  else
    cliArgs.push('--headed');
  if (args.test)
    cliArgs.push(`--grep=${args.test}`);
  if (args['no-retry'])
    cliArgs.push('--retries=0');

  let testDirFinal = testDir;
  if (args.category) {
    const candidate = path.join(testDir, args.category);
    if (fs.existsSync(candidate))
      testDirFinal = candidate;
  }
  if (testDirFinal !== testDir)
    cliArgs.push(testDirFinal);

  const env: NodeJS.ProcessEnv = {
    ...process.env,
    DATAGROK_URL: webUrl,
    DATAGROK_AUTH_TOKEN: token,
    PLAYWRIGHT_JSON_OUTPUT_NAME: reportFile,
  };
  if (token2)
    env.DATAGROK_AUTH_TOKEN_2 = token2;

  color.info(`Playwright: running ${path.relative(pkgDir, testDir) || '.'} against ${webUrl}`);

  const stdoutChunks: Buffer[] = [];
  const stderrChunks: Buffer[] = [];
  const exitCode: number = await new Promise((resolve) => {
    const isWin = process.platform === 'win32';
    const child = spawn(isWin ? 'npx.cmd' : 'npx', cliArgs, {cwd: pkgDir, env: env, shell: isWin});
    child.stdout.on('data', (d) => {
      stdoutChunks.push(d);
      if (args.gui || args.verbose)
        process.stdout.write(d);
    });
    child.stderr.on('data', (d) => {
      stderrChunks.push(d);
      process.stderr.write(d);
    });
    child.on('error', (e) => {
      color.error(`Playwright: failed to spawn npx: ${e.message}`);
      resolve(1);
    });
    child.on('close', (code) => resolve(code ?? 1));
  });

  if (args.gui) {
    return {
      ...empty,
      failed: exitCode !== 0,
      failedAmount: exitCode !== 0 ? 1 : 0,
      passedAmount: exitCode === 0 ? 1 : 0,
      verboseFailed: exitCode !== 0 ? 'Playwright (gui mode) exited non-zero\n' : '',
    };
  }

  let report: PwReport | undefined;
  if (fs.existsSync(reportFile)) {
    try { report = JSON.parse(fs.readFileSync(reportFile, 'utf-8')); }
    catch (e: any) { color.warn(`Playwright: cannot parse JSON report: ${e.message || e}`); }
  }
  if (!report) {
    const stdoutText = Buffer.concat(stdoutChunks).toString('utf-8');
    try { report = JSON.parse(stdoutText); }
    catch { /* ignore */ }
  }

  if (!report) {
    color.error('Playwright: no JSON report produced.');
    const tail = Buffer.concat(stderrChunks).toString('utf-8').slice(-2000);
    return {
      ...empty, failed: true, failedAmount: 1,
      verboseFailed: `Playwright: no JSON report. stderr tail:\n${tail}\n`,
    };
  }

  const pkgJson = JSON.parse(fs.readFileSync(path.join(pkgDir, 'package.json'), 'utf-8'));
  const owner = (pkgJson.author && (pkgJson.author.email || pkgJson.author)) || '';
  const pkgName = process.env.TARGET_PACKAGE || args.package || pkgJson.name || '';

  const rows: FlatRow[] = [];
  flattenSuites(report.suites, testDir, pkgName, typeof owner === 'string' ? owner : '', args.verbose === true, rows);

  let passedAmount = 0;
  let failedAmount = 0;
  let skippedAmount = 0;
  let verbosePassed = '';
  let verboseFailed = '';
  let verboseSkipped = '';
  for (var r of rows) {
    const line = `${r.category}: ${r.name} (${r.ms} ms)\n`;
    if (r.skipped) { skippedAmount++; verboseSkipped += line; }
    else if (r.success) { passedAmount++; verbosePassed += line; }
    else { failedAmount++; verboseFailed += `${r.category}: ${r.name} (${r.ms} ms) :  ${r.result}\n`; }
  }

  const csv = rowsToCsv(rows);
  // Persist a Playwright-only CSV so the pipeline can ship it to the Datlas
  // 'playwright' bucket, separate from the merged Puppeteer+Playwright
  // test-report.csv that feeds the legacy 'package' bucket and JUnit.
  try {
    fs.writeFileSync(path.join(pkgDir, 'test-report-playwright.csv'), csv, 'utf-8');
  } catch (e: any) {
    color.warn(`Playwright: failed to write test-report-playwright.csv: ${e.message || e}`);
  }

  return {
    failed: failedAmount > 0,
    passedAmount: passedAmount,
    failedAmount: failedAmount,
    skippedAmount: skippedAmount,
    verbosePassed: verbosePassed,
    verboseFailed: verboseFailed,
    verboseSkipped: verboseSkipped,
    csv: csv,
  };
}
