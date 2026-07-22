/* eslint-disable max-len */
/**
 * CLI side of the `grok test` Node (browserless) pass: spawns node-test-worker.js
 * with cwd = the package directory, streams its progress output through the shared
 * TestProgressReporter, and returns a ResultObject-compatible report.
 *
 * Returns undefined when the package doesn't support the Node pass (no testNode export,
 * no datagrok-api Node runtime, not published) or the worker fails — the caller then
 * runs everything in the browser as before.
 */
import {spawn} from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import readline from 'readline';
import * as color from './color-utils';
import {getDevKey, getToken, ResultObject, TestProgressReporter} from './test-utils';

const workerTimeout = 3600000;

export interface NodeTestResult extends ResultObject {
  browserTestsRemaining: number;
  nodeTestsRun: number;
}

export interface NodeTestOptions {
  packageDir: string;
  packageName: string;
  host: string;
  category?: string;
  test?: string;
  stressTest?: boolean;
  benchmark?: boolean;
  verbose?: boolean;
}

export async function runNodeTests(opts: NodeTestOptions): Promise<NodeTestResult | undefined> {
  let url: string; let key: string;
  try {
    ({url, key} = getDevKey(opts.host));
  } catch (e: any) {
    color.warn(`Node test pass skipped: ${e?.message ?? e}`);
    return undefined;
  }
  const token = await getToken(url, key);

  const workerPath = path.join(__dirname, 'node-test-worker.js');
  const outPath = path.join(os.tmpdir(), `grok-node-test-${Date.now()}-${process.pid}.json`);
  const args = [workerPath, `--out=${outPath}`];
  if (opts.category) args.push(`--category=${opts.category}`);
  if (opts.test) args.push(`--test=${opts.test}`);
  if (opts.stressTest) args.push('--stress-test');
  if (opts.benchmark) args.push('--benchmark');
  if (opts.verbose) args.push('--verbose');

  color.info('Running Node (browserless) test pass...');
  const reporter = new TestProgressReporter(opts.verbose ?? false);
  const stderrTail: string[] = [];

  const exitCode = await new Promise<number>((resolve) => {
    const child = spawn(process.execPath, args, {
      cwd: opts.packageDir,
      env: {...process.env, DG_API_URL: url, DG_API_TOKEN: token, TARGET_PACKAGE: opts.packageName},
      stdio: ['ignore', 'pipe', 'pipe'],
    });
    const killTimer = setTimeout(() => {
      color.error(`Node test pass timed out after ${workerTimeout / 1000}s. Killing worker.`);
      child.kill();
    }, workerTimeout);
    const rlOut = readline.createInterface({input: child.stdout!});
    rlOut.on('line', (line: string) => {
      if (line.startsWith('Package testing: '))
        reporter.onLine(line);
      else if (line.startsWith('Node pass unavailable: '))
        color.info(line);
      else if (opts.verbose)
        console.log(line);
    });
    const rlErr = readline.createInterface({input: child.stderr!});
    rlErr.on('line', (line: string) => {
      stderrTail.push(line);
      if (stderrTail.length > 50)
        stderrTail.shift();
      if (opts.verbose)
        console.error(line);
    });
    child.on('close', (code) => {
      clearTimeout(killTimer);
      rlOut.close();
      rlErr.close();
      resolve(code ?? 2);
    });
    child.on('error', (e) => {
      clearTimeout(killTimer);
      color.warn(`Failed to spawn Node test worker: ${e.message}`);
      resolve(2);
    });
  });
  reporter.finish();

  if (exitCode === 3) {
    color.info('Node test pass is not supported here; running all tests in the browser.');
    return undefined;
  }
  if (exitCode !== 0) {
    if (!opts.verbose && stderrTail.length > 0)
      console.error(stderrTail.join('\n'));
    color.warn(`Node test pass failed (exit code ${exitCode}); running all tests in the browser.`);
    return undefined;
  }

  try {
    const report = JSON.parse(fs.readFileSync(outPath, 'utf8'));
    fs.unlinkSync(outPath);
    color.info(`Node pass: ${report.passedAmount} passed, ${report.failedAmount} failed, ${report.skippedAmount} skipped; ` +
      `${report.browserTestsRemaining} test(s) left for the browser.`);
    return {...report, modernOutput: true};
  } catch (e: any) {
    color.warn(`Could not read Node test report: ${e?.message ?? e}; running all tests in the browser.`);
    return undefined;
  }
}
