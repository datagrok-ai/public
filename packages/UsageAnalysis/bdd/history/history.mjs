#!/usr/bin/env node
/* The run history of the bdd suites. Every `grok-bdd run` leaves its Playwright JSON report in the
   project's bdd/test-results/report.json; `record` keeps the reports of the last full runs as one
   dated record in runs/ (where they ran included: machine and stand), and `html` turns every record
   there into one standalone page. Nothing records on its own — a record is made when asked for.

     node history.mjs record [--note <text>] [--partial] [<report.json | dir>...]
     node history.mjs html [--out <file>]
*/
import {execSync} from 'node:child_process';
import {existsSync, mkdirSync, readdirSync, readFileSync, statSync, writeFileSync} from 'node:fs';
import os from 'node:os';
import {dirname, join, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const RUNS = join(HERE, 'runs');
const PACKAGES = resolve(HERE, '../../..');
// a project whose last report is older than the newest one by more than this was not part of the run
const SAME_RUN_MS = 12 * 3600 * 1000;
const GHERKIN = /^(Given|When|Then|And|But|\*)\s/;

function option(args, name) {
  const i = args.indexOf(name);
  return i < 0 ? undefined : args.splice(i, 2)[1];
}

function flag(args, name) {
  const i = args.indexOf(name);
  return i >= 0 && args.splice(i, 1).length > 0;
}

function reportFiles(sources) {
  if (sources.length === 0)
    return readdirSync(PACKAGES).map((p) => join(PACKAGES, p, 'bdd', 'test-results', 'report.json')).filter(existsSync);
  return sources.flatMap((s) => statSync(s).isDirectory()
    ? readdirSync(s).filter((f) => f.endsWith('.json')).map((f) => join(s, f)) : [s]);
}

function specFiles(dir) {
  return existsSync(dir) ? readdirSync(dir, {recursive: true}).filter((f) => String(f).endsWith('.test.ts')).length : 0;
}

const clean = (s) => s?.replace(/\u001b\[[0-9;]*m/g, '').trim();

/** The deepest step that failed: the gesture or check the test stopped at. */
function failedStep(steps) {
  for (const s of steps ?? [])
    if (s.error)
      return failedStep(s.steps) ?? s.title;
}

function testRecord(file, spec, test, t0) {
  const last = test.results.at(-1);
  const tags = spec.tags.filter((t) => !t.startsWith('realizes:'));
  const error = clean(last?.errors?.[0]?.message ?? last?.error?.message);
  const top = last?.steps ?? [];
  const groups = tags.includes('journey') ? top.filter((s) => !GHERKIN.test(s.title)) : [];
  const status = test.status === 'skipped' ? 'skipped' : test.status === 'flaky' ? 'flaky' : test.status === 'unexpected' ? 'failed' :
    groups.length === 0 && (tags.includes('known-failure') || test.expectedStatus === 'failed') ? 'known' : 'passed';
  const r = {file, title: spec.title, status, duration: Math.round(last?.duration ?? 0)};
  if (tags.length > 0)
    r.tags = tags;
  if (test.results.length > 1)
    r.retries = test.results.length - 1;
  if (last) {
    r.start = Math.max(0, Math.round(Date.parse(last.startTime) - t0));
    r.worker = last.parallelIndex;
  }
  if (groups.length > 0) {
    r.setup = Math.round(top.filter((s) => GHERKIN.test(s.title)).reduce((sum, s) => sum + s.duration, 0));
    // a failed journey lists its failed scenarios by name; an erring scenario it does not list is a
    // known failure; a passed one, most of them, leaves its status out
    r.scenarios = groups.map((s) => {
      const outcome = error?.includes(`\n\n${s.title}\n`) ? 'failed' : s.error ? 'known' : undefined;
      return {title: s.title, duration: Math.round(s.duration), ...(outcome ? {status: outcome} : {})};
    });
  }
  if (error && status !== 'passed' && status !== 'known') {
    r.error = error.length > 600 ? `${error.slice(0, 600)}…` : error;
    const failed = groups.find((s) => r.scenarios.find((x) => x.title === s.title).status === 'failed');
    r.step = failedStep(failed ? [failed] : top);
  }
  return r;
}

function grepInvert(argv = []) {
  const i = argv.findIndex((a) => a === '--grep-invert' || a.startsWith('--grep-invert='));
  return i < 0 ? undefined : argv[i] === '--grep-invert' ? argv[i + 1] : argv[i].slice('--grep-invert='.length);
}

function projectRecord(report) {
  const t0 = Date.parse(report.stats.startTime);
  const tests = [];
  const walk = (suite, file) => {
    for (const spec of suite.specs ?? [])
      for (const test of spec.tests)
        tests.push(testRecord(file, spec, test, t0));
    for (const s of suite.suites ?? [])
      walk(s, file);
  };
  for (const s of report.suites)
    walk(s, s.file.replace(/\\/g, '/').replace(/\.test\.ts$/, ''));
  const root = report.config.rootDir;
  return {
    name: /packages[\\/]([^\\/]+)[\\/]bdd/.exec(root)?.[1] ?? root,
    start: report.stats.startTime,
    duration: Math.round(report.stats.duration),
    workers: report.config.metadata?.actualWorkers ?? report.config.workers,
    // a run narrowed by a tag (`--grep-invert @full-stand` on a stand without that capability) is
    // still the full run of what the stand can run; the report keeps the option only in its argv
    ...(grepInvert(report.config.argv) ? {excluded: grepInvert(report.config.argv)} : {}),
    specs: new Set(report.suites.map((s) => s.file)).size,
    of: specFiles(root),
    tests,
  };
}

function localMachine() {
  const cpus = os.cpus();
  return {host: os.hostname(), os: `${os.platform()} ${os.release()}`, cpu: cpus[0]?.model.trim(), cores: cpus.length,
    memoryGb: Math.round(os.totalmem() / 2 ** 30), node: process.version};
}

function git() {
  const run = (cmd) => execSync(cmd, {cwd: PACKAGES, encoding: 'utf8'}).trim();
  try {
    return {branch: run('git rev-parse --abbrev-ref HEAD'), commit: run('git rev-parse --short HEAD'),
      dirty: run('git status --porcelain --untracked-files=no') !== ''};
  }
  catch {
    return undefined;
  }
}

function record(args) {
  const note = option(args, '--note');
  const partial = flag(args, '--partial');
  const reports = reportFiles(args).map((f) => ({f, report: JSON.parse(readFileSync(f, 'utf8'))}));
  if (reports.length === 0)
    throw new Error('no report to record: run the suites first (grok-bdd run leaves bdd/test-results/report.json)');
  const newest = Math.max(...reports.map((r) => Date.parse(r.report.stats.startTime)));
  const projects = [];
  for (const {report} of reports) {
    const p = projectRecord(report);
    if (newest - Date.parse(p.start) > SAME_RUN_MS)
      console.log(`skipped ${p.name}: its last run (${p.start}) is not part of this one`);
    else if (p.specs < p.of && !p.excluded && !partial)
      console.log(`skipped ${p.name}: its last run covered ${p.specs} of ${p.of} specs (--partial records it anyway)`);
    else
      projects.push(p);
  }
  if (projects.length === 0)
    throw new Error('nothing to record');
  projects.sort((a, b) => a.start.localeCompare(b.start));
  const meta = reports[0].report.config.metadata ?? {};
  const rec = {version: 1, date: projects[0].start, recorded: new Date().toISOString(), note,
    machine: meta.machine ?? localMachine(), stand: meta.stand ?? process.env.DATAGROK_URL ?? 'http://localhost:8888',
    git: git(), projects};
  mkdirSync(RUNS, {recursive: true});
  const name = `${rec.date.slice(0, 16).replace(/:/g, '-')}-${rec.machine.host}.json`;
  // one test per line: a record stays readable in a diff without doubling in size
  const text = JSON.stringify(rec, (k, v) => k === 'tests' ? v.map((t) => `\u0000${JSON.stringify(t)}\u0000`) : v, 2)
    .replace(/"\\u0000(.*)\\u0000"/g, (_, s) => JSON.parse(`"${s}"`));
  writeFileSync(join(RUNS, name), `${text}\n`);
  const count = (s) => projects.reduce((n, p) => n + p.tests.filter((t) => t.status === s).length, 0);
  console.log(`recorded runs/${name}: ${projects.map((p) => p.name).join(', ')} — ` +
    `${count('passed')} passed, ${count('failed')} failed, ${count('flaky')} flaky, ${count('known')} known-failure tests`);
}

function html(args) {
  const out = resolve(option(args, '--out') ?? join(HERE, 'history.html'));
  const records = existsSync(RUNS) ? readdirSync(RUNS).filter((f) => f.endsWith('.json'))
    .map((f) => JSON.parse(readFileSync(join(RUNS, f), 'utf8'))).sort((a, b) => a.date.localeCompare(b.date)) : [];
  const page = readFileSync(join(HERE, 'template.html'), 'utf8')
    .replace('/*HISTORY*/[]', () => JSON.stringify(records).replace(/</g, '\\u003c'));
  writeFileSync(out, page);
  console.log(`${out}: ${records.length} record(s)`);
}

const [command, ...args] = process.argv.slice(2);
if (command === 'record')
  record(args);
else if (command === 'html')
  html(args);
else {
  console.error('usage: node history.mjs record [--note <text>] [--partial] [<report.json | dir>...]\n' +
    '       node history.mjs html [--out <file>]');
  process.exitCode = 2;
}
