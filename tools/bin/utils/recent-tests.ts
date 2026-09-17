/// The plan and runner logic behind `grok test --recent`: the graph names the tests linked to the change set
/// (`grok kg tests-for --changed --output json`), this module turns that answer into a plan, filters it, runs each
/// runner row in its own directory and sums the exit codes. Pure functions apart from `runPlan`, which takes the
/// spawner as an argument so tests never launch a process.
import fs from 'fs';
import path from 'path';
import {spawn} from 'child_process';
import * as color from './color-utils';
import {printOutput} from './server-output';
import {TIERS, LINKED, Tier} from './kg/ops';

export {parseTiers, TIER_CHOICES} from './kg/ops';
/** The tests of the changed files' own units are the point of the command; `--tier linked` adds the other units and the
 * import walk, `--tier all` or `--tier feature` the feature's blast radius. */
export const DEFAULT_TIER = 'immediate';
const PLAN_ROWS_PER_TIER = 40;
const COMMAND_WIDTH = 80;
const PLAN_COMMAND_WIDTH = 120;
const PLAN_NAMES = 3;

export interface PlanTest {
  tier: string;
  framework: string;
  name: string;
  path: string;
}

export interface RunRow {
  framework: string;
  cwd: string;
  command: string;
  tests: number;
  names?: string[];
}

export interface RecentPlan {
  changes: number;
  /** The tiers the plan was asked for. */
  tiers: Tier[];
  /** Every linked test the graph counted, per tier; `tests` holds the page it returned. */
  totals: Record<string, number>;
  tests: PlanTest[];
  runs: RunRow[];
  message?: string;
}

export interface RunResult {
  runner: string;
  cwd: string;
  command: string;
  exit: number | string;
  seconds: number;
  csv?: string;
}

export interface RunContext {
  root: string;
  grokScript: string;
  passThrough: string[];
  csv?: string;
}

export type SpawnResult = {code: number | null, error?: string};
export type Spawner = (file: string, args: string[], options: {cwd: string, shell: boolean}) => Promise<SpawnResult>;

/** The monorepo root is the nearest ancestor holding both `core/` and `public/`. */
export function findMonorepoRoot(from: string): string | undefined {
  let dir = path.resolve(from);
  for (;;) {
    if (fs.existsSync(path.join(dir, 'core')) && fs.existsSync(path.join(dir, 'public')))
      return dir;
    const parent = path.dirname(dir);
    if (parent === dir)
      return undefined;
    dir = parent;
  }
}

/** The worktree's own grok.js carries the build the graph was made with; a globally installed grok falls back to itself. */
export function resolveGrokScript(root: string, running: string): string {
  const own = path.join(root, 'public', 'tools', 'bin', 'grok.js');
  return fs.existsSync(own) ? own : running;
}

export function kgArgs(base: boolean | string, tier: string = DEFAULT_TIER): string[] {
  const changed = typeof base === 'string' && base !== 'true' ? `--changed=${base}` : '--changed';
  return ['kg', 'tests-for', changed, '--tier', tier, '--output', 'json'];
}

/** The first JSON object in the child's stdout, or undefined when it printed none. */
export function extractJson(stdout: string): any | undefined {
  const start = stdout.search(/^\s*\{/m);
  if (start < 0)
    return undefined;
  try {
    return JSON.parse(stdout.slice(start));
  }
  catch {
    return undefined;
  }
}

/** Reads the `tests-for` answer; every section is optional, and a missing `run` section is an empty plan with a message. */
export function parsePlan(answer: any, tiers: Tier[] = LINKED): RecentPlan {
  const sections: any[] = Array.isArray(answer?.sections) ? answer.sections : [];
  const rows = (title: string): any[] => {
    const section = sections.find((s) => s?.title === title);
    return Array.isArray(section?.rows) ? section.rows : [];
  };
  const total = (title: string): number => Number(sections.find((s) => s?.title === title)?.total ?? rows(title).length);
  const tests: PlanTest[] = [];
  for (const tier of [...TIERS, 'tests'])
    for (const row of rows(tier)) {
      const id = String(row.test ?? '');
      const hash = id.indexOf('#');
      tests.push({
        tier: String(row.tier ?? (tier === 'tests' ? 'feature' : tier)),
        framework: String(row.framework ?? id.split(':')[1] ?? ''),
        name: String(row.name ?? (hash >= 0 ? id.slice(hash + 1) : id)),
        path: String(row.path ?? (hash >= 0 ? id.slice(id.indexOf(':', 5) + 1, hash) : '')),
      });
    }
  const runs: RunRow[] = rows('run').map((row) => ({
    framework: String(row.framework ?? ''),
    cwd: String(row.cwd ?? '.'),
    command: String(row.command ?? ''),
    tests: Number(row.tests ?? (Array.isArray(row.names) ? row.names.length : 0)),
    names: Array.isArray(row.names) ? row.names.map(String) : undefined,
  })).filter((row) => row.command.length > 0);
  const changes = total('changes');
  const totals = Object.fromEntries([...TIERS, 'tests'].map((title) => [title, total(title)]).filter(([, n]) => n));
  const hasRun = sections.some((s) => s?.title === 'run');
  const message = !sections.length ? 'the graph returned no sections' :
    !hasRun ? 'the graph answered without a run section (this grok kg build has no runner rows); nothing to run' :
    !runs.length ? (changes || tests.length ? 'no runner rows for the change set; nothing to run' : 'no changed files; nothing to run') : undefined;
  return {changes, tiers, totals, tests, runs, message};
}

/** Keeps the runner rows a `--framework a,b` or `--package <Name>` filter names; both apply after the graph call. */
export function filterRuns(runs: RunRow[], filter: {framework?: string, package?: string}): RunRow[] {
  const frameworks = filter.framework ? String(filter.framework).split(',').map((s) => s.trim().toLowerCase()).filter((s) => s.length) : [];
  const pkg = filter.package ? String(filter.package).toLowerCase() : '';
  return runs.filter((row) => {
    if (frameworks.length && !frameworks.includes(row.framework.toLowerCase()))
      return false;
    if (!pkg)
      return true;
    const cwd = row.cwd.replace(/\\/g, '/').toLowerCase();
    const command = row.command.replace(/\\/g, '/').toLowerCase();
    return new RegExp(`/packages/${pkg}(/|$)`).test(`/${cwd}`) || command.includes(`packages/${pkg}/`) ||
      new RegExp(`--package[= ]"?${pkg}"?(\\s|$)`).test(command);
  });
}

/** Splits a paste-able command into words; double and single quotes group, and are dropped. */
export function splitShellWords(command: string): string[] {
  const words: string[] = [];
  let word = '';
  let quote = '';
  let inWord = false;
  for (const ch of command) {
    if (quote) {
      if (ch === quote)
        quote = '';
      else
        word += ch;
    }
    else if (ch === '"' || ch === '\'') {
      quote = ch;
      inWord = true;
    }
    else if (/\s/.test(ch)) {
      if (inWord)
        words.push(word);
      word = '';
      inWord = false;
    }
    else {
      word += ch;
      inWord = true;
    }
  }
  if (inWord)
    words.push(word);
  return words;
}

/** The flags a parent `grok test --recent` hands to every child `grok test`; each child gets its own csv file. */
export function passThroughArgs(args: Record<string, any>): string[] {
  const out: string[] = [];
  if (args.host) out.push(`--host=${args.host}`);
  if (args['skip-build']) out.push('--skip-build');
  if (args['skip-publish']) out.push('--skip-publish');
  if (args['no-retry']) out.push('--no-retry');
  if (args.verbose) out.push('--verbose');
  return out;
}

export function childCsvPath(csv: string, index: number): string {
  const ext = path.extname(csv);
  return path.join(path.dirname(csv), `${path.basename(csv, ext)}-${index}${ext || '.csv'}`);
}

/** What to spawn for a runner row: a `grok test` row runs this grok.js with the row's words, anything else goes through the shell. */
export function childCommand(row: RunRow, index: number, ctx: RunContext): {file: string, args: string[], shell: boolean, csv?: string} {
  const words = splitShellWords(row.command);
  if (words[0] !== 'grok')
    return {file: row.command, args: [], shell: true};
  const args = [ctx.grokScript, ...words.slice(1), ...ctx.passThrough];
  const csv = ctx.csv && words[1] === 'test' ? childCsvPath(ctx.csv, index) : undefined;
  if (csv)
    args.push(`--csv=${csv}`);
  return {file: process.execPath, args, shell: false, csv};
}

export const spawnInherited: Spawner = (file, args, options) => new Promise((resolve) => {
  const child = spawn(file, args, {cwd: options.cwd, shell: options.shell, stdio: 'inherit'});
  child.on('error', (e: any) => resolve({code: null, error: e.code === 'ENOENT' ? `not found: ${file}` : e.message}));
  child.on('close', (code) => resolve({code}));
});

/** Runs every row in order, continuing on failure; a runner whose tool or directory is missing is a failed row, not a throw. */
export async function runPlan(runs: RunRow[], spawner: Spawner, ctx: RunContext): Promise<RunResult[]> {
  const results: RunResult[] = [];
  for (const [i, row] of runs.entries()) {
    const cwd = path.join(ctx.root, row.cwd);
    color.info(`\n[${i + 1}/${runs.length}] ${row.framework}: ${row.command}  (in ${row.cwd})`);
    const started = Date.now();
    const child = childCommand(row, i + 1, ctx);
    const outcome: SpawnResult = fs.existsSync(cwd) ? await spawner(child.file, child.args, {cwd, shell: child.shell}) :
      {code: null, error: `directory not found: ${cwd}`};
    const seconds = Math.round((Date.now() - started) / 100) / 10;
    results.push({runner: row.framework, cwd: row.cwd, command: row.command, exit: outcome.error ?? outcome.code ?? 'killed', seconds, csv: child.csv});
  }
  return results;
}

export function planFailed(results: RunResult[]): boolean {
  return results.some((r) => r.exit !== 0);
}

/** A cell holds one line: the first few names, and a command cut at the width a terminal shows whole. */
function planRow(r: RunRow): Record<string, unknown> {
  const names = r.names ?? [];
  return {
    framework: r.framework, cwd: r.cwd,
    command: r.command.length > PLAN_COMMAND_WIDTH ? `${r.command.slice(0, PLAN_COMMAND_WIDTH - 1)}…` : r.command,
    tests: r.tests, names: `${names.slice(0, PLAN_NAMES).join(', ')}${names.length > PLAN_NAMES ? ` (+${names.length - PLAN_NAMES} more)` : ''}`,
  };
}

export function printPlan(plan: RecentPlan): void {
  const linked = Object.values(plan.totals).reduce((sum, n) => sum + n, 0);
  console.log(`Changed files: ${plan.changes}; tiers: ${plan.tiers.join(', ')}; linked tests: ${linked}; runner rows: ${plan.runs.length}`);
  for (const tier of TIERS) {
    const rows = plan.tests.filter((t) => t.tier === tier);
    if (!rows.length)
      continue;
    const total = plan.totals[tier] ?? rows.length;
    console.log(`\n${tier} (${rows.length < total ? `${rows.length} of ${total}` : total})`);
    printOutput(rows.slice(0, PLAN_ROWS_PER_TIER), 'table');
    if (rows.length > PLAN_ROWS_PER_TIER)
      console.log(`… ${rows.length - PLAN_ROWS_PER_TIER} more`);
  }
  if (plan.runs.length) {
    console.log('\nrun');
    printOutput(plan.runs.map(planRow), 'table', PLAN_COMMAND_WIDTH);
  }
  if (plan.message)
    color.warn(`\n${plan.message}`);
}

export function printSummary(results: RunResult[]): void {
  console.log('');
  printOutput(results.map((r) => ({
    runner: r.runner, cwd: r.cwd,
    command: r.command.length > COMMAND_WIDTH ? `${r.command.slice(0, COMMAND_WIDTH - 1)}…` : r.command,
    exit: r.exit, seconds: r.seconds,
  })), 'table');
  const csvs = results.filter((r) => r.csv).map((r) => r.csv);
  if (csvs.length)
    console.log(`\nCSV reports: ${csvs.join(', ')}`);
  const failed = results.filter((r) => r.exit !== 0).length;
  if (failed)
    color.error(`\n${failed} of ${results.length} runner(s) failed`);
  else
    color.success(`\nAll ${results.length} runner(s) passed`);
}
