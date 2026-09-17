import {describe, it, expect, afterEach, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {findMonorepoRoot, parsePlan, filterRuns, splitShellWords, runPlan, planFailed, childCommand, kgArgs, extractJson, childCsvPath, printPlan, RunRow, Spawner} from '../utils/recent-tests';

const dirs: string[] = [];
function tempDir(name: string): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), `grok-recent-${name}-`));
  dirs.push(dir);
  return dir;
}

afterEach(() => {
  vi.restoreAllMocks();
  for (const dir of dirs.splice(0))
    fs.rmSync(dir, {recursive: true, force: true});
});

const todayShape = {
  op: 'tests-for',
  target: {id: 'file:public/packages/Chem/src/analysis/activity-cliffs.ts'},
  sections: [
    {title: 'features', rows: [{feature: 'domains/chem', name: 'Cheminformatics'}]},
    {title: 'tests', rows: [{framework: 'dg', level: 'unit', test: 'test:dg:public/packages/Chem/src/tests/api-based-tests.ts#chem exported/mcs', feature: 'domains/chem'}]},
    {title: 'scenarios', rows: []},
    {title: 'automations', rows: []},
  ],
};

const runRows: RunRow[] = [
  {framework: 'dg', cwd: 'public/packages/Chem', command: 'grok test --category "chem exported" --test mcs', tests: 1, names: ['chem exported/mcs']},
  {framework: 'xamgle', cwd: 'public/packages/DevTools', command: 'grok test --category "Core: xamgle" --test "Spaces | Tree"', tests: 1},
  {framework: 'dart', cwd: 'core/shared/ddt', command: 'pub run test test/data_frame_test.dart -n "sorting"', tests: 1},
  {framework: 'node', cwd: 'public/tools', command: 'npx vitest run bin/__tests__/kg.test.ts', tests: 3},
];

const contractShape = {
  op: 'tests-for',
  targets: ['file:public/packages/Chem/src/analysis/activity-cliffs.ts'],
  sections: [
    {title: 'changes', rows: [{path: 'public/packages/Chem/src/analysis/activity-cliffs.ts', known: true, repo: 'public'}]},
    {title: 'immediate', rows: [{framework: 'dg', level: 'unit', test: 'test:dg:public/packages/Chem/src/tests/menu-tests-cliffs.ts#chem/cliffs', name: 'cliffs', path: 'public/packages/Chem/src/tests/menu-tests-cliffs.ts', category: 'chem', tier: 'immediate'}]},
    {title: 'reachable', rows: []},
    {title: 'feature', rows: [{framework: 'dart', level: 'unit', test: 'test:dart:core/shared/ddt/test/data_frame_test.dart#sorting', name: 'sorting', path: 'core/shared/ddt/test/data_frame_test.dart', tier: 'feature'}]},
    {title: 'run', rows: runRows},
  ],
};

describe('findMonorepoRoot', () => {
  it('walks up to the directory holding core/ and public/', () => {
    const root = tempDir('root');
    fs.mkdirSync(path.join(root, 'core'));
    fs.mkdirSync(path.join(root, 'public', 'packages', 'Chem'), {recursive: true});
    expect(findMonorepoRoot(path.join(root, 'public', 'packages', 'Chem'))).toBe(root);
    expect(findMonorepoRoot(root)).toBe(root);
  });

  it('returns undefined when no ancestor qualifies', () => {
    const dir = tempDir('plain');
    fs.mkdirSync(path.join(dir, 'public'));
    expect(findMonorepoRoot(dir)).toBeUndefined();
  });
});

describe('kgArgs and extractJson', () => {
  it('spells --changed with and without a base, and the tier', () => {
    // immediate is the default since the unit rule (change-tests work order E): the other units' tests are --tier linked
    expect(kgArgs(true)).toEqual(['kg', 'tests-for', '--changed', '--tier', 'immediate', '--output', 'json']);
    expect(kgArgs('HEAD~3', 'immediate,feature')).toEqual(['kg', 'tests-for', '--changed=HEAD~3', '--tier', 'immediate,feature', '--output', 'json']);
  });

  it('finds the JSON object after warning lines and tolerates none', () => {
    expect(extractJson('Dart coverage partial\n{"op":"tests-for","sections":[]}\n')).toEqual({op: 'tests-for', sections: []});
    expect(extractJson('no graph found: run grok kg build')).toBeUndefined();
    expect(extractJson('{not json')).toBeUndefined();
  });
});

describe('parsePlan', () => {
  it('yields an empty plan with a message for today\'s shape without a run section', () => {
    const plan = parsePlan(todayShape);
    expect(plan.runs).toEqual([]);
    expect(plan.message).toMatch(/without a run section/);
    expect(plan.tests).toEqual([{tier: 'feature', framework: 'dg', name: 'chem exported/mcs', path: 'public/packages/Chem/src/tests/api-based-tests.ts'}]);
  });

  it('reads tiers, changes and run rows from the contract shape', () => {
    const plan = parsePlan(contractShape);
    expect(plan.changes).toBe(1);
    expect(plan.tiers).toEqual(['immediate', 'reachable']);
    expect(plan.totals).toEqual({immediate: 1, feature: 1});
    expect(parsePlan({sections: [{title: 'changes', rows: [{path: 'a'}], total: 602}, {title: 'immediate', rows: [], total: 80}, {title: 'run', rows: []}]}).totals).toEqual({immediate: 80});
    expect(plan.tests.map((t) => t.tier)).toEqual(['immediate', 'feature']);
    expect(plan.runs).toHaveLength(4);
    expect(plan.runs[0]).toEqual(runRows[0]);
    expect(plan.message).toBeUndefined();
  });

  it('says so when the run section is present but empty', () => {
    expect(parsePlan({sections: [{title: 'changes', rows: []}, {title: 'run', rows: []}]}).message).toMatch(/no changed files/);
    expect(parsePlan({sections: [{title: 'changes', rows: [{path: 'a.ts'}]}, {title: 'run', rows: []}]}).message).toMatch(/no runner rows/);
    expect(parsePlan(undefined).message).toMatch(/no sections/);
  });
});

describe('printPlan', () => {
  it('heads with the tiers, shows the first three names and cuts a long command', () => {
    const lines: string[] = [];
    vi.spyOn(console, 'log').mockImplementation((line?: unknown) => { lines.push(String(line)); });
    const long = {framework: 'dart', cwd: 'core/client/d4', command: `pub run test ${'test/legend/legend_metrics_test.dart '.repeat(5)}`, tests: 5, names: ['a', 'b', 'c', 'd', 'e']};
    printPlan({...parsePlan(contractShape, ['immediate']), runs: [runRows[0], long]});
    expect(lines[0]).toBe('Changed files: 1; tiers: immediate; linked tests: 2; runner rows: 2');
    const table = lines.join('\n');
    expect(table).toContain('a, b, c (+2 more)');
    expect(table).toContain('chem exported/mcs');
    expect(table).toContain(`${long.command.slice(0, 119)}…`);
    expect(table).not.toContain(long.command);
  });
});

describe('filterRuns', () => {
  it('keeps the frameworks listed', () => {
    expect(filterRuns(runRows, {framework: 'dart,node'}).map((r) => r.framework)).toEqual(['dart', 'node']);
    expect(filterRuns(runRows, {})).toHaveLength(4);
  });

  it('keeps the rows of a package by cwd or command', () => {
    expect(filterRuns(runRows, {package: 'Chem'}).map((r) => r.cwd)).toEqual(['public/packages/Chem']);
    expect(filterRuns(runRows, {package: 'chem'})).toHaveLength(1);
    const byCommand: RunRow = {framework: 'dg', cwd: 'public/packages', command: 'grok test --package=Chem --category x', tests: 1};
    expect(filterRuns([byCommand], {package: 'Chem'})).toHaveLength(1);
    expect(filterRuns([byCommand], {package: 'Chemistry'})).toHaveLength(0);
  });
});

describe('splitShellWords', () => {
  it('honours double quotes and drops them', () => {
    expect(splitShellWords('grok test --category "Core: xamgle" --test "Spaces | Tree"'))
      .toEqual(['grok', 'test', '--category', 'Core: xamgle', '--test', 'Spaces | Tree']);
    expect(splitShellWords('  dart test  a.dart -n \'x y\' ')).toEqual(['dart', 'test', 'a.dart', '-n', 'x y']);
    expect(splitShellWords('--test ""')).toEqual(['--test', '']);
  });
});

describe('childCommand', () => {
  const ctx = {root: '/repo', grokScript: '/repo/public/tools/bin/grok.js', passThrough: ['--host=dev', '--skip-build'], csv: path.join('out', 'report.csv')};

  it('runs grok rows through this grok.js with pass-through flags and a per-child csv', () => {
    const child = childCommand(runRows[1], 2, ctx);
    expect(child.file).toBe(process.execPath);
    expect(child.shell).toBe(false);
    expect(child.args).toEqual([ctx.grokScript, 'test', '--category', 'Core: xamgle', '--test', 'Spaces | Tree', '--host=dev', '--skip-build', `--csv=${childCsvPath(ctx.csv, 2)}`]);
    expect(child.csv).toBe(path.join('out', 'report-2.csv'));
  });

  it('runs other rows through the shell verbatim', () => {
    expect(childCommand(runRows[2], 3, ctx)).toEqual({file: runRows[2].command, args: [], shell: true});
  });
});

describe('runPlan', () => {
  it('runs rows in order, continues on failure and reports missing tools and directories', async () => {
    const root = tempDir('run');
    for (const row of runRows)
      if (row.framework !== 'node')
        fs.mkdirSync(path.join(root, row.cwd), {recursive: true});
    const calls: {file: string, cwd: string}[] = [];
    const spawner: Spawner = async (file, args, options) => {
      calls.push({file, cwd: options.cwd});
      if (file.startsWith('pub'))
        return {code: null, error: 'not found: pub'};
      return {code: args.includes('--category') && args.includes('chem exported') ? 1 : 0};
    };
    vi.spyOn(console, 'log').mockImplementation(() => {});
    const results = await runPlan(runRows, spawner, {root, grokScript: 'grok.js', passThrough: []});
    expect(calls.map((c) => c.cwd)).toEqual(runRows.slice(0, 3).map((r) => path.join(root, r.cwd)));
    expect(results.map((r) => r.exit)).toEqual([1, 0, 'not found: pub', `directory not found: ${path.join(root, 'public/tools')}`]);
    expect(results.every((r) => r.seconds >= 0)).toBe(true);
    expect(planFailed(results)).toBe(true);
    expect(planFailed(results.slice(1, 2))).toBe(false);
  });
});
