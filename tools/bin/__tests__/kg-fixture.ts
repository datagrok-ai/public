/**
 * The one harness the `grok kg` tests share (kg-codex-review-5.md step 4): a throwaway copy of a fixture
 * monorepo, the command run with its console captured, and readers over the generation a build wrote.
 * Not a test file itself, so the include glob skips it.
 */
import {expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {spawnSync} from 'child_process';
import {fileURLToPath} from 'url';
import {kg} from '../commands/kg';
import {currentDir} from '../utils/kg/generation';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';

export const FIXTURES = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg');
export const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
/** `build` is the mini monorepo with packages, Dart, docs, tests and a backlog; `good` the smaller one for the home layer. */
export type Fixture = 'build' | 'good';
const FIXTURE_DATE = '2026-01-01T00:00:00Z';

/** A throwaway copy of the fixture, so a test can break or extend one file without touching the original. */
export function copyFixture(fixture: Fixture, prepare?: (repo: string) => void): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), `grok-kg-${fixture}-`));
  fs.cpSync(path.join(FIXTURES, fixture), repo, {recursive: true});
  prepare?.(repo);
  return repo;
}

/** Runs git in [cwd] with a fixed identity and fixed dates, so a fixture commit hashes the same on every run. */
export function git(cwd: string, ...args: string[]): string {
  const env = {...process.env, GIT_AUTHOR_DATE: FIXTURE_DATE, GIT_COMMITTER_DATE: FIXTURE_DATE};
  return spawnSync('git', ['-c', 'user.email=kg@test', '-c', 'user.name=kg', '-C', cwd, ...args], {encoding: 'utf8', env}).stdout.trim();
}

export function kgRoot(repo: string): string {
  return path.join(repo, KG_DIR);
}

export function fixtureTypes(fixture: Fixture): TypeSystem {
  return loadTypeSystem(path.join(FIXTURES, fixture, KG_DIR));
}

export function write(root: string, file: string, text: string): void {
  const full = path.join(root, ...file.split('/'));
  fs.mkdirSync(path.dirname(full), {recursive: true});
  fs.writeFileSync(full, text);
}

export interface Run {
  ok: boolean;
  out: string[];
  err: string[];
  exitCode: number | undefined;
}

/** Runs the command with console captured; resets the exit code it may have set. */
export async function runKg(argv: Record<string, unknown>): Promise<Run> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    const ok = await kg(argv);
    return {ok, out: log.mock.calls.map((c) => String(c[0])), err: error.mock.calls.map((c) => String(c[0])), exitCode: process.exitCode as number | undefined};
  }
  finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

/** The rows of `nodes/<type>`, `edges/<name>` or `reports/<file>.jsonl` under a generation; none when the file is absent. */
export function readRows(genDir: string, file: string): any[] {
  const name = file.endsWith('.jsonl') ? file : `${file}.jsonl`;
  const p = path.join(genDir, file.startsWith('reports/') ? name : `data/${name}`);
  return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
}

export function readReport(genDir: string, name: string): any {
  const p = path.join(genDir, 'reports', `${name}.json`);
  return fs.existsSync(p) ? JSON.parse(fs.readFileSync(p, 'utf8')) : undefined;
}

export interface Built {
  repo: string;
  /** The output root the build wrote under (`<repo>/.kg` unless `out` or `public` said otherwise). */
  root: string;
  /** The generation directory. */
  out: string;
  manifest: any;
  rows: (file: string) => any[];
  report: (name: string) => any;
  problems: Record<string, string[]>;
}

/** Builds [repo] without an index, narrowed to [only] when given, with the fixture's `backlog/` when it has one;
 * a build that prints an error or sets an exit code fails the test here. */
export async function buildFixture(repo: string, only?: string, extra: Record<string, unknown> = {}): Promise<Built> {
  const backlog = path.join(repo, 'backlog');
  const argv: Record<string, unknown> = {_: ['kg', 'build'], kg: kgRoot(repo), db: false, output: 'json', ...extra};
  if (only !== undefined) argv.only = only;
  if (fs.existsSync(backlog) && argv.backlog === undefined) argv.backlog = backlog;
  const result = await runKg(argv);
  expect(result.err).toEqual([]);
  expect(result.exitCode).toBeUndefined();
  const root = typeof extra.out === 'string' ? extra.out : path.join(repo, ...(extra.public ? ['public', '.kg'] : ['.kg']));
  const out = currentDir(root)!;
  return {repo, root, out, manifest: JSON.parse(result.out[0]), rows: (file) => readRows(out, file), report: (name) => readReport(out, name),
    problems: readReport(out, 'problems') ?? {}};
}
