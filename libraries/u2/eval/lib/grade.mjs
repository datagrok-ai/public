// Grades each solution in a process of its own (a fresh page, stub and u2) and typechecks them all at once.
import {spawn} from 'node:child_process';
import {existsSync, readdirSync, readFileSync, statSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {typecheck} from './typecheck.mjs';

export const EVAL = fileURLToPath(new URL('..', import.meta.url));
const U2 = join(EVAL, '..');
const LIB = join(EVAL, 'lib');
const TIMEOUT_MS = 180_000;

/** [{id: the folder name, dir, meta: task.json}]: the single tasks (`tNN`), then the feature tasks (`fNN`). */
export function tasks() {
  const root = join(EVAL, 'tasks');
  const rank = (id) => (id.startsWith('t') ? '0' : '1') + id;
  return readdirSync(root).filter((d) => existsSync(join(root, d, 'task.json')))
    .sort((a, b) => rank(a).localeCompare(rank(b)))
    .map((id) => ({id, dir: join(root, id), meta: JSON.parse(readFileSync(join(root, id, 'task.json'), 'utf8'))}));
}

/** Graders run u2's test build (`npm test` compiles it into .test-build): refuses one older than its sources or CSS. */
export function assertFreshBuild() {
  const stale = [];
  const src = join(U2, 'src');
  const built = join(U2, '.test-build', 'src');
  const walk = (dir) => {
    for (const e of readdirSync(dir, {withFileTypes: true})) {
      const p = join(dir, e.name);
      if (e.isDirectory())
        walk(p);
      else if (e.name.endsWith('.ts') && !e.name.endsWith('.d.ts') && !e.name.endsWith('.g.ts')) {
        const js = join(built, p.slice(src.length)).replace(/\.ts$/, '.js');
        if (!existsSync(js) || statSync(js).mtimeMs < statSync(p).mtimeMs)
          stale.push(p);
      }
    }
  };
  walk(src);
  const styles = join(built, 'styles.g.js');
  const css = readdirSync(join(U2, 'css')).filter((f) => f.endsWith('.css')).map((f) => join(U2, 'css', f));
  if (!existsSync(styles) || css.some((f) => statSync(f).mtimeMs > statSync(styles).mtimeMs))
    stale.push(join(U2, 'css'));
  if (stale.length > 0)
    throw new Error(`${stale.length} u2 sources are newer than their test build (${stale[0]}, …): ` +
      'run `npm test` in libraries/u2 first');
}

function child(grader, solution) {
  return new Promise((resolve) => {
    const args = ['--import', join(LIB, 'register.mjs'), join(LIB, 'child.mjs'), grader, solution];
    const proc = spawn(process.execPath, args, {cwd: EVAL, stdio: ['ignore', 'pipe', 'pipe'],
      env: {...process.env, U2_STRICT: process.env.U2_STRICT ?? '1'}});
    let out = '';
    let err = '';
    proc.stdout.on('data', (d) => out += d);
    proc.stderr.on('data', (d) => err += d);
    const timer = setTimeout(() => proc.kill('SIGKILL'), TIMEOUT_MS);
    proc.on('close', (code) => {
      clearTimeout(timer);
      const marked = (tag) => out.split('\n').reverse().find((l) => l.startsWith(tag))?.slice(tag.length);
      const verdicts = marked('@@EVAL ');
      if (verdicts !== undefined)
        return resolve(JSON.parse(verdicts));
      const said = (err || out).split('\n').filter((l) => l.trim() !== '' && !l.startsWith('@@'));
      const how = code ?? `killed after ${TIMEOUT_MS / 1000} s`;
      const detail = `the grading process ended (${how}) without verdicts` +
        (said.length > 0 ? `: ${said.slice(-3).join(' | ')}` : '');
      resolve(['loads', ...JSON.parse(marked('@@CHECKS ') ?? '[]')].map((id) => ({id, pass: false, detail})));
    });
  });
}

/** Grades [{task, file}]: each answers its checks, passed, total and score (the fraction passed). */
export async function gradeAll(jobs, {concurrency = 4} = {}) {
  const types = typecheck(jobs.map((j) => j.file));
  const results = new Array(jobs.length);
  let next = 0;
  const worker = async () => {
    while (next < jobs.length) {
      const i = next++;
      const {task, file} = jobs[i];
      const checks = (await child(join(task.dir, 'grader.mjs'), file)).map((c) => c.detail === undefined ? c :
        {...c, detail: c.detail.replaceAll(pathToFileURL(U2).href + '/', 'u2/').replaceAll(U2 + '/', 'u2/')});
      const problems = types.get(file) ?? [];
      checks.push(problems.length === 0 ? {id: 'typecheck', pass: true} :
        {id: 'typecheck', pass: false, detail: problems.slice(0, 3).join(' / ')});
      const passed = checks.filter((c) => c.pass).length;
      results[i] = {task: task.id, file, checks, passed, total: checks.length, score: passed / checks.length};
    }
  };
  await Promise.all(Array.from({length: concurrency}, worker));
  return results;
}

export const failedIds = (result) => result.checks.filter((c) => !c.pass).map((c) => c.id).sort();

export function u2Build() {
  const pkg = JSON.parse(readFileSync(join(U2, 'package.json'), 'utf8'));
  const lock = JSON.parse(readFileSync(join(U2, 'css', 'styles.lock.json'), 'utf8'));
  return {version: pkg.version, stylesRevision: lock.revision};
}

export const today = () => new Date().toISOString().slice(0, 10);
