#!/usr/bin/env node
// One command for the dev loop: compile, run the unit tier, and (unless --fast) prove the
// runtime still serves a real turn end to end. Exits non-zero on any failure.
//
//   node dev/check.mjs          compile + unit + live smoke turn
//   node dev/check.mjs --fast   compile + unit only (no model call, no quota spent)

import {execFileSync, spawnSync} from 'node:child_process';
import * as path from 'node:path';
import {pathToFileURL} from 'node:url';

const DEV = import.meta.dirname;
const RUNTIME = path.join(DEV, '..', 'dockerfiles', 'claude-runtime');
const fast = process.argv.includes('--fast');

function step(label, fn) {
  const t0 = Date.now();
  process.stdout.write(`${label} … `);
  try {
    const detail = fn() ?? '';
    console.log(`ok ${detail} (${Date.now() - t0}ms)`);
  } catch (e) {
    console.log('FAILED');
    console.error(e.stdout?.toString() ?? e.message);
    process.exit(1);
  }
}

step('compile', () => {
  execFileSync('npx', ['tsc'], {cwd: RUNTIME, encoding: 'utf8', shell: true});
});

step('unit tests', () => {
  const out = execFileSync('node', ['--test', 'test/'], {cwd: RUNTIME, encoding: 'utf8'});
  const pass = out.match(/^# pass (\d+)$/m)?.[1] ?? '?';
  return `(${pass} tests)`;
});

// Cheap to check, expensive to get wrong: a malformed assert only surfaces ~20 minutes into a run.
step('benchmark suite', () => {
  const harness = path.join(DEV, 'harness');
  const out = execFileSync('node', ['lint-suite.mjs'], {cwd: harness, encoding: 'utf8'});
  return `(${out.match(/^suite: (\d+) prompts$/m)?.[1] ?? '?'} prompts)`;
});

if (fast)
  process.exit(0);

// The container has the old modules in memory until restarted, even though dist/ is bind-mounted.
step('restart runtime', () => {
  execFileSync('node', [path.join(DEV, 'runtime.mjs'), 'up'], {encoding: 'utf8'});
});

step('live turn', () => {
  // A Windows absolute path is not a valid ESM specifier — dynamic import() needs a file:// URL.
  const driveUrl = pathToFileURL(path.join(DEV, 'harness', 'drive.mjs')).href;
  const script = `
    import('${driveUrl}').then(async ({driveOnce}) => {
      const r = await driveOnce('Reply with exactly: PONG',
        {systemPromptMode: 'none', model: 'haiku', timeoutMs: 90000});
      if (!r.ok || !/PONG/.test(r.content)) {
        console.error('unexpected turn result:', JSON.stringify(r).slice(0, 400));
        process.exit(1);
      }
      console.log(r.totalMs);
      process.exit(0);
    });`;
  const r = spawnSync('node', ['-e', script], {encoding: 'utf8'});
  if (r.status !== 0)
    throw new Error(r.stderr || r.stdout);
  return `(${r.stdout.trim()}ms round trip)`;
});

console.log('\nall green');
