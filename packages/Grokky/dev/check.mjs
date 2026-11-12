#!/usr/bin/env node
// One command for the dev loop: compile and run the unit tier. Exits non-zero on any failure.
// A live turn needs model credentials, which the runtime reads only from Datagrok: run it on a stand.

import {execFileSync, execSync} from 'node:child_process';
import * as path from 'node:path';

const DEV = import.meta.dirname;
const RUNTIME = path.join(DEV, '..', 'dockerfiles', 'claude-runtime');
const MCP = path.join(DEV, '..', 'dockerfiles', 'mcp-server');
const NODE_PASS = /^(?:# |ℹ )pass (\d+)$/m;

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
  const out = execSync('npm test 2>&1', {cwd: RUNTIME, encoding: 'utf8'});
  return `(${out.match(NODE_PASS)?.[1] ?? '?'} tests)`;
});

// The MCP server is a separate image with its own build, so a runtime-only check would miss a
// broken op registry until a benchmark run failed on it 20 minutes in.
step('mcp-server', () => {
  execFileSync('npx', ['tsc'], {cwd: MCP, encoding: 'utf8', shell: true});
  const out = execFileSync('node', ['--test', 'test/'], {cwd: MCP, encoding: 'utf8'});
  return `(${out.match(NODE_PASS)?.[1] ?? '?'} tests)`;
});

// Cheap to check, expensive to get wrong: a malformed assert only surfaces ~20 minutes into a run.
step('benchmark suite', () => {
  const harness = path.join(DEV, 'harness');
  const out = execFileSync('node', ['lint-suite.mjs'], {cwd: harness, encoding: 'utf8'});
  return `(${out.match(/^suite: (\d+) prompts$/m)?.[1] ?? '?'} prompts)`;
});

console.log('\nall green');
