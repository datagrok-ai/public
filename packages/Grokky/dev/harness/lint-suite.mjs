#!/usr/bin/env node
// Pre-flight for files/benchmark/suite.yaml.
//
// A benchmark run costs ~20 minutes and real subscription quota, so a malformed YAML entry or a
// typo'd assert must be caught before it starts — not discovered from a report full of
// "assert threw". Compiles every assert exactly the way benchmark.ts does, and flags assertions
// that are too weak to carry the accuracy axis.

import * as fs from 'node:fs';
import * as path from 'node:path';
import yaml from 'yaml'; // CommonJS — no named exports
const {parse} = yaml;

const SUITE = path.join(import.meta.dirname, '..', '..', 'files', 'benchmark', 'suite.yaml');
const CATEGORIES = ['help', 'visualization', 'analysis', 'codegen', 'multitool', 'query'];
const DIFFICULTIES = ['trivial', 'standard', 'hard'];
// Mirrors the identifiers benchmark.ts injects into the assert function.
const SCOPE = ['grok', 'DG', 'view', 't', 'before', 'opened', 'openedViews', 'tools'];

const errors = [];
const warnings = [];

let suite;
try {
  suite = parse(fs.readFileSync(SUITE, 'utf8'));
} catch (e) {
  console.error(`suite.yaml is not valid YAML:\n  ${e.message}`);
  process.exit(1);
}
if (!Array.isArray(suite)) {
  console.error('suite.yaml must be a list of prompts');
  process.exit(1);
}

const at = (i, p) => `#${i + 1} "${String(p.prompt ?? '(no prompt)').slice(0, 52)}"`;

suite.forEach((p, i) => {
  if (!p || typeof p !== 'object') return errors.push(`${at(i, {})} is not a mapping`);
  if (!p.prompt) errors.push(`${at(i, p)} has no prompt`);
  if (!CATEGORIES.includes(p.category))
    errors.push(`${at(i, p)} has unknown category '${p.category}' (expected one of ${CATEGORIES.join(', ')})`);
  if (p.difficulty && !DIFFICULTIES.includes(p.difficulty))
    errors.push(`${at(i, p)} has unknown difficulty '${p.difficulty}'`);
  if (!p.assert && !p.rubric)
    warnings.push(`${at(i, p)} is latency-only (no assert, no rubric) — it cannot score`);

  if (p.assert) {
    // Exactly the wrapper benchmark.ts builds, so a syntax error here is a syntax error there.
    try {
      new Function(...SCOPE, `return (async () => (${p.assert}))();`);
    } catch (e) {
      errors.push(`${at(i, p)} assert does not compile: ${e.message}`);
    }
    // Missing values in a numeric column are NOT null: toList() yields `undefined` and get()
    // yields a 2.68e-34 sentinel, so `v !== null` counts missing rows as present and
    // `v === null` never matches one. This silently shifted two counts by exactly the number of
    // missing rows and failed correct model answers. isNone(i) is the only reliable test.
    if (/(!==|===)\s*null/.test(p.assert) && /toList\(\)|\.get\(/.test(p.assert))
      warnings.push(`${at(i, p)} tests a numeric value against null — missing numerics are ` +
        `undefined (toList) / 2.68e-34 (get). Use col.isNone(i) instead`);

    // An existence-only viewer check passes on a half-correct answer and flattens accuracy.
    // FILTERS is exempt — it is a panel, not a column-bound viewer.
    if (/\.type === DG\.VIEWER\./.test(p.assert) && !/props\./.test(p.assert) &&
        !/DG\.VIEWER\.FILTERS/.test(p.assert))
      warnings.push(`${at(i, p)} asserts a viewer exists but binds no column — consider checking props.*ColumnName`);
    if (/^t\.columns\.length > before\.cols\.length$/.test(p.assert.trim()))
      warnings.push(`${at(i, p)} only counts columns — any added column passes, including a wrong one`);
    for (const id of p.assert.matchAll(/\b(await\s+)?([a-zA-Z_$][\w$]*)\s*\(/g)) {
      // Cheap catch for a scope identifier that does not exist (e.g. `openedView` singular).
      const name = id[2];
      if (/^(openedView|toolNames|newViews)$/.test(name))
        errors.push(`${at(i, p)} references '${name}', which is not in the assert scope (${SCOPE.join(', ')})`);
    }
  }

  // A rubric that hardcodes an expected value is only as good as that value. One of these said
  // 1274 where the real answer was 878, which failed a correct model answer and read as a model
  // regression. Verify with dev/harness/introspect.mjs before trusting one.
  if (p.rubric && /\b\d{2,}\b/.test(p.rubric))
    warnings.push(`${at(i, p)} rubric hardcodes ${JSON.stringify(p.rubric.match(/\b\d{2,}\b/g))} — ` +
      `verify against the live table (introspect.mjs), a stale constant fails a correct answer`);

  if (p.table && /^SPGI$/i.test(p.table))
    errors.push(`${at(i, p)} uses table '${p.table}' — SPGI lives at 'chem/SPGI' under DemoFiles`);
});

const byKey = (key) => {
  const counts = {};
  for (const p of suite) counts[p?.[key] ?? (key === 'difficulty' ? 'standard' : '?')] =
    (counts[p?.[key] ?? (key === 'difficulty' ? 'standard' : '?')] ?? 0) + 1;
  return Object.entries(counts).map(([k, v]) => `${k} ${v}`).join(' · ');
};

console.log(`suite: ${suite.length} prompts`);
console.log(`  by difficulty : ${byKey('difficulty')}`);
console.log(`  by category   : ${byKey('category')}`);
console.log(`  scored        : ${suite.filter((p) => p?.assert).length} assert · ${suite.filter((p) => p?.rubric).length} rubric`);
console.log(`  tables        : ${[...new Set(suite.map((p) => p?.table).filter(Boolean))].join(', ') || '(none)'}`);

for (const w of warnings) console.log(`  warn: ${w}`);
for (const e of errors) console.error(`  ERROR: ${e}`);
if (errors.length) {
  console.error(`\n${errors.length} error(s) — fix before running the benchmark.`);
  process.exit(1);
}
console.log('\nsuite ok');
