#!/usr/bin/env node
// The automatic audits of a grok-bdd run, from its JSON report: which audits ran, findings by rule, by
// feature (leaks split into scopes, timers and listeners), by scenario, and the most frequent ones; with
// BDD_MEMORY_LOG files, the heap and renderer of each worker's page from its first feature to its last.
// usage: node audit-summary.mjs <report.json> [memory.jsonl ...]
import {readFileSync} from 'node:fs';

const args = process.argv.slice(2);
const reportFile = args.find((a) => a.endsWith('.json'));
const memoryFiles = args.filter((a) => a.endsWith('.jsonl'));
if (!reportFile) {
  console.error('usage: node audit-summary.mjs <report.json> [memory.jsonl ...]');
  process.exit(2);
}

const FINDING = 'bdd-audit';
const RAN = 'bdd-audits';
const SCENARIO = / \(scenario "((?:[^"\\]|\\.)*)"\)$/;

const report = JSON.parse(readFileSync(reportFile, 'utf8'));

/** Every test of the report with its spec file and the audit annotations of its last result. */
function* tests(suite, file) {
  for (const spec of suite.specs ?? []) {
    for (const test of spec.tests ?? []) {
      const last = test.results?.[test.results.length - 1];
      const ours = (list) => (list ?? []).filter((a) => a.type === FINDING || a.type === RAN);
      const annotations = ours(last?.annotations).length > 0 ? ours(last.annotations) : ours(test.annotations);
      yield {file: spec.file ?? file, title: spec.title, status: test.status, annotations};
    }
  }
  for (const child of suite.suites ?? [])
    yield* tests(child, child.file ?? file);
}

const kindOf = (detail) => /^scope /.test(detail) ? 'scopes' : / timer still pending/.test(detail) ? 'timers' :
  / listener still attached/.test(detail) ? 'listeners' : 'other';
const count = (map, key, n = 1) => map.set(key, (map.get(key) ?? 0) + n);
const sorted = (map) => [...map.entries()].sort((a, b) => b[1] - a[1] || String(a[0]).localeCompare(String(b[0])));
const pad = (s, n) => String(s).padEnd(n);

const all = [...(report.suites ?? []).flatMap((s) => [...tests(s, s.file)])];
const byRule = new Map();
const byDetail = new Map();
const ran = {layout: 0, errors: 0, leaks: 0};
const features = new Map();
const scenarios = [];
let audited = 0;
let findings = 0;
for (const t of all) {
  const ranHere = new Set(t.annotations.filter((a) => a.type === RAN).flatMap((a) => a.description.replace(SCENARIO, '').split(', ')));
  if (ranHere.size > 0)
    audited++;
  for (const name of ranHere)
    ran[name] = (ran[name] ?? 0) + 1;
  const feature = t.file.replace(/\.test\.ts$/, '');
  const f = features.get(feature) ?? {tests: 0, layout: 0, errors: 0, scopes: 0, timers: 0, listeners: 0};
  features.set(feature, f);
  f.tests++;
  const perScenario = new Map();
  for (const a of t.annotations.filter((x) => x.type === FINDING)) {
    findings++;
    const scenario = SCENARIO.exec(a.description)?.[1] ?? t.title;
    const text = a.description.replace(SCENARIO, '');
    const colon = text.indexOf(': ');
    const rule = text.slice(0, colon);
    const detail = text.slice(colon + 2);
    count(byRule, rule);
    count(byDetail, `${rule}: ${detail.replace(/#\d+/g, '#n')}`);
    const s = perScenario.get(scenario) ?? {layout: 0, errors: 0, scopes: 0, timers: 0, listeners: 0};
    perScenario.set(scenario, s);
    const key = rule === 'leak' ? kindOf(detail) : rule === 'console-error' ? 'errors' : 'layout';
    if (key in s) {
      s[key]++;
      f[key]++;
    }
  }
  for (const [scenario, s] of perScenario)
    scenarios.push({feature, scenario: scenario === t.title ? scenario : `${t.title} › ${scenario}`, ...s});
}

const st = report.stats ?? {};
console.log(`Audits of ${reportFile}`);
console.log(`tests ${all.length}: expected ${st.expected ?? '?'}, unexpected ${st.unexpected ?? '?'}, skipped ${st.skipped ?? '?'}, ` +
  `flaky ${st.flaky ?? '?'}; ${((st.duration ?? 0) / 1000).toFixed(1)} s`);
console.log(`audited ${audited} test(s): layout ${ran.layout}, console errors ${ran.errors}, leaks ${ran.leaks}; findings ${findings}`);

console.log('\nFindings by rule');
for (const [rule, n] of sorted(byRule))
  console.log(`  ${pad(rule, 18)} ${n}`);
if (byRule.size === 0)
  console.log('  none');

console.log('\nBy feature (tests | layout | console errors | leaked scopes, timers, listeners)');
for (const [feature, f] of [...features.entries()].sort((a, b) => a[0].localeCompare(b[0]))) {
  if (f.layout + f.errors + f.scopes + f.timers + f.listeners === 0)
    continue;
  console.log(`  ${pad(feature, 48)} ${pad(f.tests, 3)} | ${pad(f.layout, 3)} | ${pad(f.errors, 3)} | ${f.scopes}, ${f.timers}, ${f.listeners}`);
}
const clean = [...features.values()].filter((f) => f.layout + f.errors + f.scopes + f.timers + f.listeners === 0).length;
console.log(`  (${clean} of ${features.size} feature files without findings)`);

if (scenarios.length > 0) {
  console.log('\nBy scenario (layout | console errors | leaked scopes, timers, listeners)');
  for (const s of scenarios)
    console.log(`  ${s.feature} › ${s.scenario}: ${s.layout} | ${s.errors} | ${s.scopes}, ${s.timers}, ${s.listeners}`);
}

if (byDetail.size > 0) {
  console.log('\nMost frequent (20)');
  for (const [detail, n] of sorted(byDetail).slice(0, 20))
    console.log(`  ${pad(`${n}×`, 5)} ${detail}`);
}

for (const file of memoryFiles) {
  const lines = readFileSync(file, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)).filter((l) => !l.error);
  console.log(`\nMemory (${file}): ${lines.length} feature end(s)`);
  const workers = new Map();
  for (const l of lines) {
    const list = workers.get(l.worker) ?? [];
    workers.set(l.worker, list);
    list.push(l);
  }
  const mb = (l, key) => l?.[key]?.heapUsedMb ?? '?';
  const renderer = (l, key) => l?.[key]?.renderer?.mb ?? '?';
  for (const [worker, list] of [...workers.entries()].sort()) {
    const first = list[0];
    const last = list[list.length - 1];
    const pages = list.filter((l) => l.onPage === 1).length;
    console.log(`  worker ${worker}: ${list.length} features on ${pages} page(s); JS heap ${mb(first, 'page')} → ${mb(last, 'page')} MB, ` +
      `after GC ${mb(first, 'afterGc')} → ${mb(last, 'afterGc')} MB; renderer ${renderer(first, 'processes')} → ` +
      `${renderer(last, 'processes')} MB (after GC ${renderer(first, 'processesAfterGc')} → ${renderer(last, 'processesAfterGc')} MB)`);
  }
}
