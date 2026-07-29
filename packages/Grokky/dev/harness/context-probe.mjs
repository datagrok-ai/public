#!/usr/bin/env node
// Context-diet instrument: attributes the per-turn prompt prefix to its components.
//
// A turn's prefix is everything the model re-reads before it can emit a token: the SDK's own
// system prompt, our Datagrok prompt with its inlined skills, the plugin's skill descriptions,
// and every declared tool. `benchmark-sonnet-clean` showed a floor of ~54k tokens on a prompt
// that used no tools at all, so the prefix — not the conversation — is what TTFT is paying for.
// This measures each component by ablation instead of guessing from file sizes.
//
//   node context-probe.mjs                 # full ablation, writes files/benchmarks/context-probe.md
//   node context-probe.mjs --label after   # name the run so before/after can be diffed
//
// Reads total prefix as input + cacheRead + cacheCreation, which is cache-state independent:
// whichever way the cache breaks, those three sum to what the model actually saw.

import * as fs from 'node:fs';
import * as path from 'node:path';
import {RuntimeDriver} from './drive.mjs';

const arg = (name, dflt) => {
  const i = process.argv.indexOf(`--${name}`);
  return i >= 0 ? process.argv[i + 1] : dflt;
};

const label = arg('label', 'context-probe');
const mcpUrl = arg('mcp', 'http://localhost:32850/mcp');
const OUT_DIR = path.join(import.meta.dirname, '..', '..', 'files', 'benchmarks');

// The three view meta-tools, copied from src/ai/view-tools.ts. Duplicated rather than imported
// because that file is browser TS (imports datagrok-api); the probe only needs their wire shape.
const VIEW_TOOLS = [
  {name: 'list_view_functions', description: 'Search the functions applicable to the current view.',
    inputSchema: {type: 'object', properties: {query: {type: 'string', description: 'Search words.'}}}},
  {name: 'get_view_function_result', description: 'Invoke a READ-ONLY function of the current view.',
    inputSchema: {type: 'object', properties: {name: {type: 'string'}, parameters: {type: 'object'}}, required: ['name']}},
  {name: 'call_view_function', description: 'Invoke a state-changing function of the current view.',
    inputSchema: {type: 'object', properties: {name: {type: 'string'}, parameters: {type: 'object'}}, required: ['name']}},
];

// Each config isolates one contributor. `base` is the floor everything else is measured against.
const CONFIGS = [
  {key: 'base', title: 'SDK floor (no prompt, no plugin, no MCP, no view tools)',
    opts: {systemPromptMode: 'none'}, mcp: false},
  {key: 'mcp', title: '+ datagrok MCP server (32 tool defs)',
    opts: {systemPromptMode: 'none'}, mcp: true},
  {key: 'view', title: '+ view meta-tools (3 defs)',
    opts: {systemPromptMode: 'none', clientTools: VIEW_TOOLS}, mcp: false},
  {key: 'prompt', title: '+ Datagrok system prompt, inlined skills, plugin',
    opts: {}, mcp: false},
  {key: 'full', title: 'production (everything on)',
    opts: {clientTools: VIEW_TOOLS}, mcp: true},
];

const PROMPT = 'Reply with exactly: PROBE';

// SDK usage on the `result` message is CUMULATIVE over every API call the turn made, so a turn
// that thought and then answered counts its prefix twice. Dividing by numTurns recovers the
// per-call prefix — the thing TTFT actually pays for. Without this, enabling thinking looks
// like a context-size regression when it is really just a second call.
function prefixOf(m) {
  if (!m) return null;
  const n = (v) => (typeof v === 'number' ? v : 0);
  const total = n(m.inputTokens) + n(m.cacheReadTokens) + n(m.cacheCreationTokens);
  return {total, calls: m.numTurns || 1, perCall: Math.round(total / (m.numTurns || 1))};
}

async function measure(cfg) {
  const driver = await new RuntimeDriver(cfg.mcp ? {mcpServerUrl: mcpUrl} : {}).connect();
  try {
    // Two turns: the first may pay cache-creation on a cold prefix and the SDK occasionally
    // spends an extra round-trip warming up. The second is the steady state a user actually hits.
    let last = null;
    for (let i = 0; i < 2; i++)
      last = await driver.turn(PROMPT, {...cfg.opts, timeoutMs: 120000});
    const p = prefixOf(last.metrics);
    return {
      prefix: p?.perCall ?? null, total: p?.total ?? null, calls: p?.calls ?? null,
      ttftMs: last.ttftMs, totalMs: last.totalMs,
      error: last.error ?? (last.metrics ? null : 'no metrics — runtime too old?'),
    };
  } finally {
    driver.close();
  }
}

const only = arg('only', null)?.split(',').map((s) => s.trim());

const rows = [];
for (const cfg of CONFIGS.filter((c) => !only || only.includes(c.key))) {
  process.stdout.write(`${cfg.key.padEnd(7)} … `);
  try {
    const r = await measure(cfg);
    rows.push({...cfg, ...r});
    console.log(r.error ? `FAILED (${r.error})` :
      `${r.prefix.toLocaleString()} tok/call  (${r.total.toLocaleString()} over ${r.calls} call${r.calls === 1 ? '' : 's'})  ttft ${r.ttftMs}ms`);
  } catch (e) {
    rows.push({...cfg, prefix: null, error: e.message});
    console.log(`FAILED (${e.message})`);
  }
}

const base = rows.find((r) => r.key === 'base')?.prefix;
const full = rows.find((r) => r.key === 'full')?.prefix;
const delta = (r) => (r.prefix != null && base != null && r.key !== 'base' ? r.prefix - base : null);

const lines = [
  `# Context probe — ${label}`, '',
  `Run ${new Date().toISOString()}. Prefix = input + cacheRead + cacheCreation, so it is what the`,
  'model re-read before emitting a token, independent of how the cache happened to break.', '',
  '| Config | Prefix (tok/call) | vs floor | Calls | TTFT | What it adds |',
  '|---|---:|---:|---:|---:|---|',
];
for (const r of rows) {
  const d = delta(r);
  lines.push(`| \`${r.key}\` | ${r.prefix?.toLocaleString() ?? '—'} | ${d != null ? `+${d.toLocaleString()}` : '—'} | ` +
    `${r.calls ?? '—'} | ${r.ttftMs != null ? r.ttftMs + 'ms' : '—'} | ${r.title} |`);
}
if (base != null && full != null) {
  lines.push('', `**Floor ${base.toLocaleString()} tok · production ${full.toLocaleString()} tok** — ` +
    `Datagrok adds ${(full - base).toLocaleString()} tok (${Math.round(((full - base) / full) * 100)}% of the prefix) on every turn.`);
}
lines.push('', 'Config `full` is not the sum of the parts: components share cache blocks and the SDK',
  'reorders some of them, so read the deltas as attribution, not as an exact budget.', '');

fs.mkdirSync(OUT_DIR, {recursive: true});
const outPath = path.join(OUT_DIR, `context-probe-${label}.md`);
fs.writeFileSync(outPath, lines.join('\n'));
fs.writeFileSync(path.join(OUT_DIR, `context-probe-${label}.json`),
  JSON.stringify({label, when: new Date().toISOString(), rows}, null, 2));

console.log('\n' + lines.join('\n'));
console.log(`\nwritten: ${outPath}`);
process.exit(rows.some((r) => r.error) ? 1 : 0);
