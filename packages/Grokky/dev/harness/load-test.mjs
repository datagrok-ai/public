#!/usr/bin/env node
// Minimal load harness — N independent virtual users against the claude-runtime container.
//
// Each user is its own WebSocket connection + sessionId (exactly how distinct real users hit
// the runtime — only same-session turns queue server-side). All users start together (small
// connect stagger), then each loops: pick a prompt from the pool → send → wait for the full
// answer → think → next. Every event is printed as it happens and appended to a JSONL log,
// including the full answer text, so a run is directly observable and diffable.
//
//   node load-test.mjs --users 10                       # 10 users, full Datagrok prompt
//   node load-test.mjs --users 100 --mode bash          # 100 users, cheap echo turns (~$0.07/turn)
//   node load-test.mjs --users 50 --turns 5 --think 0-0 --yes   # max sustained pressure
//
// Options:
//   --users N        concurrent users (default 10)
//   --turns N        turns per user (default 3)
//   --mode full|bash full = prompts from the pool through the real Datagrok prompt (default);
//                    bash = `echo` via systemPromptMode:'bash' — real CLI subprocess + API call
//                    per turn at a fraction of the cost, for finding infrastructure limits
//   --prompts FILE   prompt pool, one per line, # comments (default: files/benchmark/suite.yaml)
//   --think LO-HI    seconds between a user's turns (default 3-8; 0-0 = none)
//   --stagger S      seconds over which connections open (default 2 — near-simultaneous)
//   --timeout S      per-turn timeout (default 180)
//   --seed N         PRNG seed for prompt choice + think times — same seed, same run (default 1)
//   --url WS         runtime (default ws://localhost:5355/ws — `node dev/runtime.mjs up`)
//   --host ALIAS     ~/.grok/config.yaml server whose key is sent as apiKey (per-user sync)
//   --api-key K      explicit apiKey
//   --model M        haiku | sonnet | opus (default: runtime default)
//   --container NAME docker container to sample cpu/mem/pids from (default grokky-dev-runtime;
//                    --no-docker to skip)
//   --label NAME     log file name (default <mode>-<users>u-<timestamp>)
//   --yes            skip the cost confirmation for large full-mode runs

import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import * as fs from 'node:fs';
import * as path from 'node:path';
import yaml from 'yaml';
import {RuntimeDriver, DEFAULT_URL} from './drive.mjs';

const execFileP = promisify(execFile);
const HARNESS = import.meta.dirname;

function arg(name, fallback = undefined) {
  const i = process.argv.indexOf(`--${name}`);
  return i >= 0 && process.argv[i + 1] && !process.argv[i + 1].startsWith('--') ?
    process.argv[i + 1] : fallback;
}
const flag = (name) => process.argv.includes(`--${name}`);

const users = parseInt(arg('users', '10'), 10);
const turnsPerUser = parseInt(arg('turns', '3'), 10);
const mode = arg('mode', 'full');
const [thinkLo, thinkHi] = arg('think', '3-8').split('-').map(Number);
const staggerS = parseFloat(arg('stagger', '2'));
const timeoutMs = parseInt(arg('timeout', '180'), 10) * 1000;
const seed = parseInt(arg('seed', '1'), 10);
const url = arg('url', DEFAULT_URL);
const model = arg('model');
const container = flag('no-docker') ? null : arg('container', 'grokky-dev-runtime');
const label = arg('label',
  `${mode}-${users}u-${new Date().toISOString().slice(0, 16).replace(/[:T]/g, '-')}`);

let apiKey = arg('api-key');
const host = arg('host');
if (host && !apiKey) {
  const cfgPath = path.join(process.env.USERPROFILE ?? process.env.HOME, '.grok', 'config.yaml');
  const server = yaml.parse(fs.readFileSync(cfgPath, 'utf8')).servers?.[host];
  if (!server) throw new Error(`No server '${host}' in ${cfgPath}`);
  apiKey = server.key;
}

if (mode === 'full' && users * turnsPerUser > 40 && !flag('yes')) {
  console.error(`${users * turnsPerUser} full-prompt turns ≈ ` +
    `$${(users * turnsPerUser * 0.12).toFixed(0)} — pass --yes, or use --mode bash.`);
  process.exit(1);
}

function loadPrompts() {
  const file = arg('prompts');
  if (file) {
    const lines = fs.readFileSync(file, 'utf8').split('\n')
      .map((l) => l.trim()).filter((l) => l && !l.startsWith('#'));
    if (!lines.length) throw new Error(`no prompts in ${file}`);
    return lines;
  }
  const suite = path.join(HARNESS, '..', '..', 'files', 'benchmark', 'suite.yaml');
  return yaml.parse(fs.readFileSync(suite, 'utf8')).map((e) => e.prompt);
}

// mulberry32 — deterministic runs from --seed
function rng(s) {
  return () => {
    s |= 0; s = (s + 0x6d2b79f5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const percentile = (values, p) => {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  return sorted[Math.min(sorted.length - 1, Math.ceil((p / 100) * sorted.length) - 1)];
};
const sec = (ms) => ms == null ? '—' : (ms / 1000).toFixed(1) + 's';
const oneLine = (s, n) => {
  const t = String(s ?? '').replace(/\s+/g, ' ').trim();
  return t.length > n ? t.slice(0, n - 1) + '…' : t;
};

// ---------------------------------------------------------------------------
// event log — one line to the console, one JSON object to the log file
// ---------------------------------------------------------------------------

const OUT_DIR = path.join(HARNESS, 'results');
fs.mkdirSync(OUT_DIR, {recursive: true});
const logPath = path.join(OUT_DIR, `${label}.jsonl`);
const logStream = fs.createWriteStream(logPath, {flags: 'w'});

const startedAt = Date.now();
let inFlight = 0;
let peakInFlight = 0;
const results = [];      // one entry per finished turn
let connectFailures = 0;

function log(event, consoleLine) {
  logStream.write(JSON.stringify({tMs: Date.now() - startedAt, inFlight, ...event}) + '\n');
  if (consoleLine != null)
    console.log(`[${((Date.now() - startedAt) / 1000).toFixed(1).padStart(6)}s] ${consoleLine}`);
}

const uid = (i) => 'u' + String(i).padStart(String(users - 1).length, '0');

// ---------------------------------------------------------------------------
// one virtual user
// ---------------------------------------------------------------------------

let stopRequested = false;
process.on('SIGINT', () => {
  if (stopRequested) process.exit(130);
  stopRequested = true;
  console.log('\nstopping after in-flight turns… (Ctrl+C again to force-quit)');
});

async function virtualUser(i, prompts) {
  const rand = rng(seed * 1000 + i);
  await sleep((staggerS * 1000 * i) / users);
  if (stopRequested) return;

  const d = new RuntimeDriver({url, apiKey});
  try {
    await d.connect();
  } catch (e) {
    connectFailures++;
    log({user: i, ev: 'connect_error', error: e.message}, `${uid(i)} ✗ connect: ${e.message}`);
    return;
  }
  log({user: i, ev: 'connect'}, `${uid(i)} connected`);

  const sessionId = `load-${label}-${uid(i)}`;
  try {
    for (let turn = 0; turn < turnsPerUser && !stopRequested; turn++) {
      const prompt = mode === 'bash' ? `echo ${uid(i)}-t${turn}` :
        prompts[Math.floor(rand() * prompts.length)];
      inFlight++;
      peakInFlight = Math.max(peakInFlight, inFlight);
      log({user: i, ev: 'send', turn, prompt},
        `${uid(i)} → t${turn} "${oneLine(prompt, 80)}"   (${inFlight} in flight)`);
      let r;
      try {
        r = await d.turn(prompt, {
          sessionId, timeoutMs,
          ...(model ? {model} : {}),
          ...(mode === 'bash' ? {systemPromptMode: 'bash'} : {}),
        });
      } finally {
        inFlight--;
      }
      results.push({user: i, turn, ...r});
      log({user: i, ev: 'done', turn, ok: r.ok, ttftMs: r.ttftMs, totalMs: r.totalMs,
        queued: r.queued, error: r.error, tools: r.tools, costUsd: r.metrics?.costUsd ?? null,
        answer: r.content},
      r.ok ?
        `${uid(i)} ✓ t${turn} ${sec(r.totalMs)} (ttft ${sec(r.ttftMs)})` +
          `${r.queued ? ' [queued]' : ''}  "${oneLine(r.content, 100)}"` :
        `${uid(i)} ✗ t${turn} ${sec(r.totalMs)}  ${r.error}`);
      if (turn < turnsPerUser - 1)
        await sleep((thinkLo + rand() * (thinkHi - thinkLo)) * 1000);
    }
  } finally {
    d.close();
    log({user: i, ev: 'disconnect'});
  }
}

// ---------------------------------------------------------------------------
// status ticker + docker sampling
// ---------------------------------------------------------------------------

async function dockerSample() {
  const {stdout} = await execFileP('docker',
    ['stats', '--no-stream', '--format', '{{.CPUPerc}}|{{.MemUsage}}|{{.PIDs}}', container]);
  const [cpu, mem, pids] = stdout.trim().split('|');
  const m = /^([\d.]+)\s*(B|KiB|MiB|GiB)/.exec(mem);
  const memMb = m ? Math.round(parseFloat(m[1]) *
    {B: 1 / 2 ** 20, KiB: 1 / 1024, MiB: 1, GiB: 1024}[m[2]]) : null;
  return {cpuPct: parseFloat(cpu), memMb, pids: parseInt(pids, 10)};
}

function startTicker() {
  let stopped = false;
  let containerUp = true;
  (async () => {
    while (!stopped) {
      await sleep(5000);
      if (stopped) break;
      let dock = '';
      if (container) {
        try {
          const s = await dockerSample();
          dock = ` | docker cpu ${s.cpuPct}% mem ${s.memMb}MB pids ${s.pids}`;
          log({ev: 'docker', ...s});
          if (!containerUp) { log({ev: 'container_up'}, '⚠ container back up'); containerUp = true; }
        } catch (e) {
          dock = ' | docker: unreachable';
          if (containerUp) {
            log({ev: 'container_down', error: oneLine(e.message, 200)}, `⚠ container down/unreachable`);
            containerUp = false;
          }
        }
      }
      const ok = results.filter((r) => r.ok).length;
      log({ev: 'status', ok, failed: results.length - ok + connectFailures, peakInFlight},
        `status: ${inFlight} in flight, ${ok} ok, ${results.length - ok + connectFailures} failed` +
        ` (peak concurrency ${peakInFlight})${dock}`);
    }
  })();
  return () => { stopped = true; };
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

const prompts = mode === 'bash' ? null : loadPrompts();
console.log(`${users} users × ${turnsPerUser} turns, mode=${mode}` +
  `${prompts ? ` (${prompts.length} prompts in pool)` : ''}, think ${thinkLo}-${thinkHi}s, ` +
  `seed ${seed} → ${url}${apiKey ? ' (with apiKey)' : ''}\nlog: ${logPath}`);
log({ev: 'run_start', users, turnsPerUser, mode, thinkLo, thinkHi, staggerS, seed, url,
  model: model ?? null, apiKey: !!apiKey, promptPool: prompts?.length ?? null,
  startedIso: new Date().toISOString()});

const stopTicker = startTicker();
await Promise.all(Array.from({length: users}, (_, i) => virtualUser(i, prompts)));
stopTicker();

// runtime log signals — correlate client-side timeouts with server-side kills
let watchdogKills = 0;
if (container) {
  try {
    const since = new Date(startedAt).toISOString();
    const {stdout, stderr} = await execFileP('docker', ['logs', '--since', since, container],
      {maxBuffer: 32 * 1024 * 1024});
    const lines = (stdout + '\n' + stderr).split('\n')
      .filter((l) => /watchdog|reaper|stray|error|failed/i.test(l));
    watchdogKills = lines.filter((l) => l.includes('watchdog[')).length;
    log({ev: 'runtime_logs', watchdogKills, lines: lines.slice(0, 300)});
  } catch { /* container gone — already reported by the ticker */ }
}

const ok = results.filter((r) => r.ok);
const failed = results.filter((r) => !r.ok);
const ttfts = ok.map((r) => r.ttftMs).filter((x) => x != null);
const totals = ok.map((r) => r.totalMs);
const cost = results.reduce((a, r) => a + (r.metrics?.costUsd ?? 0), 0);
const summary = {
  ev: 'summary',
  users, connectFailures, turns: results.length, ok: ok.length, failed: failed.length,
  queued: results.filter((r) => r.queued).length, peakInFlight,
  ttftP50Ms: percentile(ttfts, 50), ttftP90Ms: percentile(ttfts, 90),
  totalP50Ms: percentile(totals, 50), totalP90Ms: percentile(totals, 90),
  totalMaxMs: percentile(totals, 100), costUsd: +cost.toFixed(2), watchdogKills,
  errors: Object.fromEntries(failed.reduce((m, r) =>
    m.set(r.error ?? 'unknown', (m.get(r.error ?? 'unknown') ?? 0) + 1), new Map())),
};
log(summary);
await new Promise((r) => logStream.end(r));

console.log(`\n— summary ————————————————————————————————`);
console.log(`turns: ${ok.length}/${results.length} ok` +
  `${connectFailures ? `, ${connectFailures} users failed to connect` : ''}` +
  `${summary.queued ? `, ${summary.queued} queued` : ''}`);
console.log(`peak concurrency: ${peakInFlight}/${users} turns in flight`);
console.log(`ttft p50/p90: ${sec(summary.ttftP50Ms)}/${sec(summary.ttftP90Ms)}   ` +
  `total p50/p90/max: ${sec(summary.totalP50Ms)}/${sec(summary.totalP90Ms)}/${sec(summary.totalMaxMs)}` +
  `${cost ? `   cost $${cost.toFixed(2)}` : ''}`);
for (const [e, n] of Object.entries(summary.errors))
  console.log(`✗ ${n}× ${e}`);
if (watchdogKills)
  console.log(`⚠ ${watchdogKills} watchdog kill(s) in runtime logs — see 'runtime_logs' in the log`);
console.log(`log: ${logPath}`);
process.exit(0);
