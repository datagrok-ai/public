import {spawn} from 'node:child_process';
import type {ChildProcessWithoutNullStreams} from 'node:child_process';
import {randomUUID} from 'node:crypto';
import * as fs from 'node:fs';
import * as path from 'node:path';
import * as readline from 'node:readline';
import type {AgentCommand, AgentProcessEvent} from '../types';
import {AGENT_DIR, AGENT_PYTHON, HERMES_HOMES, PLUGIN_DIR, PROVIDER_ENV} from '../constants';
import {IDLE, killTree, raceIdle, startOrphanReaper} from '../watchdog';

// One agent process (agent/grokky_agent) per Grokky session, each with its own HERMES_HOME and process group.
// An idle one holds about 250 MB RSS: the default 12 processes, spares included, fit the container's 8 GB.
const envNumber = (name: string, fallback: number) => process.env[name] ? Number(process.env[name]) : fallback;
const MAX_PROCESSES = envNumber('HERMES_MAX_PROCESSES', 12);
const IDLE_TTL_MS = envNumber('HERMES_IDLE_TTL_MS', 15 * 60_000);
const WARM_SPARES = Math.max(0, envNumber('HERMES_WARM_SPARES', 2));
const MAX_SESSIONS = envNumber('HERMES_MAX_SESSIONS', 200);
const COMPRESSION_THRESHOLD_TOKENS = envNumber('HERMES_COMPRESSION_THRESHOLD_TOKENS', 160_000);
const readyTimeoutMs = () => envNumber('HERMES_READY_TIMEOUT_MS', 30_000);
const PARKED_DIR = path.join(HERMES_HOMES, 'parked');
// Hermes' session store (sqlite, WAL); the rest of a home is rebuilt on reload.
const SESSION_STATE = /^state\.db(-wal|-shm)?$/;
const STDERR_TAIL_LINES = 40;
// Unexpected exits in a row after which a session's store is dropped: the chat then recovers from the transcript.
const MAX_CRASHES = 2;

// The configured key and anything key-shaped, wherever an agent's output reaches the log or the browser.
const redact = (text: string) => (PROVIDER_ENV.apiKey ? text.split(PROVIDER_ENV.apiKey).join('<redacted>') : text)
  .replace(/\bsk-[\w.…-]+/g, 'sk-<redacted>');

const CHECK_TIMEOUT_MS = 120_000;

/** The startup connection test's verdict: null when the model answered, else why not; undefined while it runs. */
export let providerCheck: {error: string | null} | undefined;
let checking: ChildProcessWithoutNullStreams | undefined;

export interface AgentProcess {
  proc: ChildProcessWithoutNullStreams;
  home: string;
  ready: Promise<void>;
  /** `ready` has resolved: until then the process is neither idle nor evictable. */
  started: boolean;
  disposed: boolean;
  sessionKey?: string;
  /** The Hermes session the process holds; a parked session carries it to the process that reloads it. */
  hermesSession?: string;
  /** The parked session a reloading process has not opened yet; sent with each turn until the agent has. */
  restore?: string;
  turns: number;
  turn: {events: (ev: AgentProcessEvent) => void; settled: Promise<void>} | null;
  lastUsed: number;
  stderrTail: string[];
  /** The home is a stored session's and outlives this process (a reload that failed, or a process that died holding
   * its session). */
  keepHome?: boolean;
  /** Unexpected exits of the session's processes since a turn last completed. */
  crashes: number;
  beforeForget?: () => void;
  send(cmd: AgentCommand): void;
}

interface ParkedSession {
  home: string;
  hermesSession: string;
  crashes?: number;
  /** Resolves once the evicted process has exited and its session store is in `home`. */
  stored: Promise<void>;
}

/** The stored history of a session is gone: the turn resets the chat and the browser replays it. */
export class SessionLost extends Error {}

const bySession = new Map<string, AgentProcess>();
const all = new Set<AgentProcess>();
const spares: AgentProcess[] = [];
const parked = new Map<string, ParkedSession>();
// Spares still starting; a session that finds no ready spare takes one of these over rather than paying for a
// second cold start next to it.
const starting: {proc: Promise<AgentProcess>; claimed: boolean}[] = [];

export const processForSession = (key: string) => bySession.get(key);
export const isParked = (key: string) => parked.has(key);
export const hasSession = (key: string) => holdsSession(bySession.get(key)) || parked.has(key);
// A process restoring a stored session holds it before its first turn has run.
const holdsSession = (p?: AgentProcess) => !!p && (p.turns > 0 || p.restore !== undefined);
export const hasExited = (p: AgentProcess) => p.proc.exitCode !== null || p.proc.signalCode !== null;

export function ownedPids(): Set<number> {
  return new Set([...all, ...(checking ? [{proc: checking}] : [])].flatMap((p) => p.proc.pid ?? []));
}

async function spawnProcess(reload?: ParkedSession): Promise<AgentProcess> {
  // Only a reload waits: a spare is in `all` the moment its spawn begins, which topUpSpares counts on.
  if (reload)
    await reload.stored;
  if (reload && !fs.existsSync(reload.home))
    throw new SessionLost(`the stored home ${reload.home} is gone`);
  const spawnedAt = Date.now();
  const home = reload?.home ?? path.join(HERMES_HOMES, randomUUID());
  // Another session's agent (same uid) can write into a stored home; only the session store is kept.
  for (const f of reload ? fs.readdirSync(home).filter((n) => !SESSION_STATE.test(n)) : [])
    fs.rmSync(path.join(home, f), {recursive: true, force: true});
  fs.mkdirSync(home, {recursive: true, mode: 0o700});
  const proc = spawnAgent(home, !reload);
  const tag = path.basename(home).slice(0, 8);
  const p: AgentProcess = {
    proc, home, ready: Promise.resolve(), started: false, disposed: false, turns: 0, turn: null,
    lastUsed: Date.now(), stderrTail: [], restore: reload?.hermesSession, hermesSession: reload?.hermesSession,
    crashes: reload?.crashes ?? 0,
    send: (cmd) => {
      if (proc.stdin.writable)
        proc.stdin.write(JSON.stringify(cmd) + '\n');
    },
  };
  all.add(p);
  // Whole lines: a key split across two reads would escape the redaction.
  readline.createInterface({input: proc.stderr}).on('line', (raw) => {
    if (!raw)
      return;
    const line = redact(raw);
    p.stderrTail.push(line);
    if (p.stderrTail.length > STDERR_TAIL_LINES)
      p.stderrTail.shift();
    if (process.env['HERMES_LOG_STDERR'])
      console.log(`agent[${tag}] ${line}`);
  });
  let started!: () => void;
  let failed!: (e: Error) => void;
  const readyMs = readyTimeoutMs();
  p.ready = raceIdle(new Promise<void>((resolve, reject) => {
    started = resolve;
    failed = reject;
  }), readyMs).then((r) => {
    if (r === IDLE)
      throw new Error(`the agent startup did not finish in ${readyMs / 1000} s`);
    p.started = true;
  });
  readline.createInterface({input: proc.stdout}).on('line', (line) => {
    let ev: AgentProcessEvent;
    try {
      ev = JSON.parse(line);
    } catch {
      return console.warn(`agent[${tag}]: a stdout line is not an event`);
    }
    if (ev.type === 'error' || ev.type === 'fatal')
      ev.message = redact(ev.message);
    if (ev.type === 'ready')
      started();
    else if (ev.type === 'fatal')
      failed(new Error(ev.message));
    else if (ev.type === 'session') {
      if (p.restore !== undefined)
        console.log(`agent: reloaded session ${p.sessionKey} in ${Date.now() - spawnedAt} ms`);
      p.hermesSession = ev.id;
      p.restore = undefined;
    }
    if (ev.type === 'stopped')
      p.crashes = 0;
    p.turn?.events(ev);
  });
  proc.on('error', failed);
  proc.on('exit', (code, signal) => {
    const last = p.stderrTail.at(-1);
    failed(new Error(`the agent exited (${signal ?? code})${last ? `: ${last}` : ''}`));
    p.beforeForget?.();
    if (!p.disposed) {
      console.warn(`agent[${tag}]: exited unexpectedly (${signal ?? code}); last stderr:\n` +
        p.stderrTail.slice(-10).join('\n'));
      const stored = p.restore ?? p.hermesSession;
      if (p.sessionKey && stored && holdsSession(p))
        keepStored(p, p.sessionKey, stored);
    }
    forget(p);
  });
  return p;
}

/** An agent process in its own home, sent its init command. The credentials go over stdin, never the environment,
 * which every tool the agent runs inherits. The cwd is root-owned, so a stray relative path never lands in runtime
 * state. */
function spawnAgent(home: string, warm: boolean, args: string[] = []): ChildProcessWithoutNullStreams {
  const proc = spawn(AGENT_PYTHON, ['-m', 'grokky_agent', ...args], {
    cwd: '/', detached: true,
    env: {PATH: process.env['PATH'] ?? '/usr/bin:/bin', HOME: home, HERMES_HOME: home, TERM: 'dumb', LANG: 'C.UTF-8',
      PYTHONUNBUFFERED: '1', PYTHONPATH: AGENT_DIR},
  });
  const settings = {homes: HERMES_HOMES, skillsDir: path.join(PLUGIN_DIR, 'skills'),
    compressionThresholdTokens: COMPRESSION_THRESHOLD_TOKENS, warm};
  const init: AgentCommand = {type: 'init', credentials: PROVIDER_ENV, settings};
  proc.stdin.on('error', () => {});
  proc.stdin.write(JSON.stringify(init) + '\n');
  return proc;
}

/** The connection test: a `--check` agent makes one short call to the configured model and exits. Its verdict is
 * logged and reported by /health; turns do not wait for it. */
function checkProvider(): void {
  const home = path.join(HERMES_HOMES, `check-${randomUUID()}`);
  fs.mkdirSync(home, {recursive: true, mode: 0o700});
  const proc = checking = spawnAgent(home, false, ['--check']);
  let verdict: string | null | undefined;
  let lastStderr = '';
  readline.createInterface({input: proc.stderr}).on('line', (line) => lastStderr = line || lastStderr);
  readline.createInterface({input: proc.stdout}).on('line', (line) => {
    try {
      const ev: AgentProcessEvent = JSON.parse(line);
      if (ev.type === 'checked')
        verdict = ev.error && redact(ev.error);
      else if (ev.type === 'fatal')
        verdict = redact(ev.message);
    }
    catch {}
  });
  const timer = setTimeout(() => proc.pid && killTree(proc.pid), CHECK_TIMEOUT_MS);
  proc.on('error', (e) => verdict ??= e.message);
  proc.once('close', (code, signal) => {
    clearTimeout(timer);
    checking = undefined;
    fs.rmSync(home, {recursive: true, force: true});
    providerCheck = {error: verdict !== undefined ? verdict : signal === 'SIGKILL' ?
      `the check did not finish in ${CHECK_TIMEOUT_MS / 1000} s` :
      redact(`the check exited (${signal ?? code})${lastStderr ? `: ${lastStderr}` : ''}`)};
    const what = `${PROVIDER_ENV['provider'] ?? 'Anthropic'} ${PROVIDER_ENV['model'] ?? ''}`.trim();
    console.log(providerCheck.error ? `provider check: ${what} failed: ${providerCheck.error}` :
      `provider check: ${what} answered`);
  });
}

/** A process that died holding a session keeps its home, and the next turn reloads the store there, unless the
 * session's processes keep dying. */
function keepStored(p: AgentProcess, key: string, hermesSession: string): void {
  const crashes = p.crashes + 1;
  if (crashes >= MAX_CRASHES)
    return console.warn(`agent: dropped the stored state of session ${key} after ${crashes} unexpected exits`);
  p.keepHome = true;
  parked.set(key, {home: p.home, hermesSession, stored: Promise.resolve(), crashes});
}

function forget(p: AgentProcess, keepHome = false): void {
  all.delete(p);
  if (p.sessionKey && bySession.get(p.sessionKey) === p)
    bySession.delete(p.sessionKey);
  const i = spares.indexOf(p);
  if (i >= 0)
    spares.splice(i, 1);
  if (!keepHome && !p.keepHome)
    fs.rmSync(p.home, {recursive: true, force: true, maxRetries: 3});
}

/** Kills the agent and every process under it, tool commands in their own sessions included. */
export function dispose(p: AgentProcess, reason: string): void {
  if (p.disposed)
    return;
  p.disposed = true;
  // An exited child's pid may already belong to another process.
  const failed = p.proc.pid && !hasExited(p) ? killTree(p.proc.pid) : [];
  console.log(`agent[${path.basename(p.home).slice(0, 8)}]: dispose (${reason})` +
    (failed.length ? `; could not signal ${failed.join(', ')}` : ''));
  forget(p);
}

function topUpSpares(): void {
  const unclaimed = () => starting.filter((e) => !e.claimed).length;
  while (spares.length + unclaimed() < WARM_SPARES && all.size < MAX_PROCESSES) {
    const entry = {proc: spawnProcess(), claimed: false};
    starting.push(entry);
    entry.proc.then(async (s) => {
      try {
        await s.ready;
        if (!entry.claimed)
          spares.push(s);
      } catch (e: any) {
        console.warn('agent: warm spare failed to start:', e.message);
        dispose(s, 'spare failed');
      }
    }).catch((e) => console.warn('agent: warm spare spawn failed:', e.message))
      .finally(() => starting.splice(starting.indexOf(entry), 1));
  }
}

async function takeSpare(): Promise<AgentProcess | undefined> {
  const ready = spares.shift();
  if (ready)
    return ready;
  const entry = starting.find((e) => !e.claimed);
  if (!entry)
    return undefined;
  entry.claimed = true;
  return entry.proc.catch(() => undefined);
}

/** Ends the session's process but keeps its Hermes session store, so the next turn reloads it.
 * The store moves once the process has exited, so no write of the dying agent is lost. */
function evict(p: AgentProcess, why: string): void {
  console.log(`agent: evicting session ${p.sessionKey} (${why}); ${bySession.size - 1} session(s) remain live`);
  const key = p.sessionKey;
  if (!key || !p.hermesSession || !holdsSession(p) || p.disposed)
    return dispose(p, why);
  p.disposed = true;
  forget(p, true);
  const home = path.join(PARKED_DIR, randomUUID());
  fs.mkdirSync(home, {recursive: true, mode: 0o700});
  let stored!: () => void;
  const parkedSession: ParkedSession = {home, hermesSession: p.hermesSession,
    stored: new Promise<void>((resolve) => { stored = resolve; })};
  p.beforeForget = () => {
    p.beforeForget = undefined;
    for (const f of fs.existsSync(p.home) ? fs.readdirSync(p.home).filter((n) => SESSION_STATE.test(n)) : [])
      fs.renameSync(path.join(p.home, f), path.join(home, f));
    fs.rmSync(p.home, {recursive: true, force: true, maxRetries: 3});
    stored();
  };
  // No exit is coming for a process that has exited or could not be signalled: its store moves now.
  if (hasExited(p) || !p.proc.pid || killTree(p.proc.pid).includes(p.proc.pid))
    p.beforeForget();
  parked.set(key, parkedSession);
  trimParked();
}

function trimParked(): void {
  for (const [key, s] of parked) {
    if (parked.size + bySession.size <= MAX_SESSIONS)
      return;
    parked.delete(key);
    fs.rmSync(s.home, {recursive: true, force: true, maxRetries: 3});
    console.log(`agent: dropped the stored state of session ${key} (cap ${MAX_SESSIONS} sessions)`);
  }
}

/** Ends the session's process and keeps its history: the next turn reloads it. */
export function parkSession(key: string, reason: string): void {
  const p = bySession.get(key);
  if (p)
    evict(p, reason);
}

export function dropSession(key: string, reason: string): void {
  const p = bySession.get(key);
  if (p)
    dispose(p, reason);
  const s = parked.get(key);
  if (s) {
    parked.delete(key);
    fs.rmSync(s.home, {recursive: true, force: true, maxRetries: 3});
  }
}

/** Back under the process cap (a busy session is never evicted, so the pool can run over it). */
function evictForRoom(): void {
  while (all.size >= MAX_PROCESSES) {
    const spare = spares.shift();
    if (spare) {
      dispose(spare, 'room for a session');
      continue;
    }
    const idle = [...bySession.values()].filter((p) => p.started && !p.turn)
      .sort((a, b) => a.lastUsed - b.lastUsed)[0];
    if (!idle)
      return;
    evict(idle, `LRU, cap ${MAX_PROCESSES}`);
  }
}

export function reapIdle(now = Date.now()): void {
  for (const p of [...bySession.values()]) {
    if (p.started && !p.turn && now - p.lastUsed > IDLE_TTL_MS)
      evict(p, `idle ${Math.round((now - p.lastUsed) / 60_000)} min`);
  }
}

const startError = (e: any) => new Error(`The assistant could not start: ${e?.message ?? e}`);

export async function acquire(sessionKey: string): Promise<AgentProcess> {
  const existing = bySession.get(sessionKey);
  if (existing) {
    await existing.ready.catch((e) => {
      throw startError(e);
    });
    existing.lastUsed = Date.now();
    return existing;
  }
  // Taken before eviction makes room, which could otherwise trim this very session.
  const stored = parked.get(sessionKey);
  parked.delete(sessionKey);
  let p: AgentProcess | undefined;
  try {
    // A warm spare needs no room; a process started for this session does.
    p = stored ? undefined : await takeSpare();
    if (!p) {
      evictForRoom();
      p = await spawnProcess(stored);
    }
  } catch (e) {
    if (stored && !(e instanceof SessionLost))
      parked.set(sessionKey, stored);
    throw e;
  }
  p.sessionKey = sessionKey;
  p.lastUsed = Date.now();
  bySession.set(sessionKey, p);
  trimParked();
  topUpSpares();
  try {
    await p.ready;
  } catch (e: any) {
    // A stored session stays stored: a transient failure must not cost its history. The turn errors (no
    // session_reset, which would replay the transcript on top of that history); the next one retries.
    if (stored) {
      p.keepHome = true;
      parked.set(sessionKey, stored);
    }
    dispose(p, 'failed to start');
    throw startError(e);
  }
  return p;
}

export function disposeAll(reason: string): void {
  for (const p of [...all])
    dispose(p, reason);
}

export function startPool(): void {
  fs.mkdirSync(HERMES_HOMES, {recursive: true, mode: 0o700});
  // Homes of a previous run (container restart, crash) belong to processes that no longer exist.
  for (const entry of fs.readdirSync(HERMES_HOMES))
    fs.rmSync(path.join(HERMES_HOMES, entry), {recursive: true, force: true});
  checkProvider();
  topUpSpares();
  setInterval(() => reapIdle(), 60_000).unref();
  startOrphanReaper(ownedPids);
}
