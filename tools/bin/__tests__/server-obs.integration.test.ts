/**
 * Integration tests for the `grok s` observability examples (alerts, errors, logger, capture,
 * timeline): every example line of GROK_S.md "Alerts, errors, logging and capture" runs through the
 * real CLI against a live server and its comment is asserted.
 *
 * Environment variables:
 *   GROK_IT_HOST     Server alias (required; the suite is skipped without it)
 *   GROK_IT_HOST2    A second deployment's alias, for the multi-host examples
 *   GROK_IT_SHARED   1 on a shared stand (implied for aliases other than wt-*, obs-*, local*, localhost
 *                    URLs): no users, groups, reports or seeded errors, no alert state changes, no
 *                    long error windows; only admin-scoped, time-boxed changes, undone in afterAll
 *
 * Run: GROK_IT_HOST=wt-GROK-20884-obs GROK_IT_HOST2=obs-t npx vitest run --project integration server-obs
 *
 * Everything it creates is named ObsIt<tag> and is removed (or blocked, for users) in afterAll.
 */

import {afterAll, beforeAll, describe, expect, it} from 'vitest';
import {spawn} from 'child_process';
import * as fs from 'fs';
import * as http from 'http';
import * as os from 'os';
import * as path from 'path';
import {randomBytes, randomUUID} from 'crypto';
import {NodeApiClient} from '../utils/node-dapi';
import {getServerCredentials} from '../utils/keypair';

const HOST = process.env['GROK_IT_HOST'] ?? '';
const HOST2 = process.env['GROK_IT_HOST2'] ?? '';
const SHARED = process.env['GROK_IT_SHARED'] === '1' ||
  (HOST !== '' && !/^(wt-|obs-|local|https?:\/\/(localhost|127\.))/.test(HOST));
const GROK = path.resolve(__dirname, '..', 'grok.js');
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'obsit-'));
const TAG = `ObsIt${Date.now().toString(36)}`;
const W7 = SHARED ? '1d' : '7d';
const LONG = 180_000;
const OBS = ['alerts', 'problems', 'errors', 'logger', 'capture', 'timeline'];

interface Run {code: number | null; out: string; err: string; bytes: Buffer}

/** `grok s observe <args>` for the observability commands, else `grok s <args>`, against [host] (HOST unless the args name their own `--host`, or false for none). */
function grok(args: string[], opts: {host?: string | false; env?: Record<string, string>; timeoutMs?: number} = {}): Promise<Run> {
  const host = opts.host === undefined ? HOST : opts.host;
  const full = ['s', ...(OBS.includes(args[0]) ? ['observe'] : []), ...args, ...(host && !args.includes('--host') ? ['--host', host] : [])];
  return new Promise((resolve) => {
    const child = spawn(process.execPath, [GROK, ...full], {env: {...process.env, ...opts.env}, cwd: TMP});
    const out: Buffer[] = [];
    const err: Buffer[] = [];
    child.stdout.on('data', (d) => out.push(d));
    child.stderr.on('data', (d) => err.push(d));
    const timer = setTimeout(() => child.kill(), opts.timeoutMs ?? 150_000);
    child.on('close', (code) => {
      clearTimeout(timer);
      const bytes = Buffer.concat(out);
      resolve({code, out: bytes.toString('utf8'), err: Buffer.concat(err).toString('utf8'), bytes});
    });
  });
}

async function ok(args: string[], opts?: Parameters<typeof grok>[1]): Promise<Run> {
  const r = await grok(args, opts);
  if (r.code !== 0)
    throw new Error(`grok s ${args.join(' ')} exited ${r.code}\n${r.out}\n${r.err}`);
  return r;
}

async function json(args: string[], opts?: Parameters<typeof grok>[1]): Promise<any> {
  const r = await ok([...args, '--output', 'json'], opts);
  try {
    return JSON.parse(r.out);
  }
  catch {
    throw new Error(`grok s ${args.join(' ')} --output json printed no JSON:\n${r.out}`);
  }
}

/** A table the CLI printed: columns from the dashed separator line. */
function table(out: string): {columns: string[]; rows: Record<string, string>[]} {
  const lines = out.split(/\r?\n/);
  const sep = lines.findIndex((l) => /^-+( +-+)*\s*$/.test(l) && l.trim().length > 0);
  if (sep < 1) return {columns: [], rows: []};
  const spans: [number, number][] = [];
  const re = /-+/g;
  let m: RegExpExecArray | null;
  while ((m = re.exec(lines[sep])) !== null) spans.push([m.index, m.index + m[0].length]);
  const cell = (l: string, i: number) => l.slice(spans[i][0], i + 1 < spans.length ? spans[i + 1][0] : undefined).trim();
  const columns = spans.map((_, i) => cell(lines[sep - 1], i));
  const rows = lines.slice(sep + 1).filter((l) => l.trim()).map((l) => Object.fromEntries(columns.map((c, i) => [c, cell(l, i)])));
  return {columns, rows};
}

const CROCKFORD = '0123456789ABCDEFGHJKMNPQRSTVWXYZ';

function ulid(): string {
  let t = Date.now();
  let time = '';
  for (let i = 0; i < 10; i++) { time = CROCKFORD[t % 32] + time; t = Math.floor(t / 32); }
  return time + [...randomBytes(16)].map((b) => CROCKFORD[b % 32]).join('');
}

const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

async function until<T>(what: string, fn: () => Promise<T | undefined | null | false>, ms: number = 90_000): Promise<T> {
  const end = Date.now() + ms;
  for (;;) {
    const v = await fn();
    if (v) return v as T;
    if (Date.now() > end) throw new Error(`Timed out waiting for ${what}`);
    await sleep(3000);
  }
}

/** Client log records over `/log/stream`, as the browser posts them, with [client]'s session. */
async function postClientErrors(client: NodeApiClient, messages: {message: string; stack: string; params?: any; n?: number}[]): Promise<void> {
  const ws = new WebSocket(`${client.baseUrl.replace(/^http/, 'ws')}/log/stream`,
    {headers: {cookie: `auth=${encodeURIComponent(client.token.startsWith('Bearer ') ? client.token : `Bearer ${client.token}`)}`}} as any);
  await new Promise<void>((resolve, reject) => {
    ws.onmessage = () => resolve();
    ws.onerror = () => reject(new Error('log/stream refused the connection'));
  });
  for (const m of messages)
    for (let i = 0; i < (m.n ?? 1); i++)
      ws.send(JSON.stringify({level: 'error', message: m.message, params: {...m.params}, stackTrace: m.stack,
        time: Date.now(), severity: 'normal'}));
  await sleep(1500);
  ws.close();
}

async function request(client: NodeApiClient, method: string, p: string, requestId?: string, body?: any): Promise<Response> {
  return fetch(`${client.baseUrl}${p}`, {method, body: body === undefined ? undefined : JSON.stringify(body),
    headers: {'Authorization': client.token, 'Content-Type': 'application/json', ...(requestId ? {'x-request-id': requestId} : {})}});
}

/** A CLI home whose only server is [url] with [devKey], for a run as another user. */
function cliHomeFor(url: string, devKey: string): Record<string, string> {
  const home = fs.mkdtempSync(path.join(TMP, 'home-'));
  fs.mkdirSync(path.join(home, '.grok'));
  fs.writeFileSync(path.join(home, '.grok', 'config.yaml'), `default: obsit\nservers:\n  obsit:\n    url: ${url}\n    key: ${devKey}\n`);
  return {USERPROFILE: home, HOME: home, GROK_PRIVATE_KEY: ''};
}

const todayUtc = () => new Date().toISOString().slice(0, 10);
const daysAgoUtc = (n: number) => new Date(Date.now() - n * 86400000).toISOString().slice(0, 10);

// ─── State ──────────────────────────────────────────────────────────────────

let admin: NodeApiClient;
let keepAlive: any;
let baseUrl = '';
const seed: {
  alice?: {id: string; login: string; client: NodeApiClient; home: Record<string, string>; session?: string};
  bob?: {id: string; login: string; client: NodeApiClient; home: Record<string, string>};
  group?: {id: string; name: string};
  sigShared?: string; sigSolo?: string; sigIncident?: string; sigT?: string;
  pkg?: string; conn?: {id: string; name: string}; query?: string;
  incident?: any; report?: {id: string; number: number};
  action?: string; request?: string;
  overrides: string[]; rules: string[]; jobs: string[]; files: string[];
} = {overrides: [], rules: [], jobs: [], files: []};

async function createUser(name: string): Promise<{id: string; login: string; client: NodeApiClient; home: Record<string, string>}> {
  const login = `${TAG.toLowerCase()}.${name}`;
  const file = path.join(TMP, `${name}.json`);
  fs.writeFileSync(file, JSON.stringify({'#type': 'User', login, firstName: name, lastName: TAG, status: 'active'}));
  const user = await json(['users', 'save', '--json', file]);
  const key = await admin.get(`/users/${user.id}/dev_key?getNew=true`);
  const devKey = typeof key === 'string' ? key : key?.key ?? key?.devKey;
  return {id: user.id, login, client: await NodeApiClient.login(baseUrl, devKey), home: cliHomeFor(baseUrl, devKey)};
}

async function errorsJson(args: string[]): Promise<any[]> {
  return await json(['errors', ...args]);
}

async function signatureOf(message: string, since: string = '1h'): Promise<string | undefined> {
  const rows = await errorsJson(['list', '--since', since, '--limit', '500']);
  return rows.find((r) => String(r.error ?? '').includes(message))?.signature;
}

describe.skipIf(!HOST)('grok s observability examples', () => {
  beforeAll(async () => {
    const cred = getServerCredentials(HOST);
    baseUrl = cred.url;
    admin = await NodeApiClient.login(cred.url, cred.key ?? '', cred.privateKey);
    const act = ulid();
    seed.action = act;
    seed.request = `${act}.1`;
    await request(admin, 'GET', '/users/current', `${act}.1`);
    await request(admin, 'GET', '/info/server', `${act}.2`);
    if (SHARED) return;

    seed.pkg = `${TAG}Pkg`;
    seed.alice = await createUser('alice');
    seed.bob = await createUser('bob');
    const gfile = path.join(TMP, 'group.json');
    const gname = `${TAG}Chemists`;
    fs.writeFileSync(gfile, JSON.stringify({'#type': 'UserGroup', name: gname, friendlyName: gname}));
    seed.group = {id: (await json(['groups', 'save', '--json', gfile])).id, name: gname};
    await ok(['groups', 'add-members', gname, seed.alice.login, seed.bob.login, '--user']);
    seed.alice.session = JSON.parse(Buffer.from(seed.alice.client.token.replace(/^Bearer /, '').split('.')[1], 'base64url').toString()).id;

    const pkgParams = {packageName: seed.pkg, packageVersion: '1.0.0'};
    const loopParams = {packageName: `${TAG}LoopPkg`, packageVersion: '1.0.0'};
    await postClientErrors(seed.alice.client, [
      {message: `${TAG}Shared: grid lost its filters`, stack: `Error: ${TAG}Shared\n    at ${TAG}filters (${TAG}.js:10:1)`, params: pkgParams},
      {message: `${TAG}Solo: submit failed`, stack: `Error: ${TAG}Solo\n    at ${TAG}submit (${TAG}.js:20:1)`, params: {...pkgParams, reqId: `${act}.3`}},
    ]);
    await request(seed.alice.client, 'GET', '/users/current', `${act}.3`);
    await postClientErrors(seed.bob.client, [
      {message: `${TAG}Shared: grid lost its filters`, stack: `Error: ${TAG}Shared\n    at ${TAG}filters (${TAG}.js:10:1)`, params: pkgParams},
    ]);
    await postClientErrors(admin, [
      {message: `${TAG}Loop: render loop`, stack: `Error: ${TAG}Loop\n    at ${TAG}render (${TAG}.js:30:1)`, params: loopParams, n: 22},
    ]);

    const conn = `${TAG}DeadPg`;
    const cfile = path.join(TMP, 'conn.json');
    fs.writeFileSync(cfile, JSON.stringify({'#type': 'DataConnection', name: conn, friendlyName: conn, dataSource: 'Postgres',
      parameters: {server: 'localhost', port: 1, db: 'none'}, credentials: {parameters: {login: 'x', password: 'y'}}}));
    seed.conn = {id: (await json(['connections', 'save', '--json', cfile, '--save-credentials'])).id, name: conn};
    const query = await admin.post('/connectors/queries', {'#type': 'DataQuery', name: `${TAG}DeadQuery`, friendlyName: `${TAG}DeadQuery`,
      query: 'select 1', connection: {id: seed.conn.id}, params: []});
    seed.query = query.id;
    const run = request(admin, 'POST', `/public/v1/functions/${encodeURIComponent(`${String(query.namespace ?? '').replace(':', '.')}${TAG}DeadQuery`)}/call`, ulid() + '.1', {});

    const reportBody = {'#type': 'UserReport', type: 'report', description: `${TAG} filters lost after submit`,
      errorMessage: `${TAG}Report: filters lost`, errorStackTrace: `Error: ${TAG}Report\n    at ${TAG}report (${TAG}.js:40:1)`,
      data: {'#type': 'UserReportData', id: randomUUID(), details: {requestIds: [`${act}.3`]}}};
    const res = await request(seed.alice.client, 'POST', '/reports?sendEmail=false', undefined, reportBody);
    const report = await res.json() as any;
    seed.report = {id: report.id, number: report.number};

    seed.sigShared = await until('the shared signature', () => signatureOf(`${TAG}Shared`), 30_000);
    seed.sigSolo = await until('the solo signature', () => signatureOf(`${TAG}Solo`), 30_000);
    seed.sigIncident = await until('the loop signature', () => signatureOf(`${TAG}Loop`), 30_000);
    seed.incident = await until('the error-incident alert', async () => {
      const alerts: any[] = await json(['alerts', 'list', '--kind', 'error-incident']);
      return alerts.find((a) => a.key === seed.sigIncident!.slice(0, 12));
    }, 150_000);
    const loop = {message: `${TAG}Loop: render loop`, stack: `Error: ${TAG}Loop\n    at ${TAG}render (${TAG}.js:30:1)`, params: loopParams};
    keepAlive = setInterval(() => postClientErrors(admin, [loop]).catch(() => {}), 40_000);
    await run;
  }, 420_000);

  afterAll(async () => {
    clearInterval(keepAlive);
    for (const id of seed.overrides)
      await grok(['logger', 'revert', '--override', id, '--reason', 'obsit cleanup']);
    for (const id of seed.rules)
      await grok(['capture', 'stop', id, '--reason', 'obsit cleanup']);
    for (const id of seed.jobs)
      await admin?.request('DELETE', `/connectors/jobs/${id}`).catch(() => {});
    for (const f of seed.files)
      await grok(['files', 'delete', f]);
    if (seed.incident) {
      await grok(['alerts', 'unmute', seed.incident.id]);
      // The loop's errors stay in the detection window for a while: a reopened incident stays muted.
      await grok(['alerts', 'mute', seed.incident.id, '--for', '30m', '--reason', 'obsit cleanup']);
      await grok(['alerts', 'resolve', seed.incident.id, '--reason', 'obsit cleanup']);
    }
    if (seed.report) {
      const open: any[] = (await admin.get(`/alerts?kind=report&key=${seed.report.id}`)) ?? [];
      for (const a of open.filter((x) => x.status !== 'resolved'))
        await grok(['alerts', 'resolve', a.id, '--reason', 'obsit cleanup']);
    }
    if (seed.query) await admin.request('DELETE', `/connectors/queries/${seed.query}`).catch(() => {});
    if (seed.conn) await grok(['connections', 'delete', seed.conn.id]);
    if (seed.group) await grok(['groups', 'delete', seed.group.id]);
    for (const u of [seed.alice, seed.bob])
      if (u) await grok(['users', 'block', u.login]);
    fs.rmSync(TMP, {recursive: true, force: true});
  }, 300_000);

  // ─── Alerts ───────────────────────────────────────────────────────────────

  describe('alerts', () => {
    it('alerts list  # open, acknowledged, muted', async () => {
      const r = await ok(['alerts', 'list']);
      const t = table(r.out);
      if (t.rows.length) expect(t.columns).toEqual(['KIND', 'KEY', 'SEV', 'AUDIENCE', 'STATUS', 'OPENED', 'BY', 'SUMMARY']);
      for (const row of t.rows) expect(['open', 'acknowledged', 'muted']).toContain(row.STATUS);
      const all: any[] = await json(['alerts', 'list']);
      for (const a of all) expect(['open', 'acknowledged', 'muted']).toContain(a.status);
      if (seed.incident) expect(all.map((a) => a.id)).toContain(seed.incident.id);
    }, LONG);

    it('alerts list --status all --since 1h', async () => {
      const r = await ok(['alerts', 'list', '--status', 'all', '--since', '1h']);
      expect(r.out).toMatch(/KIND|\(no results\)/);
      const rows: any[] = await json(['alerts', 'list', '--status', 'all', '--since', '1h']);
      const hourAgo = Date.now() - 3600_000 - 60_000;
      for (const a of rows)
        expect(Math.max(Date.parse(a.openedAt), Date.parse(a.lastSeen ?? a.openedAt), Date.parse(a.resolvedAt ?? 0))).toBeGreaterThanOrEqual(hourAgo);
    }, LONG);

    it('alerts detection  # servers and their liveness', async () => {
      const t = table((await ok(['alerts', 'detection'])).out);
      expect(t.columns).toEqual(['SERVER', 'HOST NAME', 'VERSION', 'LAST SEEN', 'LIVE']);
      expect(t.rows.length).toBeGreaterThan(0);
      expect(t.rows.some((r) => r.LIVE === 'yes')).toBe(true);
      expect(Object.keys(await json(['alerts', 'detection']))).toEqual(['servers']);
    }, LONG);

    it('problems history <id>  # the records of a problem', async () => {
      const problems: any[] = await json(['problems', 'list', '--status', 'all', '--limit', '1']);
      if (!problems.length) return;
      const t = table((await ok(['problems', 'history', problems[0].id])).out);
      if (t.columns.length) expect(t.columns).toEqual(['TIME', 'RECORD', 'ALERT', 'SUMMARY']);
    }, LONG);

    it.skipIf(!HOST2)('GROK_S.md: alerts detection --all --host a --host b  # every server row, stopped ones included', async () => {
      const t = table((await ok(['alerts', 'detection', '--all', '--host', HOST, '--host', HOST2])).out);
      expect(t.columns[0]).toBe('HOST');
      expect(t.columns.slice(1)).toEqual(['SERVER', 'HOST NAME', 'VERSION', 'LAST SEEN', 'LIVE']);
      const hosts = new Set(t.rows.map((r) => r.HOST));
      expect(hosts.size).toBeGreaterThanOrEqual(1);
      const each: any = await json(['alerts', 'detection', '--all']);
      expect(t.rows.filter((r) => r.HOST.split(', ').includes(HOST)).length).toBe(each.servers.length);
    }, LONG);

    it.skipIf(SHARED)('alerts get <kind:key | id | 6+ char prefix>', async () => {
      const a = seed.incident;
      const ids: string[] = ((await json(['alerts', 'list', '--status', 'all', '--limit', '1000'])) as any[]).map((x) => x.id);
      let n = 6;
      while (ids.filter((x) => x.startsWith(String(a.id).slice(0, n))).length > 1) n++;
      if (n > 6) {
        const r = await grok(['alerts', 'get', String(a.id).slice(0, 6)]);
        expect(r.code).toBe(1);
        expect(r.out + r.err).toContain('matches several alerts');
      }
      for (const id of [`${a.kind}:${a.key}`, `${a.kind}:${String(a.key).slice(0, 6)}`, a.id, String(a.id).slice(0, n)]) {
        const r = await ok(['alerts', 'get', id]);
        expect(r.out).toContain(`${a.kind}:${a.key}`);
        expect(r.out).toMatch(/^id\s+/m);
        expect((await json(['alerts', 'get', id])).id).toBe(a.id);
      }
      const health = ((await json(['alerts', 'list'])) as any[]).find((x) => x.kind === 'health');
      if (health) expect((await ok(['alerts', 'get', `health:${health.key}`])).out).toContain(`health:${health.key}`);
    }, LONG);

    it.skipIf(SHARED)('alerts ack <kind:key> --reason', async () => {
      const a = seed.incident;
      const r = await ok(['alerts', 'ack', `${a.kind}:${a.key}`, '--reason', 'restarting grok_connect']);
      expect(r.out.trim()).toBe(`acknowledged ${a.kind}:${a.key} — restarting grok_connect`);
      expect((await json(['alerts', 'get', a.id])).status).toBe('acknowledged');
    }, LONG);

    it.skipIf(SHARED)('alerts mute --until 14:00 (bare time is UTC)', async () => {
      const a = seed.incident;
      const later = new Date(Date.now() + 2 * 3600_000);
      if (later.toISOString().slice(0, 10) !== todayUtc()) return;
      const hhmm = later.toISOString().slice(11, 16);
      const r = await ok(['alerts', 'mute', `${a.kind}:${String(a.key).slice(0, 6)}`, '--until', hhmm, '--reason', 'hotfix deploying']);
      expect(r.out.trim()).toBe(`muted ${a.kind}:${a.key} until ${hhmm}Z — hotfix deploying`);
      expect((await json(['alerts', 'get', a.id])).status).toBe('muted');
      const top: any[] = await errorsJson(['top', '--since', '1h', '--signature', seed.sigIncident!]);
      expect(top[0]?.state).toBe('muted');
      expect(top[0]?.stateUntil).toBe(`${todayUtc()}T${hhmm}:00.000Z`);
    }, LONG);

    it.skipIf(SHARED)('GROK_S.md: alerts mute --until 2026-10-04T06:00 (ISO without offset is UTC)', async () => {
      const a = seed.incident;
      const at = new Date(Date.now() + 3 * 86400_000).toISOString().slice(0, 16);
      const r = await ok(['alerts', 'mute', a.id, '--until', at, '--reason', 'monthly ELN maintenance'], {env: {TZ: 'America/Los_Angeles'}});
      expect(r.out.trim()).toBe(`muted ${a.kind}:${a.key} until ${at.slice(5, 10)} ${at.slice(11)}Z — monthly ELN maintenance`);
      const top: any[] = await errorsJson(['top', '--since', '1h', '--signature', seed.sigIncident!]);
      expect(top[0]?.stateUntil).toBe(`${at}:00.000Z`);
    }, LONG);

    it.skipIf(SHARED)('alerts unmute <kind:key>', async () => {
      const a = seed.incident;
      const r = await ok(['alerts', 'unmute', `${a.kind}:${String(a.key).slice(0, 6)}`]);
      expect(r.out.trim()).toBe(`unmuted ${a.kind}:${a.key}`);
      expect((await json(['alerts', 'get', a.id])).status).not.toBe('muted');
    }, LONG);

    it.skipIf(SHARED)('alerts mute --until-version 1.14.3 --reason', async () => {
      const a = seed.incident;
      const r = await ok(['alerts', 'mute', `${a.kind}:${String(a.key).slice(0, 6)}`, '--until-version', '1.14.3', '--reason', 'fixed in Chem 1.14.3']);
      expect(r.out.trim()).toBe(`muted ${a.kind}:${a.key} until 1.14.3 — fixed in Chem 1.14.3`);
      const top: any[] = await errorsJson(['top', '--since', '1h', '--signature', seed.sigIncident!]);
      expect(top[0]?.state).toBe('muted');
      expect(top[0]?.stateVersion).toBe('1.14.3');
      expect(table((await ok(['errors', 'top', '--since', '1h', '--signature', seed.sigIncident!])).out).rows[0].STATE).toBe('muted → 1.14.3');
      await ok(['alerts', 'unmute', a.id]);
    }, LONG);

    it.skipIf(SHARED)('alerts resolve report:<number> --reason', async () => {
      const r = await ok(['alerts', 'resolve', `report:${seed.report!.number}`, '--reason', 'duplicate of 4819']);
      expect(r.out.trim()).toBe(`resolved report:${seed.report!.id} — duplicate of 4819`);
      const rows: any[] = await admin.get(`/alerts?kind=report&key=${seed.report!.id}&status=all`);
      expect(rows[0]?.status).toBe('resolved');
      expect(rows[0]?.resolveReason).toBe('duplicate of 4819');
    }, LONG);

    it.skipIf(SHARED)('alerts resolve report:<missing number>  # not found, never a key prefix', async () => {
      const r = await grok(['alerts', 'resolve', 'report:98765432', '--reason', 'no such report']);
      expect(r.code).toBe(1);
      expect(r.out + r.err).toContain('No open alert report:98765432');
    }, LONG);

    it.skipIf(!HOST2)('alerts list --host a --host b  # several stands, HOST column', async () => {
      const r = await ok(['alerts', 'list', '--status', 'all', '--since', '7d', '--host', HOST, '--host', HOST2]);
      const t = table(r.out);
      expect(t.columns[0]).toBe('HOST');
      expect(new Set(t.rows.map((x) => x.HOST))).toEqual(new Set([HOST, HOST2].filter((h) => t.rows.some((x) => x.HOST === h))));
      const rows: any[] = await json(['alerts', 'list', '--status', 'all', '--since', '7d', '--host', HOST, '--host', HOST2]);
      for (const a of rows) expect([HOST, HOST2]).toContain(a.host);
    }, LONG);
  });

  // ─── Errors ───────────────────────────────────────────────────────────────

  describe('errors', () => {
    it('errors list --since 2h  # occurrences', async () => {
      const t = table((await ok(['errors', 'list', '--since', '2h'])).out);
      if (!SHARED || t.rows.length)
        expect(t.columns).toEqual(['TIME', 'USER', 'SOURCE', 'SIG', 'ERROR', 'PACKAGE', 'VERSION', 'ROUTE', 'SERVER', 'REQ']);
      if (!SHARED) expect(t.rows.some((r) => r.ERROR.includes(TAG))).toBe(true);
      const rows: any[] = await errorsJson(['list', '--since', '2h']);
      const cutoff = Date.now() - 2 * 3600_000 - 60_000;
      for (const r of rows) expect(Date.parse(r.time)).toBeGreaterThanOrEqual(cutoff);
    }, LONG);

    it.skipIf(SHARED)('errors list --user alice --since 1d', async () => {
      const t = table((await ok(['errors', 'list', '--user', seed.alice!.login, '--since', '1d'])).out);
      expect(t.rows.length).toBeGreaterThanOrEqual(2);
      for (const r of t.rows) expect(r.USER).toBe(seed.alice!.login);
      expect(t.rows.some((r) => r.ERROR.includes(`${TAG}Solo`))).toBe(true);
    }, LONG);

    it.skipIf(SHARED)('errors top --since 7d --by signature,package --min-users 2 --package <this run package>', async () => {
      const t = table((await ok(['errors', 'top', '--since', '7d', '--by', 'signature,package', '--min-users', '2', '--package', seed.pkg!])).out);
      expect(t.columns.slice(0, 7)).toEqual(['SIG', 'PACKAGE', 'ERROR', 'USERS', 'COUNT', 'FIRST SEEN', 'LAST']);
      const mine = t.rows.find((r) => r.ERROR.includes(`${TAG}Shared`));
      expect(mine).toBeTruthy();
      expect(mine!.USERS).toBe('2');
      expect(mine!.PACKAGE).toBe(seed.pkg);
      expect(mine!.SIG).toBe(seed.sigShared!.slice(0, 6));
      for (const r of t.rows) expect(Number(r.USERS)).toBeGreaterThanOrEqual(2);
      expect(t.rows.some((r) => r.ERROR.includes(`${TAG}Solo`))).toBe(false);
    }, LONG);

    it.skipIf(SHARED)('errors top --since 30d --group Chemists --by package', async () => {
      const t = table((await ok(['errors', 'top', '--since', '30d', '--group', seed.group!.name, '--by', 'package'])).out);
      expect(t.columns[0]).toBe('PACKAGE');
      expect(t.rows.map((r) => r.PACKAGE)).toContain(seed.pkg);
      const rows: any[] = await errorsJson(['top', '--since', '30d', '--group', seed.group!.name, '--by', 'package']);
      const mine = rows.find((r) => r.package === seed.pkg);
      expect(mine.users).toBe(2);
      expect(mine.count).toBe(3);
    }, LONG);

    it.skipIf(SHARED)('errors top --since 1h --route "POST /api/public/v1/functions/{name}/call" --by connection', async () => {
      const route = 'POST /api/public/v1/functions/{name}/call';
      const rows: any[] = await until('the failed query run', async () => {
        const r: any[] = await errorsJson(['top', '--since', '1h', '--route', route, '--by', 'connection']);
        return r.some((x) => String(x.connection).endsWith(`${TAG}DeadPg`)) ? r : undefined;
      }, 60_000);
      expect(rows.some((r) => String(r.connection).endsWith(`${TAG}DeadPg`))).toBe(true);
      const t = table((await ok(['errors', 'top', '--since', '1h', '--route', route, '--by', 'connection'])).out);
      expect(t.columns[0]).toBe('CONNECTION');
      expect(t.rows.some((r) => r.CONNECTION.endsWith(`${TAG}DeadPg`))).toBe(true);
    }, LONG);

    it.skipIf(SHARED)('errors show <6-char signature> --since 30d  # one signature in detail', async () => {
      const r = await ok(['errors', 'show', seed.sigShared!.slice(0, 6), '--since', '30d']);
      expect(r.out).toMatch(new RegExp(`^signature\\s+${seed.sigShared!.slice(0, 6)}\\s+${TAG}Shared`, 'm'));
      expect(r.out).toMatch(new RegExp(`^package\\s+${seed.pkg}\\s+first seen in 1\\.0\\.0`, 'm'));
      expect(r.out).toMatch(/^occurrences\s+2\s+users 2\s+groups .*/m);
      expect(r.out).toContain(seed.group!.name);
      const doc = await json(['errors', 'show', seed.sigShared!.slice(0, 6), '--since', '30d']);
      expect(doc.signature).toBe(seed.sigShared);
      expect(doc.users).toBe(2);
    }, LONG);

    it.skipIf(SHARED)('errors diff --before <7 days>..<yesterday> --after <today>..<today>', async () => {
      const args = ['errors', 'diff', '--before', `${daysAgoUtc(7)}..${daysAgoUtc(1)}`, '--after', `${todayUtc()}..${todayUtc()}`];
      const r = await ok(args);
      for (const label of ['NEW', 'GONE', 'RISEN', 'REGRESSED', 'incidents', 'reports'])
        expect(r.out).toMatch(new RegExp(`^${label}\\s`, 'm'));
      const doc = await json(args);
      expect(doc.categories.new.count).toBeGreaterThanOrEqual(3);
      const news = (doc.rows ?? []).filter((x: any) => x.category === 'new').map((x: any) => x.signature);
      expect(news).toEqual(expect.arrayContaining([seed.sigShared, seed.sigSolo]));
    }, LONG);

    it.skipIf(!HOST2 || SHARED)('errors diff --since 7d --host a --host b  # compare two stands', async () => {
      const cred2 = getServerCredentials(HOST2);
      const other = await NodeApiClient.login(cred2.url, cred2.key ?? '', cred2.privateKey);
      await postClientErrors(other, [{message: `${TAG}OnlyT: only on the second stand`, stack: `Error: ${TAG}OnlyT\n    at ${TAG}t (${TAG}t.js:1:1)`}]);
      const doc = await until('the second stand\'s error', async () => {
        const d = await json(['errors', 'diff', '--since', '7d', '--host', HOST, '--host', HOST2]);
        return d.hosts[1].only.some((x: any) => String(x.topError).includes(`${TAG}OnlyT`)) ? d : undefined;
      }, 60_000);
      expect(doc.hosts.map((h: any) => h.host)).toEqual([HOST, HOST2]);
      expect(doc.hosts[0].only.some((x: any) => String(x.topError).includes(`${TAG}OnlyT`))).toBe(false);
      expect(doc.hosts[0].only.some((x: any) => x.signature === seed.sigShared)).toBe(true);
      const r = await ok(['errors', 'diff', '--since', '7d', '--host', HOST, '--host', HOST2]);
      expect(r.out).toMatch(new RegExp(`^ONLY ON ${HOST} \\(.+\\)\\s+\\d+`, 'm'));
      expect(r.out).toMatch(new RegExp(`^ONLY ON ${HOST2} \\(.+\\)\\s+[1-9]\\d*`, 'm'));
      expect(r.out).toMatch(/^ON BOTH\s+\d+/m);
    }, LONG);

    it(`errors export --since ${W7} --by signature,package --format csv -O errors.csv`, async () => {
      const file = path.join(TMP, 'errors.csv');
      const r = await ok(['errors', 'export', '--since', W7, '--by', 'signature,package', '--format', 'csv', '-O', file]);
      expect(r.out).toMatch(/^Wrote \d+ bytes to /m);
      const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(Boolean);
      const header = lines[0].split(',');
      expect(header).toEqual(expect.arrayContaining(['signature', 'package']));
      if (!SHARED) {
        expect(lines.length).toBeGreaterThan(1);
        expect(lines.some((l) => l.includes(seed.sigShared!))).toBe(true);
      }
    }, LONG);

    it(`errors top --since ${W7} --by package,group --format parquet > errors.parquet`, async () => {
      const r = await ok(['errors', 'top', '--since', W7, '--by', 'package,group', '--format', 'parquet']);
      const bytes = r.bytes;
      expect(bytes.subarray(0, 4).toString()).toBe('PAR1');
      expect(bytes.subarray(bytes.length - 4).toString()).toBe('PAR1');
      const rows: any[] = JSON.parse((await ok(['errors', 'top', '--since', W7, '--by', 'package,group', '--format', 'json'])).out);
      const parquet = require('parquet-wasm');
      const arrow = require('apache-arrow');
      const t = arrow.tableFromIPC(parquet.readParquet(new Uint8Array(bytes)).intoIPCStream());
      expect(t.numRows).toBe(rows.length);
      if (rows.length) expect(t.schema.fields.map((f: any) => f.name)).toEqual(expect.arrayContaining(['package', 'group']));
    }, LONG);

    it('errors save "Errors by team, weekly" --since 7d --by group,package --schedule "MON 07:00" --to "System:AppData/Ops/errors/"', async () => {
      const name = `${TAG} errors by team, weekly`;
      const r = await ok(['errors', 'save', name, '--since', '7d', '--by', 'group,package', '--schedule', 'MON 07:00', '--to', 'System:AppData/Ops/errors/']);
      const slug = name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
      const target = `System:AppData/Ops/errors/${slug}-{date}.csv`;
      expect(r.out.trim()).toBe(`saved job "${name}" 0 7 * * 1 → ${target}`);
      const jobs: any[] = await admin.get(`/connectors/jobs?text=${encodeURIComponent(TAG)}`);
      const job = jobs.find((j) => j.friendlyName === name);
      expect(job).toBeTruthy();
      seed.jobs.push(job.id);
      const full = await admin.get(`/connectors/jobs/${job.id}?include=params`);
      expect(full.recurrence?.cronSchedule ?? full.recurrence?.cron).toBe('0 7 * * 1');
      const params = Object.fromEntries((full.params ?? []).map((p: any) => [p.name, p.defaultValue]));
      expect(params.path).toBe(target);
      expect(JSON.parse(params.spec)).toMatchObject({since: '7d', by: 'group,package'});
      const written = target.replace('{date}', todayUtc());
      seed.files.push(written);
      const pfile = path.join(TMP, 'export.json');
      fs.writeFileSync(pfile, JSON.stringify({spec: JSON.stringify({since: '1h', by: 'group,package'}), format: 'csv', path: target}));
      await ok(['functions', 'run', 'ErrorsExport', '--json', pfile]);
      const csv = (await ok(['files', 'get', written])).out;
      expect(csv.split(/\r?\n/)[0].split(',')).toEqual(expect.arrayContaining(['group', 'package']));
    }, LONG);
  });

  // ─── Logger ───────────────────────────────────────────────────────────────

  describe('logger', () => {
    let versionBefore: number;

    it('logger get server  # levels, flags, locks, overrides', async () => {
      const r = await ok(['logger', 'get', 'server']);
      expect(r.out).toMatch(/^print\s+.+post\s+/m);
      expect(r.out).toMatch(/^save\s+.+debug flags\s+/m);
      const policy = await json(['logger', 'get', 'server']);
      expect(policy).toHaveProperty('settings');
      expect(policy).toHaveProperty('locked');
      expect(policy).toHaveProperty('overrides');
      if ((policy.locked ?? []).length) expect(r.out).toMatch(/^locked\s+/m);
      const history: any[] = await json(['logger', 'history', '--limit', '1']);
      versionBefore = history[0]?.version ?? 0;
    }, LONG);

    it('logger get --scope user:admin  # effective settings, each with its source', async () => {
      const t = table((await ok(['logger', 'get', '--scope', 'user:admin'])).out);
      expect(t.columns).toEqual(['SETTING', 'VALUE', 'SOURCE']);
      expect(t.rows.map((r) => r.SETTING)).toEqual(expect.arrayContaining(['printLevels', 'saveLevels', 'debugFlags']));
      for (const r of t.rows) expect(r.SOURCE).not.toBe('');
    }, LONG);

    it('logger set server --debug-flags +query --scope user:admin --for 30m --reason "debug queries"', async () => {
      const before = await json(['logger', 'get', '--scope', 'user:admin']);
      const r = await ok(['logger', 'set', 'server', '--debug-flags', '+query', '--scope', 'user:admin', '--for', '30m', '--reason', 'debug queries']);
      expect(r.out).toMatch(/^\+ server\.debugFlags  .*query.*  scope user:admin  reverts \S+/m);
      const overrides: any[] = await json(['logger', 'overrides']);
      const mine = overrides.find((o) => o.scope === 'user' && o.reason === 'debug queries' && o.changes?.debugFlags?.includes('query'));
      expect(mine).toBeTruthy();
      seed.overrides.push(mine.id);
      const minutes = (Date.parse(mine.expiresAt) - Date.now()) / 60000;
      expect(minutes).toBeGreaterThan(28);
      expect(minutes).toBeLessThanOrEqual(30.5);
      const after = await json(['logger', 'get', '--scope', 'user:admin']);
      expect(after.settings.debugFlags).toEqual(expect.arrayContaining([...(before.settings.debugFlags ?? []), 'query']));
      expect(String(after.sources.debugFlags)).toMatch(/override/i);
    }, LONG);

    it.skipIf(SHARED)('logger set server --debug-flags +queries --scope package:Snowflake --for 30m', async () => {
      const r = await ok(['logger', 'set', 'server', '--debug-flags', '+queries', '--scope', `package:${TAG}Snowflake`, '--for', '30m']);
      expect(r.out).toMatch(new RegExp(`^\\+ server\\.debugFlags  .*query.*  scope package:${TAG}Snowflake  reverts `, 'm'));
      const overrides: any[] = await json(['logger', 'overrides']);
      const mine = overrides.find((o) => o.scope === 'package' && (o.scopeName ?? o.scopeId) === `${TAG}Snowflake`);
      expect(mine?.changes?.debugFlags).toContain('query');
      seed.overrides.push(mine.id);
      const eff = await json(['logger', 'get', '--scope', `package:${TAG}Snowflake`]);
      expect(eff.settings.debugFlags).toContain('query');
    }, LONG);

    it.skipIf(SHARED)('logger set --scope group:Chemists --print-levels error,warning --for 2h', async () => {
      const r = await ok(['logger', 'set', '--scope', `group:${seed.group!.name}`, '--print-levels', 'error,warning', '--for', '2h']);
      expect(r.out).toMatch(new RegExp(`^\\+ server\\.printLevels  error, ?warning  scope group:${seed.group!.name}  reverts `, 'm'));
      const overrides: any[] = await json(['logger', 'overrides']);
      const mine = overrides.find((o) => o.scope === 'group' && o.scopeId === seed.group!.id);
      expect(mine?.changes?.printLevels).toEqual(['error', 'warning']);
      seed.overrides.push(mine.id);
      const hours = (Date.parse(mine.expiresAt) - Date.now()) / 3600_000;
      expect(hours).toBeGreaterThan(1.9);
      const eff = await json(['logger', 'get', '--scope', `user:${seed.alice!.login}`]);
      expect(eff.settings.printLevels).toEqual(['error', 'warning']);
    }, LONG);

    it('logger overrides', async () => {
      const t = table((await ok(['logger', 'overrides'])).out);
      expect(t.columns).toEqual(['ID', 'SCOPE', 'CHANGES', 'REVERTS', 'BY', 'REASON']);
      expect(t.rows.some((r) => r.SCOPE === 'user:admin' && r.REASON === 'debug queries')).toBe(true);
    }, LONG);

    it('logger diff  # current settings vs deployment defaults', async () => {
      const r = await ok(['logger', 'diff']);
      expect(r.out).toMatch(/^\+ server\.debugFlags\s+.*query.*scope user:admin\s+reverts /m);
      const rows: any[] = await json(['logger', 'diff']);
      for (const row of rows) expect(['+', '-', '~']).toContain(row.change);
    }, LONG);

    it('logger diff --version <n>', async () => {
      if (!versionBefore) {
        const none = await grok(['logger', 'diff', '--version', '1']);
        expect(none.code).toBe(1);
        expect(none.out + none.err).toContain('Logger settings version 1 not found');
        return;
      }
      const r = await ok(['logger', 'diff', '--version', String(versionBefore)]);
      expect(r.out.length).toBeGreaterThan(0);
      const rows: any[] = await json(['logger', 'diff', '--version', String(versionBefore)]);
      expect(Array.isArray(rows)).toBe(true);
    }, LONG);

    it.skipIf(!HOST2)('GROK_S.md: logger diff --host a --host b', async () => {
      const r = await ok(['logger', 'diff', '--host', HOST, '--host', HOST2]);
      expect(r.out).toMatch(/^[+~-] server\.|\(no differences\)/m);
      const rows: any[] = await json(['logger', 'diff', '--host', HOST, '--host', HOST2]);
      for (const row of rows) expect(['+', '-', '~']).toContain(row.change);
    }, LONG);

    it.skipIf(SHARED)('GROK_S.md: logger set server --save-levels -debug --reason  # base change for All Users; logger set --set <path>=<json>', async () => {
      const base = await json(['logger', 'get', 'server']);
      const history: any[] = await json(['logger', 'history', '--limit', '1']);
      const before = history[0]?.version;
      const all = Object.entries(base.groups ?? {}).find(([, n]) => /^all users$/i.test(String(n)))![0];
      const hadDebug = (base.settings[`userGroupSettings.${all}.saveLevels`] ?? []).includes('debug');
      await sleep(6000);
      try {
        const r = await ok(['logger', 'set', 'server', '--save-levels', '-debug', '--reason', `${TAG} too much`]);
        if (hadDebug) {
          expect(r.out).toMatch(/^~ server\.saveLevels  /m);
          expect(r.out).toMatch(/^policy version \d+/m);
        }
        else
          expect(r.out.trim()).toBe('(no change: the base settings already have these values)');
        await sleep(6000);
        const after = await json(['logger', 'get', 'server']);
        expect(after.settings[`userGroupSettings.${all}.saveLevels`]).not.toContain('debug');
        const s = await ok(['logger', 'set', '--set', 'exportFlushSeconds=7', '--reason', `${TAG} faster sync`]);
        expect(s.out).toMatch(/^~ server\.exportFlushSeconds  7/m);
        expect((await json(['logger', 'get', 'server'])).settings.exportFlushSeconds).toBe(7);
        expect(((await json(['logger', 'history', '--limit', '2'])) as any[])[0].reason).toBe(`${TAG} faster sync`);
      }
      finally {
        if (((await json(['logger', 'history', '--limit', '1'])) as any[])[0]?.version !== before) {
          await sleep(6000);
          await ok(['logger', 'revert', String(before), '--reason', `${TAG} restore`]);
        }
      }
      const restored = await json(['logger', 'get', 'server']);
      expect(restored.settings).toEqual(base.settings);
    }, LONG);

    it('logger history --limit 20', async () => {
      const r = await ok(['logger', 'history', '--limit', '20']);
      if (!versionBefore) { expect(r.out.trim()).toBe('(no results)'); return; }
      const t = table(r.out);
      expect(t.columns).toEqual(['VERSION', 'CHANGED', 'BY', 'SOURCE', 'REASON', 'CHANGES']);
      expect(t.rows.length).toBeGreaterThan(0);
      expect(t.rows.length).toBeLessThanOrEqual(20);
      const versions = t.rows.map((r) => Number(r.VERSION));
      expect([...versions].sort((a, b) => b - a)).toEqual(versions);
    }, LONG);

    it('logger revert --override <id as listed>', async () => {
      const r0 = await ok(['logger', 'set', 'server', '--debug-flags', '+db', '--scope', 'user:admin', '--for', '10m', '--reason', `${TAG} revert-by-id`]);
      expect(r0.code).toBe(0);
      const listed = table((await ok(['logger', 'overrides'])).out).rows.find((x) => x.REASON === `${TAG} revert-by-id`);
      expect(listed).toBeTruthy();
      const full = ((await json(['logger', 'overrides'])) as any[]).find((o) => o.reason === `${TAG} revert-by-id`);
      seed.overrides.push(full.id);
      const r = await ok(['logger', 'revert', '--override', listed!.ID]);
      expect(r.out).toMatch(/^reverted /m);
      expect(((await json(['logger', 'overrides'])) as any[]).some((o) => o.id === full.id)).toBe(false);
    }, LONG);

    it('logger revert  # undo the latest change or override', async () => {
      await ok(['logger', 'set', 'server', '--debug-flags', '+socket', '--scope', 'user:admin', '--for', '10m', '--reason', `${TAG} latest`]);
      const created = ((await json(['logger', 'overrides'])) as any[]).find((o) => o.reason === `${TAG} latest`);
      seed.overrides.push(created.id);
      const overrides: any[] = await json(['logger', 'overrides']);
      const newest = overrides.reduce((a, b) => Date.parse(a.createdAt) >= Date.parse(b.createdAt) ? a : b);
      const lastBase = ((await json(['logger', 'history', '--limit', '1'])) as any[])[0];
      if (SHARED && (newest.id !== created.id || Date.parse(lastBase?.changedAt) > Date.parse(created.createdAt))) return;
      const r = await ok(['logger', 'revert']);
      expect(r.out).toMatch(/^reverted /m);
      expect(((await json(['logger', 'overrides'])) as any[]).some((o) => o.id === created.id)).toBe(false);
      expect(((await json(['logger', 'overrides'])) as any[]).some((o) => o.reason === 'debug queries')).toBe(true);
    }, LONG);
  });

  // ─── Capture and timelines ────────────────────────────────────────────────

  describe('capture and timeline', () => {
    let viewRule: any;
    let fnRule: any;

    it.skipIf(SHARED)('capture add --user alice --view "Queries" --capture clicks,requests,errors --for 30m --limit 2000 --reason', async () => {
      const r = await ok(['capture', 'add', '--user', seed.alice!.login, '--view', 'Queries', '--capture', 'clicks,requests,errors',
        '--for', '30m', '--limit', '2000', '--reason', 'GROK-21044: filters lost']);
      expect(r.out.trim()).toMatch(/^rule cap-\d+  active until \S+ \S+ · 1 user · 1 view · 0\/2000 events$/);
      const id = /cap-\d+/.exec(r.out)![0];
      seed.rules.push(id);
      viewRule = await json(['capture', 'show', id]);
      expect(viewRule.subject).toMatchObject({type: 'user'});
      expect(viewRule.scope).toEqual({type: 'view', value: 'Queries'});
      expect(viewRule.capture).toMatchObject({clicks: true, requests: true, errors: true, inputs: false, calls: false});
      expect(viewRule.maxEvents).toBe(2000);
      expect((Date.parse(viewRule.expiresAt) - Date.parse(viewRule.createdAt)) / 60000).toBeCloseTo(30, 0);
    }, LONG);

    it.skipIf(SHARED)('capture add --group Chemists --view "Hit Triage" --capture clicks,requests,errors --for 7d --anonymous --reason', async () => {
      const r = await ok(['capture', 'add', '--group', seed.group!.name, '--view', 'Hit Triage', '--capture', 'clicks,requests,errors',
        '--for', '7d', '--anonymous', '--reason', 'submit drop-off']);
      expect(r.out.trim()).toMatch(/^rule cap-\d+  active until .* · 1 group · 1 view · 0\/\d+ events$/);
      const id = /cap-\d+/.exec(r.out)![0];
      seed.rules.push(id);
      const rule = await json(['capture', 'show', id]);
      expect(rule.anonymous).toBe(true);
      expect(rule.subject.type).toBe('group');
      expect((Date.parse(rule.expiresAt) - Date.parse(rule.createdAt)) / 86400_000).toBeCloseTo(7, 1);
    }, LONG);

    it('capture add --user admin --function LoggingPolicy --capture calls,requests --for 10m --reason test', async () => {
      const r = await ok(['capture', 'add', '--user', 'admin', '--function', 'LoggingPolicy', '--capture', 'calls,requests', '--for', '10m', '--reason', 'test']);
      expect(r.out.trim()).toMatch(/^rule cap-\d+  active until .* · 1 user · 1 function · 0\/\d+ events$/);
      const id = /cap-\d+/.exec(r.out)![0];
      seed.rules.push(id);
      fnRule = await json(['capture', 'show', id]);
      await ok(['functions', 'run', 'LoggingPolicy()']);
      const events = await until('the captured LoggingPolicy call', async () => {
        const rows: any[] = await json(['timeline', '--rule', id]);
        return rows.some((e) => String(e.summary).includes('LoggingPolicy')) ? rows : undefined;
      }, 60_000);
      expect(events.length).toBeGreaterThan(0);
    }, LONG);

    it('GROK_S.md: capture add ... --capture clicks,inputs,requests,calls,errors,server:debug=queries,files --for 2d', async () => {
      const r = await ok(['capture', 'add', '--user', 'admin', '--view', 'Hit Triage', '--capture',
        'clicks,inputs,requests,calls,errors,server:debug=queries,files', '--for', '2d', '--limit', '2000', '--reason', `${TAG} campaign loses filters`]);
      const id = /cap-\d+/.exec(r.out)![0];
      seed.rules.push(id);
      const rule = await json(['capture', 'show', id]);
      expect(rule.capture).toMatchObject({clicks: true, inputs: true, requests: true, calls: true, errors: true, serverLevel: 'debug'});
      expect(rule.capture.debugFlags).toEqual(['query', 'storage']);
      expect((await ok(['capture', 'show', id])).out).toMatch(/^capture\s+clicks,inputs,requests,calls,errors,server:debug=query,storage$/m);
      await ok(['capture', 'stop', id, '--reason', 'obsit']);
    }, LONG);

    it('capture list --all --since 7d', async () => {
      const t = table((await ok(['capture', 'list', '--all', '--since', '7d'])).out);
      expect(t.columns).toEqual(['RULE', 'AUTHOR', 'SUBJECT', 'SCOPE', 'REASON', 'ACTIVE', 'EVENTS']);
      for (const id of seed.rules) expect(t.rows.map((r) => r.RULE)).toContain(id);
    }, LONG);

    it('capture show <cap-N>', async () => {
      const rule = viewRule ?? fnRule;
      const id = `cap-${rule.number}`;
      if (viewRule && seed.alice) {
        const res = await request(seed.alice.client, 'POST', `/logging/capture/${id}/activate`, undefined, {detail: 'Queries'});
        expect(res.status).toBe(204);
        for (let i = 1; i <= 3; i++)
          await request(seed.alice.client, 'GET', '/users/current', `${ulid()}.${i}`);
      }
      const r = await ok(['capture', 'show', id]);
      for (const label of ['rule', 'author', 'subject', 'scope', 'capture', 'reason', 'status', 'events'])
        expect(r.out).toMatch(new RegExp(`^${label}\\s+`, 'm'));
      expect(r.out).toMatch(new RegExp(`^rule\\s+${id}`, 'm'));
      if (viewRule && seed.alice) {
        const shown = await json(['capture', 'show', id]);
        expect(shown.activations.length).toBeGreaterThan(0);
        expect(r.out).toMatch(/SESSION\s+USER\s+ACTIVATED\s+UNTIL\s+TRIGGER/);
      }
    }, LONG);

    it('capture show <cap-N> --timeline --output csv > cap.csv', async () => {
      const rule = viewRule ?? fnRule;
      const id = `cap-${rule.number}`;
      const r = await ok(['capture', 'show', id, '--timeline', '--output', 'csv']);
      const file = path.join(TMP, `${id}.csv`);
      fs.writeFileSync(file, r.bytes);
      const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(Boolean);
      if (viewRule && seed.alice) {
        const rows = await until('captured requests', async () => {
          const out = (await ok(['capture', 'show', id, '--timeline', '--output', 'csv'])).out.split(/\r?\n/).filter(Boolean);
          return out.length > 1 ? out : undefined;
        }, 60_000);
        expect(rows[0]).toBe('TIME,SOURCE,SERVER,KIND,SUMMARY,STATUS,MS,REQ');
      }
      else if (lines.length) expect(lines[0]).toBe('TIME,SOURCE,SERVER,KIND,SUMMARY,STATUS,MS,REQ');
    }, LONG);

    it('timeline --rule <cap-N>  # everything the rule saw', async () => {
      const rule = viewRule ?? fnRule;
      const t = table((await ok(['timeline', '--rule', `cap-${rule.number}`])).out);
      expect(t.columns).toEqual(['TIME', 'SOURCE', 'SERVER', 'KIND', 'SUMMARY', 'STATUS', 'MS', 'REQ']);
      expect(t.rows.length).toBeGreaterThan(0);
      for (const r of t.rows) expect(r.TIME).toMatch(/^\d\d:\d\d:\d\d\.\d{3}Z$/);
      const times = t.rows.map((r) => r.TIME);
      expect([...times].sort()).toEqual(times);
    }, LONG);

    it('capture stop <cap-N> --reason "reproduced"', async () => {
      const rule = viewRule ?? fnRule;
      const id = `cap-${rule.number}`;
      const r = await ok(['capture', 'stop', id, '--reason', 'reproduced']);
      expect(r.out.trim()).toBe(`stopped ${id} — reproduced`);
      const after = await json(['capture', 'show', id]);
      expect(after.status).not.toBe('active');
      expect(table((await ok(['capture', 'list'])).out).rows.map((x) => x.RULE)).not.toContain(id);
    }, LONG);

    it('timeline --action <action id>  # one click and what it caused', async () => {
      const rows: any[] = await until('the action\'s requests', async () => {
        const r: any[] = await json(['timeline', '--action', seed.action!]);
        return r.length >= 2 ? r : undefined;
      }, 60_000);
      const ids = rows.map((e) => e.requestId).filter(Boolean);
      expect(ids).toEqual(expect.arrayContaining([`${seed.action}.1`, `${seed.action}.2`]));
      for (const id of ids) expect(String(id).startsWith(seed.action!)).toBe(true);
      const t = table((await ok(['timeline', '--action', seed.action!])).out);
      expect(t.rows.some((r) => r.REQ === `…${seed.action!.slice(-6)}.1`)).toBe(true);
    }, LONG);

    it('timeline --request <action id>.1', async () => {
      const rows: any[] = await json(['timeline', '--request', seed.request!]);
      expect(rows.length).toBeGreaterThan(0);
      for (const e of rows) if (e.requestId) expect(e.requestId).toBe(seed.request);
      expect(rows.some((e) => e.requestId === seed.request)).toBe(true);
    }, LONG);

    it.skipIf(SHARED)('timeline --session <session-id>', async () => {
      const rows: any[] = await json(['timeline', '--session', seed.alice!.session!]);
      expect(rows.length).toBeGreaterThan(0);
      expect(rows.some((e) => String(e.summary).includes(`${TAG}Solo`) || e.requestId === `${seed.action}.3`)).toBe(true);
      const t = table((await ok(['timeline', '--session', seed.alice!.session!])).out);
      expect(t.columns).toEqual(['TIME', 'SOURCE', 'SERVER', 'KIND', 'SUMMARY', 'STATUS', 'MS', 'REQ']);
    }, LONG);

    it.skipIf(SHARED)('timeline --report <number>', async () => {
      const rows: any[] = await json(['timeline', '--report', String(seed.report!.number)]);
      expect(rows.some((e) => e.requestId === `${seed.action}.3`)).toBe(true);
    }, LONG);
  });

  // ─── Tips ─────────────────────────────────────────────────────────────────

  describe('tips', () => {
    it('--output json works for scripting on each of these', async () => {
      const commands: string[][] = [
        ['alerts', 'list'], ['alerts', 'list', '--status', 'all', '--since', '1h'], ['alerts', 'detection'],
        ['errors', 'list', '--since', '2h'], ['errors', 'top', '--since', '1d', '--by', 'signature,package'],
        ['errors', 'diff', '--before', `${daysAgoUtc(1)}..${daysAgoUtc(1)}`, '--after', `${todayUtc()}..${todayUtc()}`],
        ['logger', 'get', 'server'], ['logger', 'get', '--scope', 'user:admin'], ['logger', 'overrides'],
        ['logger', 'diff'], ['logger', 'history', '--limit', '20'],
        ['capture', 'list', '--all', '--since', '7d'], ['timeline', '--action', seed.action!],
      ];
      if (seed.sigShared) commands.push(['errors', 'show', seed.sigShared.slice(0, 6), '--since', '1d']);
      if (seed.incident) commands.push(['alerts', 'get', seed.incident.id]);
      for (const c of commands) {
        const v = await json(c);
        expect(typeof v, c.join(' ')).toBe('object');
      }
    }, 600_000);

    it('times are UTC in and out', async () => {
      const r = await ok(['alerts', 'detection']);
      for (const row of table(r.out).rows) if (row['LAST SEEN']) expect(row['LAST SEEN']).toMatch(/Z$/);
      const lease = await json(['alerts', 'detection']);
      for (const s of lease.servers) expect(String(s.lastSeen)).toMatch(/Z$|\+00:00$/);
      const err = await grok(['alerts', 'mute', 'x:y', '--until', '00:00', '--reason', 'r'], {env: {TZ: 'Pacific/Kiritimati'}});
      expect(err.code).toBe(1);
      expect(err.out + err.err).toMatch(/--until 00:00 is in the past|No open alert/);
    }, LONG);

    it('a timed-out request is not repeated; a server budget refusal is shown', async () => {
      let calls = 0;
      const server = http.createServer((req, res) => {
        if (req.url?.startsWith('/users/login/dev')) {
          res.setHeader('Content-Type', 'application/json');
          res.end(JSON.stringify({token: 'Bearer stub'}));
          return;
        }
        if (req.url?.startsWith('/errors')) {
          calls++;
          if (req.url.includes('since=90d')) {
            res.statusCode = 400;
            res.setHeader('Content-Type', 'application/json');
            res.end(JSON.stringify({'#type': 'ApiError', message: 'The errors query ran over 30 s a statement or 50 s in all: narrow the window (since, from/to) or add filters', errorCode: 400}));
            return;
          }
          return;
        }
        res.statusCode = 404;
        res.end();
      });
      await new Promise<void>((resolve) => server.listen(0, '127.0.0.1', resolve));
      try {
        const port = (server.address() as any).port;
        const env = {...cliHomeFor(`http://127.0.0.1:${port}`, 'stub'), GROK_HTTP_TIMEOUT: '2000'};
        const timedOut = await grok(['errors', 'top', '--since', '7d'], {host: 'obsit', env});
        expect(timedOut.code).toBe(1);
        expect(timedOut.out + timedOut.err).toMatch(/no answer in 2000ms/);
        expect(calls).toBe(1);
        calls = 0;
        const refused = await grok(['errors', 'top', '--since', '90d'], {host: 'obsit', env});
        expect(refused.code).toBe(1);
        expect(refused.out + refused.err).toMatch(/narrow the window/);
        expect(calls).toBe(1);
      }
      finally {
        server.close();
      }
    }, LONG);

    it.skipIf(SHARED)('permissions: a user without ManageAlerts, ViewTelemetry, EditPluginsSettings is refused', async () => {
      const env = seed.bob!.home;
      const refusals: [string[], RegExp][] = [
        [['alerts', 'ack', seed.incident.id, '--reason', 'x'], /ManageAlerts/],
        [['alerts', 'mute', seed.incident.id, '--for', '1h', '--reason', 'x'], /ManageAlerts/],
        [['alerts', 'resolve', seed.incident.id, '--reason', 'x'], /ManageAlerts/],
        [['errors', 'list', '--since', '1h'], /ViewTelemetry/],
        [['errors', 'top', '--since', '1h'], /ViewTelemetry/],
        [['timeline', '--action', seed.action!], /ViewTelemetry/],
        [['logger', 'set', 'server', '--debug-flags', '+query', '--scope', `user:${seed.bob!.login}`, '--for', '5m', '--reason', 'x'], /EditPluginsSettings/],
        [['capture', 'add', '--user', seed.bob!.login, '--capture', 'requests', '--for', '5m', '--reason', 'x'], /EditPluginsSettings/],
      ];
      for (const [args, re] of refusals) {
        const r = await grok(args, {host: 'obsit', env});
        expect(r.code, args.join(' ')).toBe(1);
        expect(r.out + r.err, args.join(' ')).toMatch(re);
      }
      expect((await json(['alerts', 'get', seed.incident.id])).status).not.toBe('resolved');
    }, LONG);
  });
});
