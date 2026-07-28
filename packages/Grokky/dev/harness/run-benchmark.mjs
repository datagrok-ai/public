#!/usr/bin/env node
// Drives Grokky:runBenchmark in a real browser against a running Datagrok, then writes the report
// into files/benchmarks/ in the repo.
//
// The benchmark cannot run headlessly like the T2 driver: its asserts read live DG objects
// (viewers, filters, grid state) that only exist in a browser tab. So this logs in with the dev
// key exactly the way `grok test` does, evaluates the function in page context, and polls until it
// returns. A full arm is ~45 min, well past any default timeout, hence the explicit long waits.
//
//   node run-benchmark.mjs --label sonnet-baseline --reps 3 --model sonnet
//   node run-benchmark.mjs --label smoke --reps 1 --only trivial,visualization
//   node run-benchmark.mjs --compare sonnet-baseline,haiku-probe
//
// --headed shows the browser; --host picks a ~/.grok/config.yaml server alias (default: localhost).

import {execFileSync} from 'node:child_process';
import * as fs from 'node:fs';
import * as path from 'node:path';
import {chromium} from 'playwright';
import yaml from 'yaml';

const HARNESS = import.meta.dirname;
const OUT_DIR = path.join(HARNESS, '..', '..', 'files', 'benchmarks');
// A rep can legitimately sit at the 120 s turn timeout, and an arm is 44 prompts x reps of those.
const POLL_MS = 15000;
const MAX_RUN_MS = 3 * 60 * 60 * 1000;

function arg(name, fallback = undefined) {
  const i = process.argv.indexOf(`--${name}`);
  return i >= 0 && process.argv[i + 1] && !process.argv[i + 1].startsWith('--') ?
    process.argv[i + 1] : fallback;
}
const flag = (name) => process.argv.includes(`--${name}`);

function serverConfig(alias) {
  const cfgPath = path.join(process.env.USERPROFILE ?? process.env.HOME, '.grok', 'config.yaml');
  const cfg = yaml.parse(fs.readFileSync(cfgPath, 'utf8'));
  const server = cfg.servers?.[alias ?? cfg.default];
  if (!server) throw new Error(`No server '${alias ?? cfg.default}' in ${cfgPath}`);
  return {apiUrl: server.url.replace(/\/$/, ''), key: server.key};
}

async function api(apiUrl, route, token) {
  const r = await fetch(apiUrl + route, {headers: token ? {Authorization: token} : {}});
  if (!r.ok) throw new Error(`${route} -> ${r.status}`);
  return r.json();
}

/** Logs in via the dev key and returns a page on a fully loaded client. */
async function openClient(browser, apiUrl, token, webUrl) {
  const ctx = await browser.newContext({viewport: {width: 1600, height: 1000}});
  const page = await ctx.newPage();
  // Cookie + localStorage on the client origin, mirroring global-setup.ts in @datagrok-libraries/test.
  await page.goto(webUrl + '/oauth/', {waitUntil: 'domcontentloaded', timeout: 120000});
  await ctx.addCookies([{name: 'auth', value: token, domain: new URL(webUrl).hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(webUrl, {waitUntil: 'domcontentloaded', timeout: 120000});
  await page.waitForFunction(() => document.querySelector('#grok-preloader, .grok-preloader') == null,
    null, {timeout: 120000}).catch(() => {});
  // grok.functions is the actual gate — the DOM can be up before the JS API is bound.
  await page.waitForFunction(() => !!(window.grok && window.grok.functions && window.grok.shell),
    null, {timeout: 120000});
  // Tooltips intercept nothing here, but a stray dialog would; keep the surface quiet.
  await page.addStyleTag({content: '.d4-tooltip{display:none!important}'}).catch(() => {});
  return page;
}

/** Ensures the on-demand claude-runtime container is actually running.
 * Publishing the package removes it, and the client's first request then fails outright with
 * "Claude runtime container is not running" rather than waiting for a spawn — so a run started
 * straight after a publish dies at turn zero. Starting it explicitly makes that deterministic. */
async function ensureRuntimeContainer(page) {
  const status = await page.evaluate(async () => {
    const grok = window.grok;
    const containers = await grok.dapi.docker.dockerContainers.list();
    const c = containers.find((x) => /claude-runtime/.test(x.name ?? ''));
    if (!c) return {ok: false, detail: 'no claude-runtime container is registered'};
    // Always call run() rather than trusting the recorded status: the platform's status is its
    // own bookkeeping, and it reads "started" for a container that has since been removed
    // underneath it (an image rebuild, a docker rm, a crash). run() is idempotent and, with
    // awaitStart, returns only once the container is genuinely up.
    try {
      await grok.dapi.docker.dockerContainers.run(c.id, true);
    } catch (e) {
      return {ok: false, detail: `${c.name} failed to start: ${e?.message ?? e}`};
    }
    return {ok: true, detail: `${c.name} up (status was ${c.status})`};
  });
  console.log(`runtime: ${status.detail}`);
  if (!status.ok) throw new Error(status.detail);

  // A container that has just been (re)spawned comes up with the Dockerfile's empty `{}`
  // credentials, and an unauthenticated runtime answers every turn in ~1 s with zero tokens and
  // no error — a whole arm of plausible-looking failures. Cheap to prevent, expensive to debug.
  try {
    const out = execFileSync('node', [path.join(HARNESS, '..', 'seed-platform-creds.mjs')],
      {encoding: 'utf8'});
    for (const line of out.split('\n').filter((l) => /seeded/.test(l)))
      console.log(`creds: ${line.trim()}`);
  } catch (e) {
    console.warn(`creds: could not verify (${e.message.split('\n')[0]}) — continuing`);
  }
}

/** Starts the call in page context and polls, so a 45-minute arm never rides on one evaluate(). */
async function callLongRunning(page, fnName, params) {
  await page.evaluate(({fnName, params}) => {
    const w = window;
    w.__bench = {done: false, error: null, result: null};
    w.grok.functions.call(fnName, params)
      .then((r) => {w.__bench.result = r; w.__bench.done = true;})
      .catch((e) => {w.__bench.error = String(e?.message ?? e); w.__bench.done = true;});
  }, {fnName, params});

  const t0 = Date.now();
  let lastNote = '';
  while (Date.now() - t0 < MAX_RUN_MS) {
    await page.waitForTimeout(POLL_MS);
    const state = await page.evaluate(() => {
      const b = window.__bench ?? {};
      // The harness reports progress through balloons; surface the newest as a heartbeat.
      const notes = Array.from(document.querySelectorAll('.d4-balloon-content, .d4-balloon'));
      return {done: b.done, error: b.error, result: b.result,
        note: notes.length ? notes[notes.length - 1].textContent.trim().slice(0, 120) : ''};
    }).catch((e) => ({done: false, error: null, result: null, note: `(page busy: ${e.message})`}));
    if (state.note && state.note !== lastNote) {
      lastNote = state.note;
      console.log(`  [${Math.round((Date.now() - t0) / 1000)}s] ${state.note}`);
    }
    if (state.done) {
      if (state.error) throw new Error(state.error);
      return state.result;
    }
  }
  throw new Error(`timed out after ${Math.round(MAX_RUN_MS / 60000)} min`);
}

/** Pulls the report the run persisted to AppData into the repo, so results are reviewable in git. */
async function saveReports(page, names) {
  fs.mkdirSync(OUT_DIR, {recursive: true});
  const saved = [];
  for (const name of names) {
    const text = await page.evaluate(async (p) => {
      try {
        return await window.grok.dapi.files.readAsText(p);
      } catch {
        return null;
      }
    }, `System:AppData/Grokky/benchmarks/${name}`);
    if (text == null) {
      console.warn(`  ! ${name} not found in AppData`);
      continue;
    }
    const dest = path.join(OUT_DIR, name);
    fs.writeFileSync(dest, text);
    saved.push(dest);
  }
  return saved;
}

(async () => {
  const {apiUrl, key} = serverConfig(arg('host'));
  const {token} = await (await fetch(`${apiUrl}/users/login/dev/${key}`, {method: 'POST'})).json();
  if (!token) throw new Error(`dev-key login failed against ${apiUrl}`);
  const webUrl = (await api(apiUrl, '/admin/plugins/admin/settings', token)).settings.webRoot.replace(/\/$/, '');
  console.log(`server ${apiUrl} · client ${webUrl}`);

  const browser = await chromium.launch({headless: !flag('headed')});
  try {
    const page = await openClient(browser, apiUrl, token, webUrl);
    page.on('pageerror', (e) => console.warn('  page error:', String(e).slice(0, 160)));
    console.log('logged in, client ready');

    const compare = arg('compare');
    if (!compare)
      await ensureRuntimeContainer(page);
    if (compare) {
      console.log(`comparing: ${compare}`);
      await callLongRunning(page, 'Grokky:compareBenchmarks', {labels: compare});
      const name = `benchmark-compare-${compare.split(',').map((s) => s.trim()).join('-vs-')}.md`;
      const saved = await saveReports(page, [name]);
      console.log(saved.length ? `saved ${saved.join(', ')}` : 'comparison not retrieved');
      return;
    }

    const label = arg('label') ?? 'run';
    const reps = Number(arg('reps', '3'));
    const model = arg('model');
    const only = arg('only');
    console.log(`run "${label}" · reps ${reps} · model ${model ?? 'default'}${only ? ` · only ${only}` : ''}`);
    const t0 = Date.now();
    const msg = await callLongRunning(page, 'Grokky:runBenchmark',
      {label, reps, ...(model ? {model} : {}), ...(only ? {only} : {})});
    console.log(`\n${msg}`);
    console.log(`wall clock: ${Math.round((Date.now() - t0) / 60000)} min`);
    const saved = await saveReports(page, [`benchmark-${label}.json`, `benchmark-${label}.md`]);
    for (const s of saved) console.log(`saved ${path.relative(process.cwd(), s)}`);
  } finally {
    await browser.close();
  }
})().catch((e) => {
  console.error('FAILED:', e.message);
  process.exit(1);
});
