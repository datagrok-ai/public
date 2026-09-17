import {test, Page} from '@playwright/test';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

export const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';

export const specTestOptions = {
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
  // Stated rather than inherited: specs that read an exported artefact through
  // page.waitForEvent('download') depend on it, and Playwright's default is not
  // part of any contract this section controls.
  acceptDownloads: true,
};

export interface StepError { step: string; error: string; }

/** test.step when a test is running, a plain call otherwise (worker fixtures, global setup). */
export function phase<T>(title: string, fn: () => Promise<T>): Promise<T> {
  try { return test.step(title, fn); }
  catch (_) { return fn(); }
}

export interface LedgerEntry { kind: string; ms: number; what: string; }

/**
 * Records every page.evaluate / waitForFunction / waitForTimeout the current test makes, so the
 * JSON report can say where the time outside test.step went. Attached as annotations.
 */
export function installLedger(page: Page): LedgerEntry[] {
  const ledger: LedgerEntry[] = [];
  const p = page as any;
  if (p.__ledgerInstalled) return p.__ledger;
  p.__ledgerInstalled = true;
  p.__ledger = ledger;
  const snippet = (fn: any) => String(typeof fn === 'function' ? fn.toString() : fn)
    .replace(/\s+/g, ' ').slice(0, 140);
  for (const kind of ['evaluate', 'waitForFunction', 'waitForTimeout'] as const) {
    const orig = p[kind].bind(page);
    p[kind] = async (...args: any[]) => {
      const t0 = Date.now();
      try { return await orig(...args); }
      finally { ledger.push({kind, ms: Date.now() - t0, what: kind === 'waitForTimeout' ? String(args[0]) : snippet(args[0])}); }
    };
  }
  // the Playwright-side actions were the one unmeasured bucket (435s of a 1,923s Viewers run):
  // mouse and keyboard calls, and the locator actions, are timed the same way
  const wrap = (obj: any, kind: string, names: string[], label: (name: string, args: any[]) => string) => {
    for (const name of names) {
      const orig = obj[name]?.bind(obj);
      if (!orig) continue;
      obj[name] = async (...args: any[]) => {
        const t0 = Date.now();
        try { return await orig(...args); }
        finally { ledger.push({kind, ms: Date.now() - t0, what: label(name, args)}); }
      };
    }
  };
  wrap(page.mouse, 'mouse', ['move', 'click', 'dblclick', 'down', 'up', 'wheel'],
    (n, a) => `${n} ${typeof a[0] === 'number' ? Math.round(a[0]) + ',' + Math.round(a[1]) : ''}${a[2]?.steps ? ' steps=' + a[2].steps : ''}`);
  wrap(page.keyboard, 'keyboard', ['press', 'type', 'insertText'], (n, a) => `${n} ${String(a[0]).slice(0, 30)}`);
  const origLocator = p.locator.bind(page);
  p.locator = (...args: any[]) => {
    const loc = origLocator(...args);
    wrap(loc, 'locator', ['click', 'dblclick', 'hover', 'fill', 'press', 'pressSequentially', 'waitFor', 'scrollIntoViewIfNeeded',
      'count', 'isVisible', 'textContent', 'innerText', 'boundingBox', 'evaluate', 'evaluateAll', 'check', 'selectOption'],
    (n) => `${n} ${String(args[0]).slice(0, 80)}`);
    return loc;
  };
  return ledger;
}

export function ledgerAnnotations(ledger: LedgerEntry[]): {type: string; description: string}[] {
  const out: {type: string; description: string}[] = [];
  const byKind: Record<string, {ms: number; n: number}> = {};
  for (const e of ledger ?? []) {
    byKind[e.kind] = byKind[e.kind] ?? {ms: 0, n: 0};
    byKind[e.kind].ms += e.ms; byKind[e.kind].n++;
  }
  for (const k of Object.keys(byKind)) out.push({type: 'ledger-' + k, description: `${byKind[k].ms}ms/${byKind[k].n}`});
  for (const e of [...(ledger ?? [])].sort((a, b) => b.ms - a.ms).slice(0, 12))
    out.push({type: 'ledger-top', description: `${e.ms}ms ${e.kind} ${e.what}`});
  if (ledger) ledger.length = 0;
  return out;
}

export const stepErrors: StepError[] = [];

export async function softStep(name: string, fn: () => Promise<void>): Promise<void> {
  try { await test.step(name, fn); }
  catch (e: any) {
    // `test.skip()` inside a step signals itself by throwing TestSkipError; recording
    // that as a step error reports a deliberate skip as a failed test.
    if (e?.constructor?.name === 'TestSkipError')
      throw e;
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

// Wait until the top-menu Chem entry registers — it appears only after the
// Molecule semType is detected on the active table and the Chem package is ready.
export async function waitForChemMenu(page: Page): Promise<void> {
  await page.locator('[name="div-Chem"]').first().waitFor({state: 'attached', timeout: 15_000});
}

// Resolves when a Molecule column is typed. Both the platform's detection event
// (`SemanticTypeDetector.SEMANTIC_TYPE_DETECTED`,
// grok_shared/lib/src/semantics/semantic_type_detector.dart:289) and a poll are
// armed: a column has been observed reaching semType Molecule with no further
// global event arriving, so an event-only barrier hangs on a ready table.
export async function waitForMolecule(page: Page, timeoutMs = 45_000): Promise<void> {
  await page.evaluate(({timeout}) => new Promise<void>((resolve, reject) => {
    const g = (window as any).grok;
    const typed = () => [g?.shell?.t, (window as any).__df]
      .filter(Boolean)
      .some((t: any) => t.columns.toList().some((c: any) => c.semType === 'Molecule'));

    if (typed())
      return resolve();

    const sub = g.events.onEvent('ddt-semantic-type-detected').subscribe(() => {
      if (typed())
        done();
    });
    const poll = setInterval(() => { if (typed()) done(); }, 500);
    const timer = setTimeout(
      () => done(new Error('waitForMolecule: no Molecule column detected within ' + timeout + 'ms')), timeout);

    function done(err?: Error) {
      clearTimeout(timer);
      clearInterval(poll);
      sub.unsubscribe();
      err ? reject(err) : resolve();
    }

    if (typed())
      done();
  }), {timeout: timeoutMs});
}

// A page switched to the second user must not be mistaken for an up-and-logged-in primary page.
const secondUserPages = new WeakSet<Page>();

async function injectToken(page: Page, token: string, opts: {hideTooltips?: boolean} = {}): Promise<void> {
  // Navigate to the origin first so the cookie/localStorage entries are
  // attached to the right host. The `/oauth/` path matches what `grok test`
  // does for Puppeteer (test-utils.ts:135).
  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl);
  // Cold CI Datlas keeps #grok-preloader in the DOM past the timeout even though
  // the Browse tree is already interactive (confirmed via CI page snapshot). Make
  // the wait best-effort and neutralise the preloader's click-interception so a
  // lingering preloader can neither hard-fail login nor intercept later clicks.
  await page.waitForFunction(() => document.querySelector('#grok-preloader, .grok-preloader') == null, null, {timeout: 30_000})
    .catch(() => { /* tolerate a lingering preloader — neutralised below */ });
  // Specs that assert on tooltips (the TestTrack lanes) must keep them; the rest hide them so a
  // tooltip left under the pointer cannot intercept a click.
  const css = '#grok-preloader, .grok-preloader { pointer-events: none !important; }' +
    (opts.hideTooltips === false ? '' : ' .d4-tooltip { display: none !important; }');
  await page.addStyleTag({content: css}).catch(() => {});
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

// Local mode (`?mode=local`, core/docs/features/ui2/LOCAL_MODE.md): the client boots with no
// authenticated session and answers every API call from static files. A spec whose subject is
// client behaviour runs identically there, without the token exchange, the boot round-trips or
// the per-spec dataset read — set DATAGROK_MODE=local to take that lane.
export const localMode = process.env.DATAGROK_MODE === 'local';

// Which client a page is running. The lane is declared per spec by the fixture it imports
// (`test` = server, `localTest` = local), so one run can hold both; DATAGROK_MODE=local is
// kept as a run-wide override for measuring the same spec both ways.
const lanes = new WeakMap<Page, 'local' | 'server'>();

export function setLane(page: Page, lane: 'local' | 'server'): void {
  lanes.set(page, lane);
}

/** Whether the page belongs to a shared-page lane fixture (as opposed to a plain per-test page). */
export function hasLane(page: Page): boolean {
  return lanes.has(page);
}

export function laneOf(page: Page): 'local' | 'server' {
  // DATAGROK_MODE=server forces every lane onto a real server: local mode is a fast lane for
  // building, and the discipline that goes with it is proving the same specs against a stand.
  if (process.env.DATAGROK_MODE === 'server') return 'server';
  return lanes.get(page) ?? (localMode ? 'local' : 'server');
}

// Datasets a local-mode run must not fetch from the server, mapped to the checked-in copy under
// `packages/`. `ApiTests/files/datasets/demog.csv` is byte-identical to System:DemoFiles/demog.csv,
// so category counts, tooltips and legend labels assert the same values in both lanes.
const LOCAL_DATASETS: Record<string, string> = {
  'System:DemoFiles/demog.csv': 'ApiTests/files/datasets/demog.csv',
  'System:DemoFiles/demog-1000.csv': 'ApiTests/files/datasets/demog-1000.csv',
  'System:AppData/Chem/tests/spgi-100.csv': 'UITests/files/SPGI_v2_100.csv',
};

// The lib runs from its own checkout, from a consumer's node_modules copy, or through the pnpm
// workspace link, so the repo root is searched for upwards from both the run dir and this file.
function findPackagesFile(rel: string): string | null {
  for (const start of [process.cwd(), __dirname]) {
    let dir = path.resolve(start);
    for (;;) {
      const candidate = path.join(dir, 'packages', rel);
      if (fs.existsSync(candidate)) return candidate;
      const parent = path.dirname(dir);
      if (parent === dir) break;
      dir = parent;
    }
  }
  return null;
}

let localCsvCache: Record<string, string> | null = null;

/**
 * Installs `__readCsv(path)`, the one seam a spec needs to run in either lane: on a server it is
 * `dapi.files.readAsText` + `DataFrame.fromCsv`, in local mode it parses a CSV shipped into the page
 * from LOCAL_DATASETS. A path with no local copy falls back to the server read, and throws by name
 * when that fails too, rather than resolving to an empty table.
 */
export async function installCsvBridge(page: Page): Promise<void> {
  const local = laneOf(page) === 'local';
  let texts: Record<string, string> = {};
  if (local) {
    if (!localCsvCache) {
      localCsvCache = {};
      for (const serverPath of Object.keys(LOCAL_DATASETS)) {
        const file = findPackagesFile(LOCAL_DATASETS[serverPath]);
        if (file) localCsvCache[serverPath] = fs.readFileSync(file, 'utf8');
      }
    }
    texts = localCsvCache;
  }
  await page.evaluate(({csv, local}) => {
    const w = window as any;
    w.__csv = csv;
    // a server read of demog.csv costs 1-5s on dev; the bytes do not change between the tests
    // of one worker, so the text is fetched once and parsed per call
    w.__csvText = w.__csvText ?? {};
    const serverRead = async (p: string) => {
      // a dev stall in the read must fail by name, not hold the spec to its timeout
      // (forms-spec once sat 560s in readAsText of curves.csv)
      if (!(p in w.__csvText))
        w.__csvText[p] = await Promise.race([w.grok.dapi.files.readAsText(p),
          new Promise<string>((_, rej) => setTimeout(() => rej(new Error(`readAsText("${p}") timed out after 30s`)), 30_000))]);
      return w.DG.DataFrame.fromCsv(w.__csvText[p]);
    };
    w.__readCsv = async (p: string) => {
      if (!local)
        return serverRead(p);
      if (p in w.__csv)
        return w.DG.DataFrame.fromCsv(w.__csv[p]);
      try {
        return await serverRead(p);
      }
      catch (_) {
        throw new Error(`No local copy of "${p}" — add it to LOCAL_DATASETS or run this spec on a server`);
      }
    };
  }, {csv: texts, local});
}

/**
 * Console noise a local-mode boot on dev produces that no spec caused: the deployed
 * `web/local/api.json` lists a staged package whose bundle was never copied under
 * `web/local/pkg/`, so the client 404s on it. Specs that assert a zero console-error count
 * must not be charged for it; anything else still fails them.
 */
export function isLocalBootNoise(text: string): boolean {
  // A local-client defect, not the viewer under test: grid_editors.dart:199 _initCellEditing
  // calls Node.remove on an already-detached node, raising a removeChild NotFoundError plus a
  // companion multi-line "Stack trace <id>". It fires at an arbitrary moment in ~1 run in 3 on
  // ?mode=local and never on the authenticated client, so whichever step happens to be open
  // when it lands fails its no-error floor.
  return /Failed to load resource/.test(text) || /local\/pkg\//.test(text) ||
    /removeChild.*no longer a child/.test(text) || /^Stack trace [A-Za-z0-9]+/.test(text.trim());
}

export async function openLocalDatagrok(page: Page): Promise<void> {
  const alreadyUp = await page.evaluate(() =>
    !!(window as any).grok?.shell && document.querySelector('#grok-preloader, .grok-preloader') == null,
  ).catch(() => false);
  if (alreadyUp) return;
  await page.goto(`${baseUrl}/?mode=local`);
  // A stand that does not serve local mode ignores the parameter and returns the LOGIN page
  // (verified against public.datagrok.ai): `grok.shell` exists there, so waiting on the shell
  // alone burns the full timeout and then reports a bare Playwright timeout. Race the two
  // outcomes instead and name the real cause — local mode ships in the client, so the target
  // stand has to be built from a revision that carries it.
  const outcome = await page.waitForFunction(() => {
    const w = window as any;
    if (document.querySelector('input[type="password"], .grok-login')) return 'login';
    return document.querySelector('#grok-preloader, .grok-preloader') == null && !!w.grok?.shell ? 'ready' : false;
  }, null, {timeout: 120_000}).then((h) => h.jsonValue());
  if (outcome === 'login')
    throw new Error(`${baseUrl} does not serve local mode: ?mode=local returned the login page. ` +
      'The client must be built from a revision that carries it (core/docs/features/ui2/LOCAL_MODE.md).');
  // The deployed fixture registers a debugging package, so a local boot raises a sticky
  // "Debugging packages" balloon — and a sticky balloon's container eats clicks (shared-page.ts).
  // It arrives after the shell is up, i.e. after the first spec's resetShell has already run,
  // which is what made that spec's column-selector pick fail 3 runs in 4. Wait it out here,
  // once per boot, so no spec starts under it.
  await page.waitForFunction(() => document.querySelectorAll('.d4-balloon').length > 0,
    null, {timeout: 5_000}).catch(() => {});
  await page.evaluate(() => {
    for (const b of Array.from(document.querySelectorAll('.d4-balloon'))) b.remove();
    for (const c of Array.from(document.querySelectorAll('.d4-balloon-container')))
      (c as HTMLElement).innerHTML = '';
  });
}

/** Boots whichever client the run asked for. */
export async function openDatagrok(page: Page): Promise<void> {
  await (laneOf(page) === 'local' ? openLocalDatagrok(page) : loginToDatagrok(page, {hideTooltips: false}));
  await installCsvBridge(page);
}

// The session behind DATAGROK_AUTH_TOKEN is shared by every spec, and the relogin scenarios
// POST /users/logout, which invalidates it server-side for everyone who logs in afterwards.
// Mint a throwaway session per context so one spec's logout cannot lock the rest out.
// Returns undefined whenever it cannot mint (no dev key — e.g. a keypair CI login —, no matching
// server in the config, a refused connection, a rejected key), and the caller then falls back to
// the shared DATAGROK_AUTH_TOKEN.
export async function mintToken(): Promise<string | undefined> {
  let apiUrl = process.env.DATAGROK_API_URL;
  let devKey = process.env.DATAGROK_DEV_KEY;
  if (!apiUrl || !devKey) {
    // with DATAGROK_URL set, only a config entry for that same host may be used: the config's
    // default server can be a different stand, and its token would not log in here
    const cfg = readDevKeyFromConfig('key', !!process.env.DATAGROK_URL);
    if (!cfg)
      return undefined;
    apiUrl = cfg.apiUrl;
    devKey = cfg.key;
  }
  try {
    const response = await fetch(`${apiUrl.replace(/\/$/, '')}/users/login/dev`,
      {method: 'POST', headers: {'Authorization': `Dev ${devKey}`}, signal: AbortSignal.timeout(20_000)});
    const json = await response.json().catch(() => null) as any;
    return json?.token || undefined;
  }
  catch (_) {
    return undefined;
  }
}

export async function loginToDatagrok(page: Page, opts: {hideTooltips?: boolean} = {}): Promise<void> {
  // Idempotent so a spec running on the worker-scoped booted page (shared-page.ts) can keep
  // its own login call: re-injecting would re-navigate and pay the ~10s boot this exists to
  // avoid. A page that is not up yet reports false and takes the full path.
  const alreadyUp = await page.evaluate(() =>
    !!(window as any).grok?.shell && document.querySelector('.grok-preloader') == null,
  ).catch(() => false);
  if (alreadyUp && !secondUserPages.has(page)) return;
  const token = (await mintToken()) ?? process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  await injectToken(page, token, opts);
  secondUserPages.delete(page);
}

export function fileBrowseUrl(relPath: string): string {
  return `${baseUrl}/file/${relPath}?browse=files`;
}

// Authenticate, then open the dataset into a table view. The TestTrack original
// deep-linked to the dot-namespaced file-browse URL (/file/System.AppData/...?browse=files),
// but that route 404s on the minimal ui_tests CI stack (the xamgle SPA resolves it as a
// static asset). Instead load the platform normally and open the file via the file API —
// the CI-safe pattern used across the suite (cf. Bio openBioDataset). The dot-namespaced
// connector prefix (System.AppData) is converted to the API colon form (System:AppData).
export async function loginAndOpenFile(page: Page, relPath: string): Promise<void> {
  await loginToDatagrok(page);
  const apiPath = relPath.replace(/^System\.AppData\//, 'System:AppData/');
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    const df = await g.dapi.files.readCsv(p);
    g.shell.addTableView(df);
  }, apiPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60_000});
}

// Read a dev key from ~/.grok/config.yaml for the server whose url matches the current
// DATAGROK_URL (falling back to the configured default server unless `sameHostOnly`). `field` is
// `key:` for the test user and `key2:` for the second user. js-yaml is loaded lazily: it is not a
// dependency of this lib, and without it (or without a config) there is simply no config key.
function readDevKeyFromConfig(field: 'key' | 'key2', sameHostOnly = false): {apiUrl: string; key: string} | null {
  try {
    const confPath = path.join(os.homedir(), '.grok', 'config.yaml');
    if (!fs.existsSync(confPath)) return null;
    const yaml = require('js-yaml');
    const cfg = yaml.load(fs.readFileSync(confPath, 'utf8')) as any;
    const servers = cfg?.servers ?? {};
    let wantHost: string | null = null;
    try { wantHost = new URL(baseUrl).host; }
    catch (_) { wantHost = null; }
    for (const name of Object.keys(servers)) {
      const s = servers[name];
      if (!s?.url || !s?.[field]) continue;
      let h: string | null = null;
      try { h = new URL(s.url).host; }
      catch (_) { h = null; }
      if (h && wantHost && h === wantHost)
        return {apiUrl: String(s.url).replace(/\/$/, ''), key: String(s[field])};
    }
    if (sameHostOnly) return null;
    const def = cfg?.default;
    if (def && servers[def]?.[field] && servers[def]?.url)
      return {apiUrl: String(servers[def].url).replace(/\/$/, ''), key: String(servers[def][field])};
    return null;
  }
  catch (_) {
    return null;
  }
}

async function exchangeDevKeyForToken(apiUrl: string, key: string): Promise<string> {
  const resp = await fetch(`${apiUrl}/users/login/dev`, {method: 'POST', headers: {'Authorization': `Dev ${key}`}});
  const json = await resp.json() as any;
  if (json?.isSuccess === true && json?.token) return json.token;
  throw new Error(`Second-user dev-key login failed at ${apiUrl}: ${JSON.stringify(json).slice(0, 200)}`);
}

// Resolve the second-user token: env first (DATAGROK_AUTH_TOKEN_2, which the CI runner exports after
// provisioning the `test2` user), then a second-user dev key (DATAGROK_DEV_KEY_2), then a `key2:` in
// ~/.grok/config.yaml. Throws when none is available — a two-user spec MUST NOT silently pass without
// its second user. Cached so the login claim can be read (getSecondUserLogin) without a second exchange.
let _secondTokenCache: string | null = null;
export const hasSecondUser = (): boolean =>
  !!(process.env.DATAGROK_AUTH_TOKEN_2 || process.env.DATAGROK_DEV_KEY_2 || readDevKeyFromConfig('key2'));

export async function resolveSecondUserToken(): Promise<string> {
  if (_secondTokenCache) return _secondTokenCache;
  const envTok = process.env.DATAGROK_AUTH_TOKEN_2;
  if (envTok && envTok.length > 0) return (_secondTokenCache = envTok);
  const key2 = process.env.DATAGROK_DEV_KEY_2;
  if (key2 && key2.length > 0) {
    const apiUrl = (process.env.DATAGROK_URL ?? baseUrl).replace(/\/$/, '') + '/api';
    return (_secondTokenCache = await exchangeDevKeyForToken(apiUrl, key2));
  }
  const cfg = readDevKeyFromConfig('key2');
  if (!cfg)
    throw new Error(
      'No second-user credentials available. Set DATAGROK_AUTH_TOKEN_2 / DATAGROK_DEV_KEY_2, ' +
      'or add a `key2:` (second-user dev key) to the matching server in ~/.grok/config.yaml.');
  return (_secondTokenCache = await exchangeDevKeyForToken(cfg.apiUrl, cfg.key));
}

// Resolve the second user's login so a two-user spec can learn WHO the recipient is in order to share
// the project with them. Prefer decoding the JWT `sub` (/`usr.login`) claim (no network round-trip);
// fall back to /api/users/current with the token when it isn't a decodable JWT (e.g. a plain session
// token from the CI runner's /users/login provisioning of `test2`).
export async function getSecondUserLogin(): Promise<string> {
  const token = await resolveSecondUserToken();
  try {
    const payload = token.replace(/^Bearer\s+/i, '').split('.')[1];
    if (payload) {
      const claims = JSON.parse(Buffer.from(payload, 'base64').toString('utf8'));
      const login = claims?.sub ?? claims?.usr?.login;
      if (login) return login;
    }
  }
  catch (_) { /* not a JWT — fall through to the REST lookup */ }
  const apiUrl = (process.env.DATAGROK_URL ?? baseUrl).replace(/\/$/, '') + '/api';
  const resp = await fetch(`${apiUrl}/users/current`, {headers: {Authorization: token}});
  const user = await resp.json() as any;
  if (!user?.login)
    throw new Error('Could not resolve second-user login (JWT claim and /users/current both failed)');
  return user.login;
}

export async function loginAsSecondUser(page: Page): Promise<void> {
  const token2 = await resolveSecondUserToken();
  await injectToken(page, token2);
  secondUserPages.add(page);
}
