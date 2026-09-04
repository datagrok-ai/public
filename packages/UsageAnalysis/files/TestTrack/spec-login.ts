import {test, Page} from '@playwright/test';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import * as yaml from 'js-yaml';

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

export const stepErrors: StepError[] = [];

export async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

export async function waitForChemMenu(page: Page) {
  await page.locator('[name="div-Chem"]').first().waitFor({state: 'attached', timeout: 15_000});
}

// Resolves when a Molecule column is typed. Both the platform's detection event
// (`SemanticTypeDetector.SEMANTIC_TYPE_DETECTED`,
// grok_shared/lib/src/semantics/semantic_type_detector.dart:289) and a poll are
// armed: a column has been observed reaching semType Molecule with no further
// global event arriving, so an event-only barrier hangs on a ready table.
export async function waitForMolecule(page: Page, timeoutMs = 45_000) {
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

async function injectToken(page: Page, token: string) {

  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl);
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
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

export function setLane(page: Page, lane: 'local' | 'server') {
  lanes.set(page, lane);
}

export function laneOf(page: Page): 'local' | 'server' {
  // DATAGROK_MODE=server forces every lane onto a real server: local mode is a fast lane for
  // building, and the discipline that goes with it is proving the same specs against a stand.
  if (process.env.DATAGROK_MODE === 'server') return 'server';
  return lanes.get(page) ?? (localMode ? 'local' : 'server');
}

// Datasets a local-mode run must not fetch from the server, mapped to the checked-in copy.
// `packages/ApiTests/files/datasets/demog.csv` is byte-identical to System:DemoFiles/demog.csv,
// so category counts, tooltips and legend labels assert the same values in both lanes.
const LOCAL_DATASETS: Record<string, string> = {
  'System:DemoFiles/demog.csv': '../../../ApiTests/files/datasets/demog.csv',
  'System:DemoFiles/demog-1000.csv': '../../../ApiTests/files/datasets/demog-1000.csv',
  'System:AppData/Chem/tests/spgi-100.csv': '../../../UITests/files/SPGI_v2_100.csv',
};

/**
 * Installs `__readCsv(path)`, the one seam a spec needs to run in either lane: on a server it is
 * `dapi.files.readCsv`, in local mode it parses a CSV shipped into the page from LOCAL_DATASETS.
 * An unmapped path throws by name rather than resolving to an empty table, which is how local
 * mode degrades everything server-backed.
 */
let localCsvCache: Record<string, string> | null = null;

export async function installCsvBridge(page: Page) {
  const local = laneOf(page) === 'local';
  let texts: Record<string, string> = {};
  if (local) {
    if (!localCsvCache) {
      localCsvCache = {};
      for (const serverPath of Object.keys(LOCAL_DATASETS)) {
        const file = path.resolve(__dirname, LOCAL_DATASETS[serverPath]);
        if (fs.existsSync(file)) localCsvCache[serverPath] = fs.readFileSync(file, 'utf8');
      }
    }
    texts = localCsvCache;
  }
  await page.evaluate(({csv, local}) => {
    const w = window as any;
    w.__csv = csv;
    w.__readCsv = async (p: string) => {
      if (!local) return w.grok.dapi.files.readCsv(p);
      if (!(p in w.__csv))
        throw new Error(`No local copy of "${p}" — add it to LOCAL_DATASETS or run this spec on a server`);
      return w.DG.DataFrame.fromCsv(w.__csv[p]);
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

export async function openLocalDatagrok(page: Page) {
  const alreadyUp = await page.evaluate(() =>
    !!(window as any).grok?.shell && document.querySelector('.grok-preloader') == null,
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
    return document.querySelector('.grok-preloader') == null && !!w.grok?.shell ? 'ready' : false;
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
export async function openDatagrok(page: Page) {
  await (laneOf(page) === 'local' ? openLocalDatagrok(page) : loginToDatagrok(page));
  await installCsvBridge(page);
}

export async function loginToDatagrok(page: Page) {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  // Idempotent so a spec running on the worker-scoped booted page (shared-page.ts) can keep
  // its own login call: re-injecting would re-navigate and pay the ~10s boot this exists to
  // avoid. A page that is not up yet reports false and takes the full path.
  const alreadyUp = await page.evaluate(() =>
    !!(window as any).grok?.shell && document.querySelector('.grok-preloader') == null,
  ).catch(() => false);
  if (alreadyUp) return;
  await injectToken(page, token);
}

export function fileBrowseUrl(relPath: string): string {
  return `${baseUrl}/file/${relPath}?browse=files`;
}

export async function loginAndOpenFile(page: Page, relPath: string) {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');

  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);

  await page.goto(fileBrowseUrl(relPath));
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60_000});
}

function readSecondUserDevKeyFromConfig(): {apiUrl: string; key2: string} | null {
  try {
    const confPath = path.join(os.homedir(), '.grok', 'config.yaml');
    if (!fs.existsSync(confPath)) return null;
    const cfg = yaml.load(fs.readFileSync(confPath, 'utf8')) as any;
    const servers = cfg?.servers ?? {};
    let wantHost: string | null = null;
    try { wantHost = new URL(baseUrl).host; } catch (_) { wantHost = null; }
    for (const name of Object.keys(servers)) {
      const s = servers[name];
      if (!s?.url || !s?.key2) continue;
      let h: string | null = null;
      try { h = new URL(s.url).host; } catch (_) { h = null; }
      if (h && wantHost && h === wantHost)
        return {apiUrl: String(s.url).replace(/\/$/, ''), key2: String(s.key2)};
    }
    const def = cfg?.default;
    if (def && servers[def]?.key2 && servers[def]?.url)
      return {apiUrl: String(servers[def].url).replace(/\/$/, ''), key2: String(servers[def].key2)};
    return null;
  } catch (_) {
    return null;
  }
}

async function exchangeDevKeyForToken(apiUrl: string, key: string): Promise<string> {
  const resp = await fetch(`${apiUrl}/users/login/dev/${key}`, {method: 'POST'});
  const json = await resp.json() as any;
  if (json?.isSuccess === true && json?.token) return json.token;
  throw new Error(`Second-user dev-key login failed at ${apiUrl}: ${JSON.stringify(json).slice(0, 200)}`);
}

let _secondTokenCache: string | null = null;
export async function resolveSecondUserToken(): Promise<string> {
  if (_secondTokenCache) return _secondTokenCache;
  const envTok = process.env.DATAGROK_AUTH_TOKEN_2;
  if (envTok && envTok.length > 0) return (_secondTokenCache = envTok);
  const cfg = readSecondUserDevKeyFromConfig();
  if (!cfg)
    throw new Error(
      'No second-user credentials available. Set DATAGROK_AUTH_TOKEN_2 / DATAGROK_DEV_KEY_2, ' +
      'or add a `key2:` (second-user dev key) to the matching server in ~/.grok/config.yaml.');
  return (_secondTokenCache = await exchangeDevKeyForToken(cfg.apiUrl, cfg.key2));
}

export async function getSecondUserLogin(): Promise<string> {
  const token = await resolveSecondUserToken();
  try {
    const payload = token.replace(/^Bearer\s+/i, '').split('.')[1];
    const claims = JSON.parse(Buffer.from(payload, 'base64').toString('utf8'));
    const login = claims?.sub ?? claims?.usr?.login;
    if (!login) throw new Error('no sub/usr.login claim');
    return login;
  } catch (e: any) {
    throw new Error(`Could not read second-user login from token claim: ${e?.message ?? e}`);
  }
}

export async function loginAsSecondUser(page: Page) {
  const token2 = await resolveSecondUserToken();
  await injectToken(page, token2);
}
