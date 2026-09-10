/* Playwright global setup: the same storage-state login every Datagrok suite uses
   (@datagrok-libraries/test), with one convenience for local runs — when `grok test` has not
   provided DATAGROK_AUTH_TOKEN, the token is minted from the dev key in ~/.grok/config.yaml
   (server DATAGROK_SERVER, `localhost` by default). */
import {randomUUID} from 'node:crypto';
import {existsSync, mkdirSync, readFileSync, writeFileSync} from 'node:fs';
import {homedir} from 'node:os';
import {dirname, join} from 'node:path';
import {chromium, FullConfig} from '@playwright/test';
import * as libSetupModule from '@datagrok-libraries/test/src/playwright/global-setup.js';

type Setup = (config: FullConfig) => Promise<void>;
// the library is CommonJS: Node's ESM interop hands the whole `module.exports` over as `default`
const exported = libSetupModule as unknown as {default?: Setup | {default?: Setup}};
const libSetup: Setup = typeof exported.default === 'function' ? exported.default :
  (exported.default as {default?: Setup} | undefined)?.default ?? (libSetupModule as unknown as Setup);

export const DEFAULT_URL = 'http://localhost:8888';

export interface GrokServer {
  url: string;
  key: string;
}

/** The `servers:` block of ~/.grok/config.yaml (datagrok-tools); a minimal reader, the file is flat. */
export function grokServers(file = join(homedir(), '.grok', 'config.yaml')): Record<string, GrokServer> {
  if (!existsSync(file))
    return {};
  const servers: Record<string, GrokServer> = {};
  let current: string | undefined;
  for (const line of readFileSync(file, 'utf8').split(/\r?\n/)) {
    const name = /^ {2}([\w-]+):\s*$/.exec(line);
    if (name) {
      current = name[1];
      servers[current] = {url: '', key: ''};
      continue;
    }
    const field = /^ {4}(url|key):\s*(.+?)\s*$/.exec(line);
    if (field && current)
      servers[current][field[1] as 'url' | 'key'] = field[2];
  }
  return servers;
}

export async function mintToken(apiUrl: string, key: string): Promise<string> {
  const response = await fetch(`${apiUrl.replace(/\/$/, '')}/users/login/dev/${key}`,
    {method: 'POST', body: '', signal: AbortSignal.timeout(10000)});
  const json = await response.json().catch(() => ({})) as {token?: string; isSuccess?: boolean};
  if (!json.token)
    throw new Error(`dev-key login failed at ${apiUrl} (status ${response.status})`);
  return json.token;
}

/** Session cookie + localStorage after a real login-form sign-in, saved where the shared base
 * config expects the storage state — the fallback for a stand without a usable dev key
 * (DATAGROK_LOGIN / DATAGROK_PASSWORD, admin/admin by default). */
async function formLogin(url: string): Promise<void> {
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';
  const browser = await chromium.launch();
  try {
    const context = await browser.newContext();
    const page = await context.newPage();
    // a cold local client compiles on first request: the form can take a minute to appear
    await page.goto(`${url}/login.html`, {waitUntil: 'domcontentloaded', timeout: 180000});
    const loginInput = page.locator('input[placeholder="Login or Email"]:visible').first();
    await loginInput.waitFor({state: 'visible', timeout: 180000});
    await loginInput.fill(login);
    await page.locator('input[type="password"]:visible').first().fill(password);
    await page.locator('button:visible', {hasText: /^Login$/i}).first().click();
    await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
    const state = JSON.stringify(await context.storageState());
    const root = process.cwd();
    mkdirSync(join(root, 'e2e'), {recursive: true});
    for (const file of ['e2e/.auth.json', 'e2e/.auth.public.json', '.auth.json'])
      writeFileSync(join(root, file), state);
  } finally {
    await browser.close();
  }
}

export default async function globalSetup(config: FullConfig): Promise<void> {
  const url = (process.env.DATAGROK_URL ?? DEFAULT_URL).replace(/\/$/, '');
  process.env.DATAGROK_URL = url;
  // the shared setup writes e2e/.auth.json under cwd; the config reads it under the project root
  const testDir = config.projects[0]?.testDir;
  if (testDir)
    process.chdir(dirname(testDir));
  if (!process.env.DATAGROK_AUTH_TOKEN) {
    const name = process.env.DATAGROK_SERVER ?? 'localhost';
    const server = grokServers()[name];
    const candidates = [process.env.DATAGROK_API_URL, `${url}/api`, server?.url].filter((x): x is string => !!x);
    let token: string | undefined;
    let api: string | undefined;
    for (api of candidates) {
      if (!server?.key)
        break;
      token = await mintToken(api, server.key).catch(() => undefined);
      if (token)
        break;
    }
    if (!token) {
      console.log(`bdd: no dev-key token for "${name}" — signing in through the login form`);
      await formLogin(url);
      return;
    }
    process.env.DATAGROK_AUTH_TOKEN = token;
    process.env.DATAGROK_SHARING_LOGIN ??= await ensureSecondAccount(api!, token);
  }
  await warmClient(url);
  await libSetup(config);
}

export const SECOND_LOGIN = 'bddsecond';

/** The account a sharing feature shares with. DATAGROK_SHARING_LOGIN names one; without it the
 * setup makes sure a "bddsecond" user exists on the stand — a user cannot be deleted, so it is
 * created once and kept — and the workers inherit the variable pointing at it. A missing login
 * answers 200 with an ApiError body, hence the type check. */
async function ensureSecondAccount(apiUrl: string, token: string): Promise<string> {
  const base = apiUrl.replace(/\/$/, '');
  const headers = {Authorization: token, 'Content-Type': 'application/json'};
  const found = await fetch(`${base}/public/v1/users/${SECOND_LOGIN}`, {headers, signal: AbortSignal.timeout(10000)})
    .then((r) => r.json()).catch(() => ({})) as {'#type'?: string; id?: string; firstName?: string; lastName?: string};
  if (found['#type'] === 'User' && found.firstName === SECOND_LOGIN && !found.lastName)
    return SECOND_LOGIN;
  // the typeahead shows a user by name, and the step looks for the login there: the name is the
  // login. An update sends the user back whole: without its `project` the server saves the
  // personal root project under a new id and refuses the duplicate name.
  const user = {...found, '#type': 'User', id: found.id ?? randomUUID(), login: SECOND_LOGIN, firstName: SECOND_LOGIN, lastName: '', status: 'active'};
  const saved = await fetch(`${base}/public/v1/users`, {method: 'POST', headers, body: JSON.stringify(user), signal: AbortSignal.timeout(30000)});
  if (!saved.ok)
    throw new Error(`could not create the "${SECOND_LOGIN}" user the sharing features share with (status ${saved.status}): ${await saved.text()}`);
  console.log(`bdd: created the "${SECOND_LOGIN}" user for the sharing features`);
  return SECOND_LOGIN;
}

/** Asks the stand for the client bundle before a browser does. On a dev stand the Dart client is
 * served by `pub serve`, which recompiles the whole bundle after any source edit: the first
 * request then blocks for minutes while the browser's own load — and the 60 s the shell wait
 * allows — expires on a page that has not been served yet. One plain request outside the browser
 * pays that wait once, with no cap, and is a few milliseconds when nothing is compiling. A stand
 * serving a built client (a deployment, CI) answers at once and nothing is lost. */
async function warmClient(url: string): Promise<void> {
  const start = Date.now();
  for (const path of ['login.dart.js', 'login.dart.js_1.part.js'])
    await fetch(`${url}/${path}`).then((r) => r.arrayBuffer()).catch(() => undefined);
  const seconds = Math.round((Date.now() - start) / 1000);
  if (seconds >= 5)
    console.log(`bdd: the stand took ${seconds} s to serve its client (a dev stand recompiles it after a source change)`);
}
