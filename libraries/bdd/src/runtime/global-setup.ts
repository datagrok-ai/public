/* Playwright global setup: the same storage-state login every Datagrok suite uses
   (@datagrok-libraries/test), with one convenience for local runs — when `grok test` has not
   provided DATAGROK_AUTH_TOKEN, the token is minted from the dev key in ~/.grok/config.yaml
   (server DATAGROK_SERVER, `localhost` by default). */
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
    for (const api of candidates) {
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
  }
  await libSetup(config);
}
