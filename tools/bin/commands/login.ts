import crypto from 'crypto';
import http from 'http';
import os from 'os';
import {spawn} from 'child_process';
import {AddressInfo} from 'net';
import * as color from '../utils/color-utils';
import {Config, Indexable} from '../utils/utils';
import * as kp from '../utils/keypair';

/**
 * `grok login <server>` - generates a keypair, registers its public half on the
 * server, and stores the private half locally. Replaces the developer key: no
 * reusable secret is copied out of the UI, and the key can be given an expiry
 * and revoked on its own.
 */
export async function login(args: LoginArgs): Promise<boolean> {
  if (args['_'].length > 2)
    return false;
  const target = (args['_'][1] ?? '').toString();

  let url: string;
  let alias: string;
  const config = kp.readConfig();
  try {
    ({url, alias} = resolveTarget(target, config, args.alias));
    url = await resolveApiRoot(url);
  } catch (error: any) {
    color.error(error.message);
    return false;
  }

  const name = args.name ?? `${os.userInfo().username}-${os.hostname()}`;
  const expires = parseExpiry(args.expires);
  if (args.expires && !expires) {
    color.error('--expires takes a number of days or an ISO date (2027-01-31)');
    return false;
  }

  const {publicKey, privateKey} = kp.generateKeyPair();
  console.log(`Registering key "${name}" on ${url}`);

  let registered: Indexable;
  try {
    registered = args.code
      ? await kp.enrollWithCode(url, args.code, publicKey, name, expires)
      : await enrollInBrowser(url, publicKey, name, expires);
  } catch (error: any) {
    color.error(error.message ?? String(error));
    return false;
  }
  if (registered.login == null) {
    color.error(registered.comment ?? registered.message ?? 'The server did not accept the key');
    return false;
  }

  const keyFile = kp.savePrivateKey(alias, privateKey);
  config.servers ??= {};
  config.servers[alias] = {...(config.servers[alias] ?? {}), url, key: config.servers[alias]?.key ?? '',
    keyFile, login: registered.login};
  config.default ??= alias;
  kp.writeConfig(config);

  // Proves the round trip before the user walks away, rather than at the next
  // `grok publish`: a key that registered but cannot sign is worse than none.
  try {
    await kp.keyLogin(url, privateKey);
  } catch (error: any) {
    color.error(`The key was registered but the test login failed: ${error.message ?? error}`);
    return false;
  }

  color.success(`Logged in to ${url} as ${registered.login}`);
  console.log(`  key         ${name} (${kp.fingerprint(publicKey)})`);
  console.log(`  private key ${keyFile}`);
  console.log(`  alias       ${alias} - use it as \`grok publish ${alias}\`, \`grok test --host ${alias}\``);
  if (expires)
    console.log(`  expires     ${expires}`);
  return true;
}

function resolveTarget(target: string, config: Config, aliasArg?: string): {url: string, alias: string} {
  if (target === '') {
    const alias = aliasArg ?? config.default;
    if (!alias || !config.servers?.[alias])
      throw new Error('Which server? Pass a URL or a configured alias: grok login https://dev.datagrok.ai/api');
    return {url: config.servers[alias].url, alias};
  }
  const configured = config.servers?.[target];
  if (configured)
    return {url: configured.url, alias: aliasArg ?? target};
  let parsed: URL;
  try {
    parsed = new URL(target);
  } catch {
    throw new Error(`"${target}" is neither a URL nor a server in your config`);
  }
  return {url: parsed.href.replace(/\/$/, ''), alias: aliasArg ?? defaultAlias(parsed)};
}

/**
 * The API base for [url]. A stand behind nginx serves it at `<origin>/api`, a bare Datlas at
 * the origin itself, and there is no telling which from the URL alone - so ask, rather than
 * guess and fail at the first call.
 *
 * A 200 is not the answer: nginx serves the single-page app for anything it does not route,
 * so the origin of a real stand answers `/info/server` with the app's HTML. Only a JSON body
 * that names the server counts.
 */
async function resolveApiRoot(url: string): Promise<string> {
  const candidates = /\/api$/.test(url) ? [url] : [`${url}/api`, url];
  for (const candidate of candidates) {
    try {
      const response = await fetch(`${candidate}/info/server`);
      if (!response.ok)
        continue;
      const info = JSON.parse(await response.text());
      if (info?.webRoot != null || info?.Version != null)
        return candidate;
    } catch { /* not JSON, or unreachable: try the next shape */ }
  }
  throw new Error(`${url} does not answer as a Datagrok API (tried ${candidates.join(' and ')})`);
}

function defaultAlias(url: URL): string {
  const host = url.hostname;
  if (host === 'localhost' || host === '127.0.0.1' || host === '::1')
    return 'local';
  return host.split('.')[0];
}

/** Days from now, or an ISO date, as the ISO instant the server stores. */
function parseExpiry(value?: string | number): string | undefined {
  if (value == null || value === '')
    return undefined;
  const days = Number(value);
  if (!isNaN(days) && days > 0)
    return new Date(Date.now() + days * 86400000).toISOString();
  const date = new Date(String(value));
  return isNaN(date.getTime()) ? undefined : date.toISOString();
}

/**
 * Opens the server's enrollment page in a browser and waits for it to call back
 * on a loopback listener. The user authenticates however that stand does -
 * password, SSO, SAML - because the approval happens in the platform UI.
 */
async function enrollInBrowser(url: string, publicKey: kp.Jwk, name: string,
  expires?: string): Promise<Indexable> {
  const origin = new URL(url).origin;
  // The page asks for this before it registers anything. It exists only in this terminal, so a
  // link someone else sent has nothing for the user to type — which is the difference between
  // approving one's own `grok login` and handing an attacker a key to one's account.
  const verify = crypto.randomBytes(4).toString('hex').toUpperCase().slice(0, 6);
  return await new Promise<Indexable>((resolve, reject) => {
    const timeout = setTimeout(() => {
      server.close();
      reject(new Error('Timed out waiting for the browser. Use `grok login <server> --code <code>` ' +
        'with a code from your profile page instead.'));
    }, 5 * 60 * 1000);

    const server = http.createServer((req, res) => {
      const params = new URL(req.url ?? '/', 'http://127.0.0.1').searchParams;
      const status = params.get('status');
      res.writeHead(200, {'content-type': 'text/html; charset=utf-8'});
      res.end(status === 'ok'
        ? '<h3>Key registered. You can close this tab and return to the terminal.</h3>'
        : `<h3>Key registration was cancelled.</h3><p>${escapeHtml(params.get('message') ?? '')}</p>`);
      clearTimeout(timeout);
      server.close();
      if (status === 'ok')
        resolve({login: params.get('login'), fingerprint: params.get('fingerprint')});
      else
        reject(new Error(params.get('message') ?? 'Key registration was cancelled in the browser'));
    });

    server.listen(0, '127.0.0.1', () => {
      const port = (server.address() as AddressInfo).port;
      const enrollUrl = `${origin}/enroll-key?` + new URLSearchParams({
        pk: JSON.stringify(publicKey),
        name,
        ...(expires ? {expires} : {}),
        verify,
        cb: `http://127.0.0.1:${port}`,
      }).toString();
      console.log(`Opening ${origin} in your browser to approve the key...`);
      console.log(`If it does not open, visit:\n  ${enrollUrl}`);
      console.log(`\n  Verification code: ${verify}\n`);
      openBrowser(enrollUrl);
    });
  });
}

function escapeHtml(s: string): string {
  return s.replace(/[&<>"]/g, (c) => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c]!);
}

function openBrowser(url: string): void {
  try {
    // cmd splits an unquoted argument at `&`, and the enrollment URL is all query parameters;
    // rundll32 takes the whole thing and hands it to the default browser.
    const [command, args] = process.platform === 'win32'
      ? ['rundll32', ['url.dll,FileProtocolHandler', url]]
      : process.platform === 'darwin' ? ['open', [url]] : ['xdg-open', [url]];
    spawn(command, args, {detached: true, stdio: 'ignore'}).unref();
  } catch {
    // The URL is printed above; a headless box just uses that.
  }
}

interface LoginArgs {
  _: string[],
  /** One-shot enrollment code from the profile page; skips the browser step. */
  code?: string,
  /** Key name shown in the profile. Defaults to user@host. */
  name?: string,
  /** Days from now, or an ISO date. Without it the key never expires. */
  expires?: string | number,
  /** Config alias to write. Defaults to the server's first host label. */
  alias?: string,
}
