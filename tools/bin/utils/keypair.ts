import crypto from 'crypto';
import fs from 'fs';
import os from 'os';
import path from 'path';
import yaml from 'js-yaml';
import {Config, Indexable} from './utils';

const grokDir = path.join(os.homedir(), '.grok');
const confPath = path.join(grokDir, 'config.yaml');
const keysDir = path.join(grokDir, 'keys');

/** Prefix the server wraps a nonce in before checking the signature. */
const SIGNATURE_PREFIX = 'datagrok-login:';

export interface Jwk extends Indexable { kty: string }

export interface KeyPair { publicKey: Jwk, privateKey: Jwk }

export interface ServerCredentials {
  /** API root, e.g. `https://dev.datagrok.ai/api`. */
  url: string,
  alias?: string,
  /** Developer key — deprecated, kept as a fallback while stands migrate. */
  key?: string,
  /** Private key in JWK form, from the key file or from `GROK_PRIVATE_KEY`. */
  privateKey?: Jwk,
  /** Where the private key came from, for error messages. */
  privateKeySource?: string,
  fingerprint?: string,
  login?: string,
}

/** Generates the keypair `grok login` registers: EC P-256, the JOSE ES256 curve. */
export function generateKeyPair(): KeyPair {
  const {publicKey, privateKey} = crypto.generateKeyPairSync('ec', {namedCurve: 'prime256v1'});
  return {
    publicKey: publicKey.export({format: 'jwk'}) as unknown as Jwk,
    privateKey: privateKey.export({format: 'jwk'}) as unknown as Jwk,
  };
}

/** The public half of a private JWK — what gets registered on the server. */
export function publicPart(privateKey: Jwk): Jwk {
  const {d, p, q, dp, dq, qi, ...pub} = privateKey as Indexable;
  return pub as Jwk;
}

/**
 * RFC 7638 JWK thumbprint, base64url without padding. The server computes the
 * same value from the stored key, so this is what identifies a key at login.
 */
export function fingerprint(jwk: Jwk): string {
  const members: Indexable = {
    EC: ['crv', 'kty', 'x', 'y'],
    RSA: ['e', 'kty', 'n'],
  };
  const order: string[] = members[jwk.kty];
  if (!order)
    throw new Error(`Unsupported key type "${jwk.kty}" - expected EC or RSA`);
  const canonical = '{' + order.map((m) => `${JSON.stringify(m)}:${JSON.stringify(jwk[m])}`).join(',') + '}';
  return crypto.createHash('sha256').update(canonical).digest('base64url');
}

/**
 * Signs the login nonce. [audience] is the API root this client dialed: it is part of what is
 * signed, so a server cannot relay the signature to a second stand where the same key is
 * enrolled. ECDSA signatures go out in the raw r||s form (JOSE's, not DER's), which is what
 * the server's verifier expects.
 */
export function signNonce(privateKey: Jwk, audience: string, nonce: string): string {
  const key = crypto.createPrivateKey({key: privateKey as any, format: 'jwk'});
  const options: Indexable = {key};
  if (privateKey.kty === 'EC')
    options.dsaEncoding = 'ieee-p1363';
  return crypto.sign('sha256', Buffer.from(`${SIGNATURE_PREFIX}${audience}:${nonce}`), options as any)
    .toString('base64url');
}

/** Exchanges a keypair for a session token: ask for a nonce, sign it, log in. */
export async function keyLogin(url: string, privateKey: Jwk): Promise<string> {
  const fp = fingerprint(publicPart(privateKey));
  const challenge = await postJson(`${url}/users/login/key/challenge`, {fingerprint: fp});
  const response = await postJson(`${url}/users/login/key`, {
    fingerprint: fp,
    audience: url,
    nonce: challenge.nonce,
    signature: signNonce(privateKey, url, challenge.nonce),
  });
  if (response.isSuccess !== true)
    throw new Error(response.comment ?? 'Key login failed');
  return response.token;
}

/** Registers [publicKey] with a one-shot enrollment code from the user profile. */
export async function enrollWithCode(url: string, code: string, publicKey: Jwk,
  name: string, expires?: string): Promise<Indexable> {
  return await postJson(`${url}/users/keys/enroll`,
    {code, name, expires, publicKey: JSON.stringify(publicKey)});
}

async function postJson(url: string, body: Indexable): Promise<Indexable> {
  const response = await fetch(url, {
    method: 'POST',
    headers: {'content-type': 'application/json'},
    body: JSON.stringify(body),
  });
  const text = await response.text();
  try {
    return JSON.parse(text);
  } catch {
    throw new Error(`Unexpected response from ${url} (status ${response.status}): ${text.slice(0, 200)}`);
  }
}

export function keyFilePath(alias: string): string {
  return path.join(keysDir, `${alias}.json`);
}

/** Writes the private key readable only by its owner, the way ssh-keygen does. */
export function savePrivateKey(alias: string, privateKey: Jwk): string {
  fs.mkdirSync(keysDir, {recursive: true, mode: 0o700});
  const file = keyFilePath(alias);
  fs.writeFileSync(file, JSON.stringify(privateKey, null, 2), {mode: 0o600});
  // mkdirSync/writeFileSync ignore `mode` when the path already exists.
  try {
    fs.chmodSync(keysDir, 0o700);
    fs.chmodSync(file, 0o600);
  } catch { /* Windows has no POSIX modes; ACLs already keep it in the profile. */ }
  return file;
}

export function readConfig(): Config {
  if (!fs.existsSync(confPath))
    return {servers: {}, default: ''};
  return (yaml.load(fs.readFileSync(confPath, {encoding: 'utf-8'})) as Config) ?? {servers: {}, default: ''};
}

export function writeConfig(config: Config): void {
  fs.mkdirSync(grokDir, {recursive: true});
  fs.writeFileSync(confPath, yaml.dump(config));
}

/**
 * The key at [file], or `undefined` when there is none. A config can name a key file the
 * deployment has not filled in yet — CI writes the entry and the secret separately — and that
 * has to mean "fall back to the developer key". A file that exists but cannot be read or parsed
 * is a different thing, and says so.
 */
function tryLoadKeyFile(file: string): Jwk | undefined {
  const expanded = file.startsWith('~') ? path.join(os.homedir(), file.slice(1)) : file;
  let text: string;
  try {
    text = fs.readFileSync(expanded, {encoding: 'utf-8'});
  } catch (error: any) {
    if (error?.code === 'ENOENT')
      return undefined;
    throw new Error(`cannot read the private key at ${expanded}: ${error?.message ?? error}`);
  }
  try {
    return parseKey(text);
  } catch (error: any) {
    throw new Error(`${expanded} is not a private key in JWK form: ${error?.message ?? error}`);
  }
}

/**
 * The credentials for one server: its config entry, with `GROK_PRIVATE_KEY`
 * taking precedence so a CI job can hold the key in a secret rather than on
 * disk. [hostKey] is an alias, a URL, or empty for the default server.
 */
export function getServerCredentials(hostKey: string): ServerCredentials {
  const config = readConfig();
  let host = (hostKey === '' || hostKey == null ? config.default : hostKey).trim();
  let entry: Indexable | undefined;
  let alias: string | undefined;
  let url: string;
  try {
    url = new URL(host).href;
    if (url.endsWith('/')) url = url.slice(0, -1);
    // Several aliases can name the same server. Prefer one that has a keypair: a
    // dev-key-only entry matching first would silently downgrade the login.
    const matches = Object.keys(config.servers ?? {}).filter((name) => config.servers[name].url === url);
    alias = matches.find((name) => config.servers[name].keyFile || fs.existsSync(keyFilePath(name))) ?? matches[0];
    entry = alias == null ? undefined : config.servers[alias];
  } catch (error) {
    entry = config.servers?.[host];
    if (entry == null)
      throw new Error(`Unknown server alias. Please add it to ${confPath}`);
    alias = host;
    url = entry.url;
  }

  const cred: ServerCredentials = {url: url!, alias, key: entry?.key, login: entry?.login};
  if (process.env.GROK_PRIVATE_KEY) {
    cred.privateKey = parseKey(process.env.GROK_PRIVATE_KEY);
    cred.privateKeySource = 'GROK_PRIVATE_KEY';
  }
  else if (entry?.keyFile && (cred.privateKey = tryLoadKeyFile(entry.keyFile)) != null)
    cred.privateKeySource = entry.keyFile;
  else if (alias && (cred.privateKey = tryLoadKeyFile(keyFilePath(alias))) != null)
    cred.privateKeySource = keyFilePath(alias);
  if (process.env.GROK_LOGIN)
    cred.login = process.env.GROK_LOGIN;
  return cred;
}

/** A JWK, raw or base64-encoded — CI secret stores mangle multi-line values. */
function parseKey(value: string): Jwk {
  const text = value.trim().startsWith('{') ? value : Buffer.from(value.trim(), 'base64').toString('utf-8');
  return JSON.parse(text) as Jwk;
}

export function hasKeypair(cred: ServerCredentials): boolean {
  return cred.privateKey != null;
}

/** A private JWK held in an environment variable, raw or base64-encoded. */
export function keyFromEnv(name: string): Jwk | undefined {
  const value = process.env[name];
  return value ? parseKey(value) : undefined;
}

/**
 * The private key to use for [url], or `undefined` to fall back to the developer
 * key. An explicitly supplied [devKey] that differs from the one configured for
 * this server means the caller wants *that* identity - a second CI user, say -
 * so the configured keypair must not silently take over.
 */
export function keypairFor(url: string, devKey?: string): Jwk | undefined {
  let cred: ServerCredentials;
  try {
    cred = getServerCredentials(url);
  } catch {
    return undefined;
  }
  if (cred.privateKey == null)
    return undefined;
  if (devKey && devKey !== cred.key)
    return undefined;
  return cred.privateKey;
}
