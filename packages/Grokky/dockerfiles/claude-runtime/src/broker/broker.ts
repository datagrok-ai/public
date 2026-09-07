import * as http from 'node:http';
import * as https from 'node:https';
import * as fs from 'node:fs';
import {spawn, ChildProcess} from 'node:child_process';
import {randomUUID} from 'node:crypto';
import {URL} from 'node:url';
import aws4 from 'aws4';
import {BROKER_HOST, BROKER_PORT, BrokerMode} from './provider-env';

const OAUTH_BETA = 'oauth-2025-04-20';
const OAUTH_CLIENT_ID = '9d1c250a-e61b-44d9-88ed-5944d1962f5e';
const OAUTH_TOKEN_URL = 'https://platform.claude.com/v1/oauth/token';
const AUTH_HEADERS = ['authorization', 'x-api-key', 'api-key'];
const AUTH_LOGIN_TIMEOUT_MS = 2 * 60_000;
const MAX_MCP_TARGETS = 2000;

interface McpTarget {targetUrl: string; apiKey?: string; apiUrl?: string}
interface OAuthCreds {accessToken?: string; refreshToken?: string; expiresAt?: number}

export interface BrokerConfig {
  mode: BrokerMode;
  credentialsPath: string;
  anthropicKey?: string;
  awsBearerToken?: string;
  awsAccessKeyId?: string;
  awsSecretAccessKey?: string;
  awsSessionToken?: string;
  region?: string;
  foundryApiKey?: string;
  foundryResource?: string;
  models: {opus?: string; sonnet?: string; haiku?: string};
  upstreams: {anthropic: string; bedrock?: string; foundry?: string};
}

export function resolveConfig(env: NodeJS.ProcessEnv = process.env): BrokerConfig {
  const provider = env['provider'] || 'Anthropic';
  const credentialsPath = env['CREDENTIALS_PATH'] || '/home/broker/.claude/.credentials.json';

  let mode: BrokerMode;
  if (provider === 'Bedrock')
    mode = 'bedrock';
  else if (provider === 'Microsoft Foundry')
    mode = 'foundry';
  else if (env['apiKey'])
    mode = 'anthropic';
  else
    mode = 'subscription';

  return {
    mode,
    credentialsPath,
    anthropicKey: env['apiKey'],
    awsBearerToken: env['awsBearerToken'],
    awsAccessKeyId: env['awsAccessKeyId'],
    awsSecretAccessKey: env['awsSecretAccessKey'],
    awsSessionToken: env['awsSessionToken'],
    region: env['region'],
    foundryApiKey: env['foundryApiKey'],
    foundryResource: env['foundryResource'],
    models: {opus: env['opusModel'], sonnet: env['sonnetModel'], haiku: env['haikuModel']},
    upstreams: {
      anthropic: env['ANTHROPIC_UPSTREAM'] || 'https://api.anthropic.com',
      bedrock: env['BEDROCK_UPSTREAM'] ||
        (env['region'] ? `https://bedrock-runtime.${env['region']}.amazonaws.com` : undefined),
      foundry: env['FOUNDRY_UPSTREAM'] ||
        (env['foundryResource'] ? `https://${env['foundryResource']}.services.ai.azure.com/anthropic` : undefined),
    },
  };
}

function readSubscription(path: string): OAuthCreds | undefined {
  try {
    const j = JSON.parse(fs.readFileSync(path, 'utf8'));
    const o = j?.claudeAiOauth ?? j;
    return o?.accessToken ? o : undefined;
  } catch {
    return undefined;
  }
}

function writeSubscription(path: string, creds: OAuthCreds): void {
  fs.writeFileSync(path, JSON.stringify({claudeAiOauth: creds}), {mode: 0o600});
}

function send(res: http.ServerResponse, status: number, body: unknown): void {
  res.writeHead(status, {'content-type': 'application/json'});
  res.end(JSON.stringify(body));
}

function sendAuthError(res: http.ServerResponse, message: string): void {
  res.writeHead(401, {'content-type': 'application/json'});
  res.end(JSON.stringify({type: 'error', error: {type: 'authentication_error', message}}));
}

function readBody(req: http.IncomingMessage): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    const chunks: Buffer[] = [];
    req.on('data', (c) => chunks.push(c as Buffer));
    req.on('end', () => resolve(Buffer.concat(chunks)));
    req.on('error', reject);
  });
}

function baseHeaders(
  req: http.IncomingMessage, target: URL, inject: Record<string, string>,
): Record<string, string | string[]> {
  const headers: Record<string, string | string[]> = {};
  for (const [k, v] of Object.entries(req.headers)) {
    if (v !== undefined && k !== 'host' && k !== 'content-length' && !AUTH_HEADERS.includes(k))
      headers[k] = v;
  }
  Object.assign(headers, inject);
  headers['host'] = target.host;
  return headers;
}

function sendUpstream(
  target: URL, method: string, headers: http.OutgoingHttpHeaders, body: Buffer,
): Promise<http.IncomingMessage> {
  return new Promise((resolve, reject) => {
    const lib = target.protocol === 'https:' ? https : http;
    const upstream = lib.request({
      protocol: target.protocol,
      hostname: target.hostname,
      port: target.port || (target.protocol === 'https:' ? 443 : 80),
      method,
      path: target.pathname + target.search,
      headers,
    }, resolve);
    upstream.on('error', reject);
    if (body.length)
      upstream.write(body);
    upstream.end();
  });
}

interface ProxyOpts {
  sign?: {accessKeyId: string; secretAccessKey: string; sessionToken?: string; region: string};
  onUnauthorized?: () => Promise<Record<string, string> | null>;
}

async function proxy(
  req: http.IncomingMessage, res: http.ServerResponse, targetUrl: string,
  inject: Record<string, string>, opts: ProxyOpts = {},
): Promise<void> {
  const target = new URL(targetUrl);
  const body = await readBody(req);
  const method = req.method || 'GET';

  const attempt = (headersInject: Record<string, string>) => {
    let headers = baseHeaders(req, target, headersInject);
    if (opts.sign) {
      const signed = aws4.sign({
        host: target.host, method, path: target.pathname + target.search,
        service: 'bedrock', region: opts.sign.region, headers, body,
      }, opts.sign);
      headers = signed.headers as Record<string, string | string[]>;
    }
    return sendUpstream(target, method, headers, body);
  };

  try {
    let up = await attempt(inject);
    if (up.statusCode === 401 && opts.onUnauthorized) {
      const retryInject = await opts.onUnauthorized();
      if (retryInject) {
        up.resume();
        up = await attempt(retryInject);
      }
    }
    res.writeHead(up.statusCode || 502, up.headers);
    up.pipe(res);
  } catch (e: any) {
    console.error('[broker] upstream error:', e.message);
    if (!res.headersSent)
      send(res, 502, {error: 'broker upstream error'});
    else
      res.end();
  }
}

async function refreshSubscription(cfg: BrokerConfig): Promise<string | null> {
  const creds = readSubscription(cfg.credentialsPath);
  if (!creds?.refreshToken)
    return null;
  try {
    const resp = await fetch(OAUTH_TOKEN_URL, {
      method: 'POST',
      headers: {'content-type': 'application/json'},
      body: JSON.stringify({
        grant_type: 'refresh_token', refresh_token: creds.refreshToken, client_id: OAUTH_CLIENT_ID,
      }),
    });
    if (!resp.ok)
      return null;
    const j: any = await resp.json();
    const next: OAuthCreds = {
      accessToken: j.access_token,
      refreshToken: j.refresh_token || creds.refreshToken,
      expiresAt: j.expires_in ? Date.now() + j.expires_in * 1000 : undefined,
    };
    writeSubscription(cfg.credentialsPath, next);
    return next.accessToken || null;
  } catch (e: any) {
    console.error('[broker] subscription refresh failed:', e.message);
    return null;
  }
}

function mergeBeta(existing: string | string[] | undefined): string {
  const parts = new Set<string>();
  const raw = Array.isArray(existing) ? existing.join(',') : (existing || '');
  for (const part of raw.split(',')) {
    if (part.trim())
      parts.add(part.trim());
  }
  parts.add(OAUTH_BETA);
  return [...parts].join(',');
}

abstract class Provider {
  constructor(protected cfg: BrokerConfig) {}
  abstract readonly mode: BrokerMode;
  abstract hasCredential(): boolean;
  authRequired(): boolean { return false; }
  abstract forward(
    req: http.IncomingMessage, res: http.ServerResponse, pathname: string, search: string): Promise<void>;
}

class AnthropicProvider extends Provider {
  readonly mode = 'anthropic' as const;
  hasCredential(): boolean { return !!this.cfg.anthropicKey; }
  async forward(req: http.IncomingMessage, res: http.ServerResponse, pathname: string, search: string): Promise<void> {
    if (!this.cfg.anthropicKey)
      return send(res, 500, {error: 'no anthropic credential configured'});
    return proxy(req, res, this.cfg.upstreams.anthropic + pathname + search, {'x-api-key': this.cfg.anthropicKey});
  }
}

class SubscriptionProvider extends Provider {
  readonly mode = 'subscription' as const;
  hasCredential(): boolean { return !!readSubscription(this.cfg.credentialsPath)?.accessToken; }
  authRequired(): boolean { return !this.hasCredential(); }
  async forward(req: http.IncomingMessage, res: http.ServerResponse, pathname: string, search: string): Promise<void> {
    const creds = readSubscription(this.cfg.credentialsPath);
    if (!creds?.accessToken)
      return sendAuthError(res, 'no subscription token available — sign in required');
    const beta = mergeBeta(req.headers['anthropic-beta']);
    return proxy(req, res, this.cfg.upstreams.anthropic + pathname + search,
      {'Authorization': `Bearer ${creds.accessToken}`, 'anthropic-beta': beta},
      {onUnauthorized: async () => {
        const token = await refreshSubscription(this.cfg);
        return token ? {'Authorization': `Bearer ${token}`, 'anthropic-beta': beta} : null;
      }});
  }
}

class BedrockProvider extends Provider {
  readonly mode = 'bedrock' as const;
  hasCredential(): boolean {
    return !!(this.cfg.awsBearerToken || (this.cfg.awsAccessKeyId && this.cfg.awsSecretAccessKey));
  }
  async forward(req: http.IncomingMessage, res: http.ServerResponse, pathname: string, search: string): Promise<void> {
    if (!this.cfg.upstreams.bedrock)
      return send(res, 500, {error: 'bedrock upstream not configured'});
    const dest = this.cfg.upstreams.bedrock + pathname.replace(/^\/bedrock/, '') + search;
    if (this.cfg.awsBearerToken)
      return proxy(req, res, dest, {'Authorization': `Bearer ${this.cfg.awsBearerToken}`});
    if (this.cfg.awsAccessKeyId && this.cfg.awsSecretAccessKey && this.cfg.region) {
      return proxy(req, res, dest, {}, {sign: {
        accessKeyId: this.cfg.awsAccessKeyId, secretAccessKey: this.cfg.awsSecretAccessKey,
        sessionToken: this.cfg.awsSessionToken, region: this.cfg.region,
      }});
    }
    return send(res, 500, {error: 'bedrock selected but no credentials configured'});
  }
}

class FoundryProvider extends Provider {
  readonly mode = 'foundry' as const;
  hasCredential(): boolean { return !!this.cfg.foundryApiKey; }
  async forward(req: http.IncomingMessage, res: http.ServerResponse, pathname: string, search: string): Promise<void> {
    if (!this.cfg.upstreams.foundry || !this.cfg.foundryApiKey)
      return send(res, 500, {error: 'foundry upstream not configured'});
    return proxy(req, res, this.cfg.upstreams.foundry + pathname.replace(/^\/foundry/, '') + search,
      {'api-key': this.cfg.foundryApiKey});
  }
}

function makeProvider(cfg: BrokerConfig): Provider {
  switch (cfg.mode) {
    case 'bedrock': return new BedrockProvider(cfg);
    case 'foundry': return new FoundryProvider(cfg);
    case 'anthropic': return new AnthropicProvider(cfg);
    default: return new SubscriptionProvider(cfg);
  }
}

function createAuthRelay() {
  let proc: ChildProcess | null = null;
  let timer: NodeJS.Timeout | null = null;

  const kill = () => {
    proc?.removeAllListeners('exit');
    proc?.kill();
    if (timer)
      clearTimeout(timer);
    proc = null;
    timer = null;
  };

  const start = (): Promise<{url?: string; error?: string}> => new Promise((resolve) => {
    kill();
    proc = spawn('claude', ['auth', 'login'], {env: {...process.env, TERM: 'dumb', FORCE_COLOR: '0'}});
    timer = setTimeout(() => {
      kill();
      resolve({error: 'authentication timed out'});
    }, AUTH_LOGIN_TIMEOUT_MS);
    let done = false;
    const onData = (buf: Buffer) => {
      if (done)
        return;
      const m = buf.toString().match(/https:\/\/claude\.com\/cai\/oauth[^\s]+/);
      if (m) {
        done = true;
        resolve({url: m[0]});
      }
    };
    proc.stdout?.on('data', onData);
    proc.stderr?.on('data', onData);
    proc.on('error', (e) => {
      kill();
      resolve({error: e.message});
    });
  });

  const code = (value: string): Promise<{status: string; message?: string}> => new Promise((resolve) => {
    if (!proc?.stdin)
      return resolve({status: 'error', message: 'no authentication in progress'});
    proc.on('exit', (exitCode) => {
      kill();
      resolve(exitCode === 0 ? {status: 'done'} : {status: 'error', message: `login exited with code ${exitCode}`});
    });
    proc.stdin.write(value + '\n');
  });

  return {start, code};
}

export function createBroker(cfg: BrokerConfig): http.Server {
  const provider = makeProvider(cfg);
  const mcpTargets = new Map<string, McpTarget>();
  const auth = createAuthRelay();

  return http.createServer(async (req, res) => {
    try {
      const url = new URL(req.url || '/', `http://${BROKER_HOST}`);
      const p = url.pathname;

      if (p === '/health')
        return send(res, 200, {status: 'ok'});

      if (p === '/status') {
        return send(res, 200, {mode: provider.mode, authRequired: provider.authRequired(),
          region: cfg.region, foundryResource: cfg.foundryResource, models: cfg.models});
      }

      if (p === '/mcp-session' && req.method === 'POST') {
        const body = JSON.parse((await readBody(req)).toString() || '{}');
        if (!body.targetUrl)
          return send(res, 400, {error: 'targetUrl is required'});
        if (mcpTargets.size >= MAX_MCP_TARGETS)
          mcpTargets.delete(mcpTargets.keys().next().value as string);
        const token = randomUUID();
        mcpTargets.set(token, {targetUrl: body.targetUrl, apiKey: body.apiKey, apiUrl: body.apiUrl});
        return send(res, 200, {token});
      }

      if (p === '/auth/start' && req.method === 'POST')
        return send(res, 200, await auth.start());

      if (p === '/auth/code' && req.method === 'POST') {
        const body = JSON.parse((await readBody(req)).toString() || '{}');
        return send(res, 200, await auth.code(body.code ?? ''));
      }

      if (p.startsWith('/mcp/')) {
        const t = mcpTargets.get(p.split('/')[2]);
        if (!t)
          return send(res, 404, {error: 'unknown mcp session'});
        const inject: Record<string, string> = {};
        if (t.apiKey) {
          inject['Authorization'] = t.apiKey;
          inject['x-user-api-key'] = t.apiKey;
        }
        if (t.apiUrl)
          inject['x-datagrok-api-url'] = t.apiUrl;
        return proxy(req, res, t.targetUrl, inject);
      }

      return await provider.forward(req, res, p, url.search);
    } catch (e: any) {
      console.error('[broker] request error:', e.message);
      if (!res.headersSent)
        send(res, 500, {error: 'broker error'});
    }
  });
}

export function startBroker(): http.Server {
  const cfg = resolveConfig();
  const server = createBroker(cfg);
  server.listen(BROKER_PORT, BROKER_HOST, () =>
    console.log(`[broker] listening on ${BROKER_HOST}:${BROKER_PORT} (mode: ${cfg.mode})`));
  return server;
}

if (require.main === module)
  startBroker();
