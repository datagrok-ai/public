/// Docs: [Grok Dapi](/docs/plans/grok-dapi/)
import {randomUUID} from 'crypto';

export function ensureBodyId<T extends {id?: string} | Record<string, any>>(body: T): T {
  if (body && typeof body === 'object' && !(body as any).id)
    (body as any).id = randomUUID();
  return body;
}

export interface BatchOperation {
  id?: string;
  action: string;
  params: Record<string, any> | Array<Record<string, any>>;
  dependsOn?: Array<string | {id: string; allowPartial?: boolean}>;
}

export interface BatchRequest {
  operations: BatchOperation[];
  options?: {
    concurrency?: number;
    stopOnError?: boolean;
    transaction?: boolean;
  };
}

export interface BatchResponse {
  summary: {total: number; succeeded: number; partial: number; failed: number; skipped: number};
  results: Array<{
    id?: string;
    action: string;
    status: 'success' | 'error' | 'skipped' | 'partial';
    result?: any;
    error?: {error: string; errorCode?: number};
    reason?: string;
    summary?: {total: number; succeeded: number; failed: number};
    results?: Array<{index: number; status: string; result?: any; error?: any}>;
  }>;
}

export interface NodeApiError {
  error: string;
  source?: string;
  errorCode?: number;
  stackTrace?: string;
  innerError?: NodeApiError;
}

const setting = (name: string, fallback: number): number => {
  const value = Number(process.env[`GROK_HTTP_${name}`]);
  return Number.isFinite(value) && value >= 0 ? value : fallback;
};

const BYTES_TIMEOUT_MS = setting('BYTES_TIMEOUT', 600000);

/** Load shedding, not a verdict on the request: the same call succeeds once the queue drains. */
const RETRIABLE_STATUS = new Set([429, 502, 503, 504]);

/**
 * Without a deadline one unresponsive entity stalls a whole pull — `GET /projects/{id}` on a
 * space holding tens of thousands of children never answers. A request that hung or dropped is
 * retried, since a deadline is as often a server busy with this very pull as a dead one; a reply
 * the server actually sent is not. The deadline covers the body too, so a transfer that is slow
 * by nature rather than stuck (`.d42` table data) asks for a longer one.
 */
async function fetchOrRetry(url: string, opts: RequestInit, retriable: boolean,
                            timeoutMs: number = setting('TIMEOUT', 60000)): Promise<Response> {
  const retries = setting('RETRIES', 3);
  for (let attempt = 0; ; attempt++) {
    const last = !retriable || attempt >= retries;
    try {
      const res = await fetch(url, {...opts, signal: AbortSignal.timeout(timeoutMs)});
      if (last || !RETRIABLE_STATUS.has(res.status))
        return res;
      await res.body?.cancel();
    } catch (err: any) {
      if (last)
        throw new Error(`${opts.method ?? 'GET'} ${url}: ` +
          (err?.name === 'TimeoutError' ? `no answer in ${timeoutMs}ms` : err?.message ?? err));
    }
    await new Promise((resolve) => setTimeout(resolve, setting('BACKOFF', 1000) * Math.pow(2, attempt)));
  }
}

export class NodeApiClient {
  /** Set by `createClient` when the run asked for an admin session, so a re-login restores it. */
  adminMode: boolean = false;

  constructor(public baseUrl: string, public token: string, private devKey?: string) {}

  static async login(baseUrl: string, devKey: string): Promise<NodeApiClient> {
    const res = await fetch(`${baseUrl}/users/login/dev/${devKey}`, {method: 'POST'});
    const json = await res.json() as any;
    if (!json.token)
      throw new Error('Login failed. Check your developer key.');
    return new NodeApiClient(baseUrl, json.token, devKey);
  }

  /**
   * A stand serving several isolates can reject a session one of them does not know, and an
   * hour-long walk has no way to ask the operator to log in again. The developer key is good
   * for a new session, so one is taken rather than losing the run.
   */
  private async reauthenticate(): Promise<boolean> {
    if (!this.devKey)
      return false;
    const fresh = await NodeApiClient.login(this.baseUrl, this.devKey).catch(() => null);
    if (!fresh)
      return false;
    this.token = fresh.token;
    if (this.adminMode)
      this.token = (await fresh.post('/users/sessions/current/admin'))?.token ?? this.token;
    return true;
  }

  async request(method: string, path: string, body?: any, headers?: Record<string, string>,
                reauthed: boolean = false): Promise<any> {
    const url = `${this.baseUrl}${path}`;
    const opts: RequestInit = {
      method,
      headers: {
        'Authorization': this.token,
        'Content-Type': 'application/json',
        ...headers,
      },
    };
    if (body !== undefined)
      opts.body = JSON.stringify(body);

    const res = await fetchOrRetry(url, opts, method === 'GET');

    if (res.status === 401 && !reauthed && await this.reauthenticate())
      return this.request(method, path, body, headers, true);
    if (!res.ok)
      await throwHttpError(res);

    if (res.status === 204 || res.headers.get('content-length') === '0')
      return null;

    const ct = res.headers.get('content-type') ?? '';
    if (ct.includes('application/json'))
      return res.json();
    return res.text();
  }

  get(path: string): Promise<any> { return this.request('GET', path); }
  post(path: string, body?: any): Promise<any> { return this.request('POST', path, body); }
  del(path: string): Promise<any> { return this.request('DELETE', path); }

  /**
   * POST raw bytes — used for file/table uploads where the body must be the content
   * itself, not JSON. Defaults to `application/octet-stream`; pass `text/csv` (or
   * similar) when the server demands a specific content type.
   */
  async putBytes(path: string, bytes: Uint8Array | Buffer,
                 contentType: string = 'application/octet-stream'): Promise<any> {
    const res = await fetchOrRetry(`${this.baseUrl}${path}`, {
      method: 'POST',
      headers: {
        'Authorization': this.token,
        'Content-Type': contentType,
      },
      body: bytes as any,
    }, false, BYTES_TIMEOUT_MS);
    if (!res.ok)
      await throwHttpError(res);
    const ct = res.headers.get('content-type') ?? '';
    if (ct.includes('application/json'))
      return res.json();
    return res.text();
  }

  /** GET raw bytes — d42 table data, file content, model blobs. */
  async getBytes(path: string): Promise<Buffer> {
    const res = await fetchOrRetry(`${this.baseUrl}${path}`, {headers: {'Authorization': this.token}}, true, BYTES_TIMEOUT_MS);
    if (!res.ok)
      await throwHttpError(res);
    return Buffer.from(await res.arrayBuffer());
  }
}

// Read as text first to avoid "Body has already been read" when JSON.parse fails
async function throwHttpError(res: Response): Promise<never> {
  const rawText = await res.text();
  let errBody: any;
  // A gateway answers overload with an HTML page, where the status is the only real information.
  const markup = rawText.trimStart().startsWith('<') || rawText.length > 200;
  try { errBody = JSON.parse(rawText); }
  catch { errBody = {error: markup || !rawText ? `HTTP ${res.status} ${res.statusText}`.trim() : rawText}; }
  const err: NodeApiError = {
    error: errBody?.message ?? errBody?.error ?? `HTTP ${res.status}`,
    source: errBody?.source ?? 'Server',
    errorCode: errBody?.errorCode ?? res.status,
    stackTrace: errBody?.stackTrace,
  };
  throw Object.assign(new Error(err.error), {apiError: err});
}

/** Empty values are sent verbatim — an empty `namespace` selects the root namespace. */
export function buildQuery(params: Record<string, any>): string {
  const entries = Object.entries(params).filter(([, v]) => v !== undefined && v !== null);
  if (!entries.length) return '';
  return '?' + entries.map(([k, v]) => `${encodeURIComponent(k)}=${encodeURIComponent(String(v))}`).join('&');
}

export class NodeHttpDataSource<T = any> {
  protected _filter: string = '';
  protected _limit: number = 50;
  protected _page: number = 0;
  protected _order: string = '';

  constructor(protected client: NodeApiClient, protected path: string) {}

  filter(w: string): this { this._filter = w; return this; }
  by(n: number): this { this._limit = n; return this; }
  page(n: number): this { this._page = n; return this; }
  order(field: string, desc: boolean = false): this { this._order = desc ? `-${field}` : field; return this; }

  async list(): Promise<T[]> {
    const q = buildQuery({
      text: this._filter || undefined,
      limit: this._limit,
      page: this._page || undefined,
      order: this._order || undefined,
    });
    return this.client.get(`/public/v1/${this.path}${q}`);
  }

  async find(id: string): Promise<T> {
    return this.client.get(`/public/v1/${this.path}/${encodeURIComponent(id.replace(':', '.'))}`);
  }

  async count(): Promise<number> {
    const q = buildQuery({text: this._filter || undefined});
    const res = await this.client.get(`/public/v1/${this.path}/count${q}`);
    return typeof res === 'number' ? res : (res?.count ?? 0);
  }

  async delete(idOrEntity: string | {id?: string}): Promise<void> {
    const id = typeof idOrEntity === 'string' ? idOrEntity : (idOrEntity?.id ?? '');
    await this.client.del(`/public/v1/${this.path}/${encodeURIComponent(id)}`);
  }
}

/**
 * Generic client for the internal entity routers (`/projects`, `/scripts`,
 * `/connectors/queries`, ...) the browser itself uses. Unlike `NodeDapi.raw` it goes
 * through `client.request`, so a non-2xx response throws instead of returning the
 * error body as data.
 */
export class InternalDataSource {
  constructor(private client: NodeApiClient, public route: string) {}

  /**
   * The internal routers answer a missing or rejected entity with HTTP 200 and an
   * `ApiError` body, so success has to be decided from the payload, not the status.
   * Only a router that says so is a 404 — everything else is a server-side failure and
   * must not be mistaken for an absent entity.
   */
  private async call(method: string, path: string, body?: any): Promise<any> {
    const res = await this.client.request(method, path, body);
    const parsed = typeof res === 'string' ? tryParseJson(res) : res;
    if (parsed?.['#type'] === 'ApiError') {
      const err: NodeApiError = {error: parsed.message ?? 'Request failed', source: 'Server',
        errorCode: parsed.errorCode ?? 500, stackTrace: parsed.stackTrace};
      throw Object.assign(new Error(err.error), {apiError: err});
    }
    return parsed ?? res;
  }

  list(params: Record<string, any> = {}): Promise<any[]> {
    return this.call('GET', `${this.route}${buildQuery(params)}`);
  }

  /** Server paging is 1-based (`repository_query.dart` `paging`), so page 0 would repeat page 1. */
  async listAll(params: Record<string, any> = {}, pageSize: number = 500): Promise<any[]> {
    const all: any[] = [];
    for (let page = 1; ; page++) {
      const batch: any[] = await this.list({...params, limit: pageSize, page}) ?? [];
      all.push(...batch);
      if (batch.length < pageSize)
        return all;
    }
  }

  async find(id: string, include?: string): Promise<any> {
    try {
      return await this.call('GET', `${this.route}/${encodeURIComponent(id)}${buildQuery({include})}`);
    } catch (err: any) {
      if (err?.apiError?.errorCode === 404 || err?.message === 'Not Found')
        return null;
      throw err;
    }
  }

  save(json: any, query?: string): Promise<any> {
    return this.call('POST', `${this.route}${query ? '?' + query : ''}`, ensureBodyId(json));
  }

  delete(id: string): Promise<any> {
    return this.call('DELETE', `${this.route}/${encodeURIComponent(id)}`);
  }
}

export type MemberAddStatus = 'added' | 'updated' | 'noop' | 'error';
export type MemberRemoveStatus = 'removed' | 'not-member' | 'error';

export interface MemberAddResult {
  member: string;
  status: MemberAddStatus;
  error?: string;
}

export interface MemberRemoveResult {
  member: string;
  status: MemberRemoveStatus;
  error?: string;
}

const UUID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

export class NodeGroupsDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'groups'); }

  async save(group: any, saveRelations: boolean = false): Promise<any> {
    const q = buildQuery({saveRelations: saveRelations ? 'true' : undefined});
    return this.client.post(`/public/v1/groups${q}`, ensureBodyId(group));
  }

  async lookup(name: string): Promise<any[]> {
    const q = buildQuery({query: name});
    return this.client.get(`/public/v1/groups/lookup${q}`);
  }

  async resolve(idOrName: string, opts: {personalOnly?: boolean} = {}): Promise<any> {
    if (UUID_RE.test(idOrName))
      return this.find(idOrName);
    const matches: any[] = await this.lookup(idOrName);
    let candidates = matches;
    if (opts.personalOnly)
      candidates = matches.filter((g) => g?.personal === true);
    if (!candidates.length) {
      const suffix = opts.personalOnly ? ' (personal)' : '';
      throw new Error(`No group matching '${idOrName}'${suffix}`);
    }
    if (candidates.length > 1) {
      const list = candidates.map((g) => `  ${g.id}  ${g.friendlyName ?? g.name ?? ''}`).join('\n');
      throw new Error(`Multiple groups match '${idOrName}':\n${list}\nUse the ID to disambiguate.`);
    }
    return candidates[0];
  }

  async addMembers(group: string, members: string[], isAdmin: boolean = false, personalOnly: boolean = false): Promise<MemberAddResult[]> {
    // Always fetch via find() so parent.children comes back expanded; lookup() returns a
    // pruned projection and replacing that empty list on save would drop existing members.
    const resolved = await this.resolve(group);
    const parent = await this.find(resolved.id);
    const children: any[] = Array.isArray(parent.children) ? parent.children : [];
    const results: MemberAddResult[] = [];
    let mutated = false;

    for (const m of members) {
      let child: any;
      try {
        child = await this.resolve(m, {personalOnly});
      } catch (err: any) {
        results.push({member: m, status: 'error', error: err?.message ?? String(err)});
        continue;
      }
      const existing = children.find((r) => r?.child?.id === child.id);
      if (existing) {
        // Server returns isAdmin as null/undefined for non-admin relations; normalize
        // the comparison so re-runs report `noop` instead of `updated`.
        if ((existing.isAdmin ?? false) === isAdmin) {
          results.push({member: m, status: 'noop'});
        } else {
          existing.isAdmin = isAdmin;
          mutated = true;
          results.push({member: m, status: 'updated'});
        }
      } else {
        // Each GroupRelation row needs a non-null id; the server rejects the save otherwise.
        children.push({id: randomUUID(), parent: {id: parent.id}, child: {id: child.id}, isAdmin});
        mutated = true;
        results.push({member: m, status: 'added'});
      }
    }

    if (mutated) {
      parent.children = children;
      await this.save(parent, true);
    }
    return results;
  }

  async removeMembers(group: string, members: string[], personalOnly: boolean = false): Promise<MemberRemoveResult[]> {
    const resolved = await this.resolve(group);
    const parent = await this.find(resolved.id);
    const results: MemberRemoveResult[] = [];
    const children: any[] = Array.isArray(parent.children) ? parent.children : [];
    let mutated = false;
    for (const m of members) {
      let child: any;
      try {
        child = await this.resolve(m, {personalOnly});
      } catch (err: any) {
        results.push({member: m, status: 'error', error: err?.message ?? String(err)});
        continue;
      }
      const idx = children.findIndex((r) => r?.child?.id === child.id);
      if (idx === -1) {
        results.push({member: m, status: 'not-member'});
      } else {
        children.splice(idx, 1);
        mutated = true;
        results.push({member: m, status: 'removed'});
      }
    }

    if (mutated) {
      parent.children = children;
      await this.save(parent, true);
    }
    return results;
  }

  async getMembers(group: string, admin?: boolean): Promise<any[]> {
    const parent = await this.resolve(group);
    const q = buildQuery({admin: admin === undefined ? undefined : String(admin)});
    return this.client.get(`/public/v1/groups/${encodeURIComponent(parent.id)}/members${q}`);
  }

  async getMemberships(group: string, admin?: boolean): Promise<any[]> {
    const parent = await this.resolve(group);
    const q = buildQuery({admin: admin === undefined ? undefined : String(admin)});
    return this.client.get(`/public/v1/groups/${encodeURIComponent(parent.id)}/memberships${q}`);
  }
}

export class NodeSharesDataSource {
  constructor(private client: NodeApiClient) {}

  async share(entity: string, groups: string, access: string = 'View'): Promise<any> {
    const name = encodeURIComponent(entity.replace(':', '.'));
    const q = buildQuery({groups, access});
    return this.client.post(`/public/v1/entities/${name}/shares${q}`);
  }

  async list(entityId: string): Promise<any[]> {
    const q = buildQuery({entityId});
    return this.client.get(`/privileges/permissions${q}`);
  }
}

export class NodeUsersDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'users'); }

  async save(user: any): Promise<any> {
    return this.client.post('/public/v1/users', ensureBodyId(user));
  }

  async block(user: any): Promise<void> {
    await this.client.post('/public/v1/users/block', user);
  }

  async unblock(user: any): Promise<void> {
    await this.client.post('/public/v1/users/unblock', user);
  }
}

export class NodeConnectionsDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'connections'); }

  async save(conn: any, saveCredentials: boolean = false): Promise<any> {
    const q = buildQuery({saveCredentials: saveCredentials ? 'true' : undefined});
    return this.client.post(`/public/v1/connections${q}`, conn);
  }

  async test(conn: any): Promise<void> {
    const result = await this.client.post(`/public/v1/connections/test`, conn);
    const text = typeof result === 'string' ? result.replace(/^"|"$/g, '') : String(result ?? '');
    if (text !== 'ok')
      throw new Error(text || 'Connection test failed');
  }
}

export class NodeFuncsDataSource extends NodeHttpDataSource {
  async run(name: string, params?: Record<string, any>): Promise<any> {
    const normalizedName = name.replace(':', '.');
    const result = await this.client.post(`/public/v1/functions/${encodeURIComponent(normalizedName)}/call`, params ?? {});
    // Datagrok returns HTTP 200 with an ApiError body when the function doesn't exist
    const parsed = typeof result === 'string' ? tryParseJson(result) : result;
    if (parsed?.['#type'] === 'ApiError') {
      const err: NodeApiError = {error: parsed.message ?? 'Function call failed', errorCode: parsed.errorCode, stackTrace: parsed.stackTrace};
      throw Object.assign(new Error(err.error), {apiError: err});
    }
    return result;
  }
}

function tryParseJson(s: string): any {
  try { return JSON.parse(s); } catch { return null; }
}

export type PackageOpStatus = 'installed' | 'noop' | 'error';

export interface PackageOpResult {
  package: string;
  version: string;
  status: PackageOpStatus;
  id?: string;
  note?: string;
  error?: string;
}

export class NodePackagesDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'packages'); }

  /** GET /packages (internal API) — includes publishedVersions and repository. */
  async listFull(text?: string): Promise<any[]> {
    return this.client.get(`/packages${buildQuery({text})}`);
  }

  /** Resolve by UUID, name, or friendlyName (case-insensitive). Returns null when
      not found so callers can pass the raw string through and let the server's own
      'Package not found' surface. */
  async resolve(idOrName: string): Promise<any> {
    const all: any[] = await this.listFull();
    if (UUID_RE.test(idOrName))
      return all.find((p) => p?.id === idOrName) ?? null;
    const want = idOrName.toLowerCase();
    const matches = all.filter((p) =>
      (p?.name ?? '').toLowerCase() === want || (p?.friendlyName ?? '').toLowerCase() === want);
    if (matches.length > 1) {
      const list = matches.map((p) => `  ${p.id}  ${p.name}`).join('\n');
      throw new Error(`Multiple packages match '${idOrName}':\n${list}\nUse the ID to disambiguate.`);
    }
    return matches[0] ?? null;
  }

  /** Install/activate via the DeployPackageVersion server func: pulls the version
      from the package repository (npm) when needed, then makes it current.
      'latest' marks the package for server-side auto-update. Synchronous — returns
      the new published-package id, or null when nothing changed. */
  async install(name: string, desiredVersion: string = 'latest'): Promise<string | null> {
    const result = await new NodeFuncsDataSource(this.client, 'functions')
      .run('DeployPackageVersion', {name, desiredVersion});
    const value = (result && typeof result === 'object') ? (result.result ?? result.id ?? null) : result;
    const text = typeof value === 'string' ? value.replace(/^"|"$/g, '') : '';
    return text && text !== 'null' ? text : null;
  }

  async uninstall(idOrName: string): Promise<{id: string; name: string; repoBacked: boolean}> {
    const pkg = await this.resolve(idOrName);
    if (!pkg)
      throw new Error(`Package '${idOrName}' not found`);
    await this.client.del(`/packages/${encodeURIComponent(pkg.id)}`);
    return {id: pkg.id, name: pkg.name, repoBacked: !!pkg.repository};
  }

  async versions(idOrName: string): Promise<any[]> {
    const pkg = await this.resolve(idOrName);
    if (!pkg)
      throw new Error(`Package '${idOrName}' not found`);
    return Array.isArray(pkg.publishedVersions) ? pkg.publishedVersions : [];
  }

  async outdated(): Promise<Array<{name: string; installed: string; latest: string; desiredVersion: string}>> {
    const all: any[] = await this.listFull();
    const rows: Array<{name: string; installed: string; latest: string; desiredVersion: string}> = [];
    for (const p of all) {
      const versions: any[] = Array.isArray(p?.publishedVersions) ? p.publishedVersions : [];
      // isLocal is not reliably serialized, so 'installed' = has a current non-debug version
      const current = versions.find((v) => v?.isCurrent && !v?.debug);
      const latest = versions.find((v) => v?.isLatest);
      // Numeric compare so locally published builds (1.2.3.X-suffix) don't count
      // as outdated against an equal or older registry version.
      if (current && latest && NodePackagesDataSource.compareVersions(latest.version, current.version) > 0)
        rows.push({name: p.name, installed: current.version, latest: latest.version, desiredVersion: p.desiredVersion ?? ''});
    }
    return rows;
  }

  /** Compare the leading numeric dot-groups of two version strings (suffixes like
      '.X-4e32bb91' are ignored): positive when a > b. */
  static compareVersions(a: string, b: string): number {
    const nums = (s: string) => (String(s).match(/^\d+(\.\d+)*/)?.[0] ?? '').split('.').map(Number);
    const av = nums(a), bv = nums(b);
    for (let i = 0; i < Math.max(av.length, bv.length); i++) {
      const d = (av[i] ?? 0) - (bv[i] ?? 0);
      if (d !== 0) return d;
    }
    return 0;
  }

  /** The public API has no DELETE for packages — route to the internal endpoint. */
  async delete(idOrEntity: string | {id?: string}): Promise<void> {
    const id = typeof idOrEntity === 'string' ? idOrEntity : (idOrEntity?.id ?? '');
    await this.client.del(`/packages/${encodeURIComponent(id)}`);
  }
}

export class NodeFilesDataSource {
  constructor(private client: NodeApiClient) {}

  /**
   * Split a user-facing file path into {connector, path}.
   *
   * Input format: `<connector>/<file-path>` where `<connector>` is the connection's
   * full name — including namespace — e.g. `System:DemoFiles/smiles_1M.csv`.
   * The connector can contain colons (the namespace separator); the file path
   * starts after the first `/`. Colons in the connector segment are converted
   * to `.` so it forms a single URL path segment (the Dart server reverses
   * this with `replaceAll('.', ':')`).
   */
  private splitPath(filePath: string): {connector: string; path: string} {
    const slashIdx = filePath.indexOf('/');
    if (slashIdx === -1)
      return {connector: filePath.replace(/:/g, '.'), path: ''};
    const connector = filePath.slice(0, slashIdx).replace(/:/g, '.');
    const path = filePath.slice(slashIdx + 1);
    return {connector, path};
  }

  async list(filePath: string, recursive: boolean = false): Promise<any[]> {
    const {connector, path} = this.splitPath(filePath);
    const q = buildQuery({recursive: recursive ? 'true' : undefined});
    const seg = path ? `${connector}/${path}` : connector;
    return this.client.get(`/public/v1/files/${seg}${q}`);
  }

  async get(filePath: string): Promise<any> {
    const {connector, path} = this.splitPath(filePath);
    const seg = path ? `${connector}/${path}` : connector;
    return this.client.get(`/public/v1/files/${seg}`);
  }

  async delete(filePath: string): Promise<void> {
    const {connector, path} = this.splitPath(filePath);
    const seg = path ? `${connector}/${path}` : connector;
    await this.client.del(`/public/v1/files/${seg}`);
  }

  /**
   * Upload a local file to a Datagrok file share.
   * Streams raw bytes to POST `/public/v1/files/<connector>/<path>` — no base64,
   * no JSON wrapping — so it handles large files without blowing up memory.
   */
  async put(localPath: string, remotePath: string): Promise<any> {
    const fs = require('fs') as typeof import('fs');
    const {connector, path} = this.splitPath(remotePath);
    if (!path) throw new Error(`Remote path must include a file name after the connector: got '${remotePath}'`);
    const bytes = fs.readFileSync(localPath);
    const res = await this.client.putBytes(`/public/v1/files/${connector}/${path}`, bytes);
    return {path: remotePath, size: bytes.length, response: res};
  }
}

export class NodeTablesDataSource {
  constructor(private client: NodeApiClient) {}

  /** GET /public/v1/tables/{name} — returns CSV text. Name accepts UUID or `namespace:name`. */
  async download(name: string): Promise<string> {
    const seg = encodeURIComponent(name.replace(/:/g, '.'));
    const result = await this.client.get(`/public/v1/tables/${seg}`);
    // Datagrok returns HTTP 200 + ApiError JSON when the table isn't found.
    const parsed = typeof result === 'string' ? tryParseJson(result) : result;
    if (parsed && typeof parsed === 'object' && parsed['#type'] === 'ApiError') {
      const err: NodeApiError = {error: parsed.message ?? 'Table download failed', errorCode: parsed.errorCode, stackTrace: parsed.stackTrace};
      throw Object.assign(new Error(err.error), {apiError: err});
    }
    return result as string;
  }

  /**
   * POST /public/v1/tables/{name} with raw bytes. Returns `{ID, Grok name, Markup, URL}`.
   * Defaults to `text/csv`; pass `application/octet-stream` to upload a `.d42`
   * binary blob — the server content-negotiates on the header and persists either form.
   */
  async upload(name: string, localPath: string, contentType: string = 'text/csv'): Promise<any> {
    const fs = require('fs') as typeof import('fs');
    const bytes = fs.readFileSync(localPath);
    const seg = encodeURIComponent(name.replace(/:/g, '.'));
    return this.client.putBytes(`/public/v1/tables/${seg}`, bytes, contentType);
  }
}

export class NodeDapi {
  constructor(public client: NodeApiClient) {}

  get users(): NodeUsersDataSource { return new NodeUsersDataSource(this.client); }
  get groups(): NodeGroupsDataSource { return new NodeGroupsDataSource(this.client); }
  get functions(): NodeFuncsDataSource { return new NodeFuncsDataSource(this.client, 'functions'); }
  get connections(): NodeConnectionsDataSource { return new NodeConnectionsDataSource(this.client); }
  get queries(): InternalDataSource { return this.internal('/connectors/queries'); }
  get scripts(): InternalDataSource { return this.internal('/scripts'); }
  get packages(): NodePackagesDataSource { return new NodePackagesDataSource(this.client); }
  get reports(): InternalDataSource { return this.internal('/reports'); }
  get files(): NodeFilesDataSource { return new NodeFilesDataSource(this.client); }
  get shares(): NodeSharesDataSource { return new NodeSharesDataSource(this.client); }
  get tables(): NodeTablesDataSource { return new NodeTablesDataSource(this.client); }

  internal(route: string): InternalDataSource { return new InternalDataSource(this.client, route); }

  async serverInfo(): Promise<{version: string; commit?: string}> {
    const raw = await this.client.get('/info/server');
    const info = typeof raw === 'string' ? JSON.parse(raw) : raw;
    return {version: info?.Version ?? '', commit: info?.Commit};
  }

  async raw(method: string, path: string, body?: any): Promise<any> {
    // Raw paths are relative to server root (e.g. /api/users/current).
    // Strip the trailing /api from baseUrl to avoid double prefix.
    const serverRoot = this.client.baseUrl.replace(/\/api\/?$/, '');
    const url = `${serverRoot}${path}`;
    const opts: RequestInit = {
      method: method.toUpperCase(),
      headers: {
        'Authorization': this.client.token,
        'Content-Type': 'application/json',
      } as Record<string, string>,
    };
    if (body !== undefined) opts.body = JSON.stringify(body);
    const res = await fetch(url, opts);
    const ct = res.headers.get('content-type') ?? '';
    if (ct.includes('application/json')) return res.json();
    return res.text();
  }

  async batch(request: BatchRequest): Promise<BatchResponse> {
    return this.client.post('/public/v1/batch', request);
  }

  async describe(entityType: string): Promise<any> {
    try {
      return await this.client.get(`/public/v1/entity-types/${encodeURIComponent(entityType)}`);
    } catch {
      return await this.client.get(`/entities/types${buildQuery({name: entityType})}`);
    }
  }
}
