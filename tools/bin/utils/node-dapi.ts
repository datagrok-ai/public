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
  /** The decoded error envelope, for callers that read structured fields (per-row reports, plans). */
  body?: any;
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
    const json = (res.headers.get('content-type') ?? '').includes('application/json') ? await res.json() as any : null;
    if (!json)
      throw new Error(`Login failed at ${baseUrl} (HTTP ${res.status}): not a Datagrok API URL — it should end with /api`);
    if (!json.token)
      throw new Error(`Login failed at ${baseUrl}: ${json.message ?? 'check your developer key'}`);
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
                reauthed: boolean = false, timeoutMs?: number): Promise<any> {
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

    const res = await fetchOrRetry(url, opts, method === 'GET', timeoutMs);

    if (res.status === 401 && !reauthed && await this.reauthenticate())
      return this.request(method, path, body, headers, true, timeoutMs);
    if (!res.ok)
      await throwHttpError(res);

    if (res.status === 204 || res.headers.get('content-length') === '0')
      return null;

    const ct = res.headers.get('content-type') ?? '';
    return throwIfApiError(ct.includes('application/json') ? await res.json() : await res.text());
  }

  get(path: string): Promise<any> { return this.request('GET', path); }
  post(path: string, body?: any, timeoutMs?: number): Promise<any> {
    return this.request('POST', path, body, undefined, false, timeoutMs);
  }
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
    return throwIfApiError(ct.includes('application/json') ? await res.json() : await res.text());
  }

  /** POST a JSON body and read the response as raw bytes (d42 query results). */
  async postForBytes(path: string, body: any): Promise<Buffer> {
    const res = await fetch(`${this.baseUrl}${path}`, {
      method: 'POST',
      headers: {'Authorization': this.token, 'Content-Type': 'application/json'},
      body: JSON.stringify(body),
    });
    if (!res.ok)
      await throwHttpError(res);
    return Buffer.from(await res.arrayBuffer());
  }

  /** GET raw bytes — d42 table data, file content, model blobs. */
  async getBytes(path: string): Promise<Buffer> {
    const res = await fetchOrRetry(`${this.baseUrl}${path}`, {headers: {'Authorization': this.token}}, true, BYTES_TIMEOUT_MS);
    if (!res.ok)
      await throwHttpError(res);
    const bytes = Buffer.from(await res.arrayBuffer());
    // A table whose data file is missing answers 200 — as text/plain — with an ApiError body.
    // Writing that into the bundle ships an error message as a table and only fails on the far
    // stand at push time; d42 never starts with `{`, so an envelope here is an error, not data.
    if (bytes[0] === 0x7B)
      throwIfApiError(bytes.toString('utf8'));
    return bytes;
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
    body: errBody,
  };
  throw Object.assign(new Error(err.error), {apiError: err});
}

/**
 * The server answers many failures with HTTP 200 and an `ApiError` body — as an object, or as
 * a JSON string when the handler returned text. Every response passes through here so a
 * missing entity, a rejected save or an unknown function throws instead of printing as data.
 */
export function throwIfApiError<T>(payload: T): T {
  const parsed: any = typeof payload === 'string' ? tryParseJson(payload) : payload;
  if (parsed?.['#type'] === 'ApiError') {
    const err: NodeApiError = {error: parsed.message ?? 'Request failed', source: 'Server',
      errorCode: parsed.errorCode, stackTrace: parsed.stackTrace, body: parsed};
    throw Object.assign(new Error(err.error), {apiError: err});
  }
  return payload;
}

/**
 * `grok s raw` paths are API-relative; a leading `/api` is accepted and dropped, so the same
 * path works against an nginx-fronted `.../api` base and a bare Datlas root.
 */
export function apiPath(path: string): string {
  const p = path.startsWith('/') ? path : `/${path}`;
  return p === '/api' ? '/' : p.startsWith('/api/') ? p.slice(4) : p;
}

/** Empty values are sent verbatim — an empty `namespace` selects the root namespace. */
export function buildQuery(params: Record<string, any>): string {
  const entries = Object.entries(params).filter(([, v]) => v !== undefined && v !== null);
  if (!entries.length) return '';
  return '?' + entries.map(([k, v]) => `${encodeURIComponent(k)}=${encodeURIComponent(String(v))}`).join('&');
}

/**
 * A public-API entity. `find`/`save`/`delete` go to `/public/v1/<path>`; `list` and `count` go to
 * [listRoute] when given — the internal router that pages (`limit`, 1-based `page`, `order`) and
 * has `/count`, which the public list routes of users, groups, connections and functions do not.
 */
export class NodeHttpDataSource<T = any> {
  protected _filter: string = '';
  protected _limit: number = 50;
  protected _page: number = 0;
  protected _order: string = '';

  constructor(protected client: NodeApiClient, protected path: string, protected listRoute?: string) {}

  filter(w: string): this { this._filter = w; return this; }
  by(n: number): this { this._limit = n; return this; }
  /** Zero-based page of `by(n)` rows. */
  page(n: number): this { this._page = n; return this; }
  /** Smart-order syntax: `!field` is descending. */
  order(field: string, desc: boolean = false): this { this._order = desc ? `!${field}` : field; return this; }

  async list(): Promise<T[]> {
    const q = buildQuery({
      text: this._filter || undefined,
      limit: this._limit,
      page: this._page + 1,
      order: this._order || undefined,
    });
    return this.client.get(`${this.listRoute ?? `/public/v1/${this.path}`}${q}`);
  }

  async find(id: string): Promise<T> {
    return this.client.get(`/public/v1/${this.path}/${encodeURIComponent(id.replace(':', '.'))}`);
  }

  async count(): Promise<number> {
    const q = buildQuery({text: this._filter || undefined});
    const res = await this.client.get(`${this.listRoute ?? `/public/v1/${this.path}`}/count${q}`);
    return typeof res === 'number' ? res : Number(res?.count ?? res ?? 0);
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
    return (typeof res === 'string' ? tryParseJson(res) : res) ?? res;
  }

  list(params: Record<string, any> = {}): Promise<any[]> {
    return this.call('GET', `${this.route}${buildQuery(params)}`);
  }

  async count(params: Record<string, any> = {}): Promise<number> {
    return Number(await this.call('GET', `${this.route}/count${buildQuery(params)}`));
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
  constructor(client: NodeApiClient) { super(client, 'groups', '/groups'); }

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
    // lookup is a substring search ('admin' also finds 'Administrators'): an exact name wins
    const want = idOrName.toLowerCase();
    const exact = candidates.filter((g) => [g?.name, g?.friendlyName].some((n) => (n ?? '').toLowerCase() === want));
    if (exact.length) candidates = exact;
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

  async getMembers(group: string, admin?: boolean, personalOnly: boolean = false): Promise<any[]> {
    const parent = await this.resolve(group, {personalOnly});
    const q = buildQuery({admin: admin === undefined ? undefined : String(admin)});
    return this.client.get(`/public/v1/groups/${encodeURIComponent(parent.id)}/members${q}`);
  }

  async getMemberships(group: string, admin?: boolean, personalOnly: boolean = false): Promise<any[]> {
    const parent = await this.resolve(group, {personalOnly});
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

  /** `all` walks the entity's project links the way the browser's sharing dialog does; a bare `entityId` answers nothing for most entities. */
  async list(entityId: string): Promise<any[]> {
    const q = buildQuery({entityId, all: 'true'});
    return this.client.get(`/privileges/permissions${q}`);
  }
}

export class NodeUsersDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'users', '/users'); }

  async save(user: any): Promise<any> {
    return this.client.post('/public/v1/users', ensureBodyId(user));
  }

  /**
   * Datagrok has no user deletion: the server offers no route, the UI only blocks, and removing
   * the entity record (`DELETE /entities/{id}`, what the batch API does) leaves the `users` row,
   * the personal group and the root project behind, so the login can never be re-created.
   */
  async delete(idOrLogin: string | {id?: string}): Promise<void> {
    const login = typeof idOrLogin === 'string' ? idOrLogin : (idOrLogin?.id ?? '');
    throw new Error(`Users cannot be deleted through the API; block the account instead: grok s users block ${login}`);
  }

  async block(user: any): Promise<void> {
    await this.client.post('/public/v1/users/block', user);
  }

  async unblock(user: any): Promise<void> {
    await this.client.post('/public/v1/users/unblock', user);
  }
}

export class NodeConnectionsDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'connections', '/connectors/connections'); }

  async save(conn: any, saveCredentials: boolean = false): Promise<any> {
    const q = buildQuery({saveCredentials: saveCredentials ? 'true' : undefined});
    return this.client.post(`/public/v1/connections${q}`, conn);
  }

  /** The route answers 200 for an unknown id, so look the connection up first. */
  async delete(idOrName: string | {id?: string}): Promise<void> {
    const conn = typeof idOrName === 'string' ? await this.find(idOrName) : idOrName;
    await super.delete(conn?.id ?? '');
  }

  async test(conn: any): Promise<void> {
    const result = await this.client.post(`/public/v1/connections/test`, conn);
    const text = typeof result === 'string' ? result.replace(/^"|"$/g, '') : String(result ?? '');
    if (text !== 'ok')
      throw new Error(text || 'Connection test failed');
  }
}

export class NodeFuncsDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'functions', '/log/funcs'); }

  async run(name: string, params?: Record<string, any>): Promise<any> {
    const normalizedName = name.replace(':', '.');
    return this.client.post(`/public/v1/functions/${encodeURIComponent(normalizedName)}/call`, params ?? {});
  }

  /** The public API has no DELETE for functions; scripts and queries go to their own routers. */
  async delete(idOrName: string | {id?: string}): Promise<void> {
    const func: any = typeof idOrName === 'string' ? await this.find(idOrName) : idOrName;
    const route = func?.['#type'] === 'Script' ? '/scripts' : func?.['#type'] === 'DataQuery' ? '/connectors/queries' : null;
    if (!route)
      throw new Error(`Only scripts and queries can be deleted; '${func?.name ?? idOrName}' is a ${func?.['#type'] ?? 'function'} (package functions go away with their package)`);
    await this.client.del(`${route}/${encodeURIComponent(func.id)}`);
  }
}

/**
 * The server binds a call's arguments by parameter name, so positional values have to be
 * mapped onto the function's inputs (`parameterInfos`, declared order; outputs carry
 * `isInput: false`) before the call.
 */
export function mapPositionalParams(params: Record<string, any>, parameterInfos: Record<string, any> | undefined,
                                    funcName: string): Record<string, any> {
  const positional = Object.keys(params).filter((k) => /^\d+$/.test(k)).sort((a, b) => Number(a) - Number(b));
  if (!positional.length) return params;
  const inputs = Object.values(parameterInfos ?? {}).filter((p: any) => p?.isInput !== false).map((p: any) => p.name);
  if (positional.length > inputs.length)
    throw new Error(`${funcName} takes ${inputs.length} input${inputs.length === 1 ? '' : 's'} (${inputs.join(', ')}), got ${positional.length}`);
  const mapped: Record<string, any> = {};
  for (const [k, v] of Object.entries(params))
    mapped[/^\d+$/.test(k) ? inputs[Number(k)] : k] = v;
  return mapped;
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

  /** No `/count` route for packages; the catalog is small enough to count client-side. */
  async count(): Promise<number> {
    return (await this.listFull(this._filter || undefined)).length;
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
    const result = await new NodeFuncsDataSource(this.client)
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

  /**
   * The public files route only downloads, so listing goes through the internal
   * `/connectors/connections/{id}/files/<path>` the file browser uses (a trailing `/`
   * is what makes it a listing). Resolves to `FileInfo` records (`path`, `isFile`, `size`).
   */
  async list(filePath: string, recursive: boolean = false): Promise<any[]> {
    const {connector, path} = this.splitPath(filePath);
    const conn = await new NodeConnectionsDataSource(this.client).find(connector);
    const q = buildQuery({recursive: recursive ? 'true' : undefined});
    const seg = path ? `${path.replace(/\/+$/, '')}/` : '';
    const res = await this.client.get(`/connectors/connections/${encodeURIComponent(conn.id)}/files/${seg}${q}`);
    return Array.isArray(res) ? res : [];
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

/**
 * Tables live under a project namespace (`Admin:MyTable:MyTable`), so a table is addressed by
 * UUID, by that full name, or by its bare name when only one table carries it.
 */
export class NodeTablesDataSource extends NodeHttpDataSource {
  constructor(client: NodeApiClient) { super(client, 'tables', '/tables'); }

  async find(idOrName: string): Promise<any> {
    if (UUID_RE.test(idOrName))
      return this.client.get(`/tables/${encodeURIComponent(idOrName)}`);
    const all: any[] = await this.client.get(`/tables${buildQuery({text: idOrName, limit: 100})}`);
    const matches = all.filter((t) => t?.name === idOrName || `${t?.namespace ?? ''}${t?.name ?? ''}` === idOrName);
    if (!matches.length)
      throw new Error(`No table named '${idOrName}'`);
    if (matches.length > 1) {
      const list = matches.map((t) => `  ${t.id}  ${t.namespace ?? ''}${t.name}`).join('\n');
      throw new Error(`Multiple tables match '${idOrName}':\n${list}\nUse the full name or the ID.`);
    }
    return matches[0];
  }

  async delete(idOrName: string | {id?: string}): Promise<void> {
    const table = typeof idOrName === 'string' ? await this.find(idOrName) : idOrName;
    await this.client.del(`/tables/${encodeURIComponent(table?.id ?? '')}`);
  }

  /** GET /public/v1/tables/{id} — returns CSV text. */
  async download(idOrName: string): Promise<string> {
    const table = await this.find(idOrName);
    return this.client.get(`/public/v1/tables/${encodeURIComponent(table.id)}`) as Promise<string>;
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

export type DomainAddress = {schema: string; table?: string};

/** `'<schema>'` or `'<schema>.<table>'` — the address every `grok s domains` verb takes. */
export function parseDomainAddress(s: string, opts: {table?: boolean} = {}): DomainAddress {
  const str = String(s ?? '');
  const dot = str.indexOf('.');
  const address = dot === -1 ? {schema: str} : {schema: str.slice(0, dot), table: str.slice(dot + 1)};
  if (!address.schema || (dot !== -1 && !address.table))
    throw new Error(`Invalid domain address '${str}': expected <schema> or <schema>.<table>`);
  if (opts.table === true && !address.table)
    throw new Error(`'${str}' names a schema; this command needs a table: <schema>.<table>`);
  if (opts.table === false && address.table)
    throw new Error(`'${str}' names a table; this command needs a schema`);
  return address;
}

export interface DomainBatchOptions {
  mode?: 'insert' | 'upsert';
  allOrNothing?: boolean;
  errorOnDuplicate?: boolean;
}

/**
 * Client for entity-mapped domain tables (`/domains/...`, DomainsRouter). Rows are plain
 * JSON objects; the server validates, permission-checks and audits every write. Errors
 * arrive as HTTP 4xx with a JSON envelope, so `client.request` throws with `apiError.body`
 * carrying the structured fields (per-row `rows`, a dry-run `plan`, version numbers).
 */
export class NodeDomainsDataSource {
  constructor(private client: NodeApiClient) {}

  private rows(schema: string, table: string): string {
    return `/domains/${encodeURIComponent(schema)}/${encodeURIComponent(table)}`;
  }

  /** Registered schemas with their tables (`GET /domains/schemas`); [text] is a smart filter. */
  schemas(text?: string): Promise<any[]> {
    return this.client.get(`/domains/schemas${buildQuery({text: text || undefined})}`);
  }

  async schema(name: string): Promise<any> {
    const all: any[] = await this.schemas();
    const s = all.find((x) => x?.name === name);
    if (!s) throw new Error(`Domain schema '${name}' not found`);
    return s;
  }

  /** Registry entity id of a schema or a table — the target of the grants endpoints. */
  async entityId(address: DomainAddress): Promise<string> {
    const s = await this.schema(address.schema);
    if (!address.table) return s.id;
    const t = (s.tables ?? []).find((x: any) => x?.name === address.table);
    if (!t) throw new Error(`Domain table '${address.schema}.${address.table}' not found`);
    return t.id;
  }

  manifest(schema: string): Promise<any> {
    return this.client.get(`/domains/schemas/${encodeURIComponent(schema)}/manifest`);
  }

  createSchema(name: string, friendlyName?: string, description?: string): Promise<any> {
    return this.client.post('/domains/schemas', {name, friendlyName, description});
  }

  applySchema(schema: string, body: any, dryRun: boolean = false): Promise<any> {
    const q = buildQuery({dryRun: dryRun ? 'true' : undefined});
    return this.client.post(`/domains/schemas/${encodeURIComponent(schema)}/apply${q}`, body);
  }

  deleteSchema(schema: string): Promise<any> {
    return this.client.del(`/domains/schemas/${encodeURIComponent(schema)}`);
  }

  schemaAudit(schema: string, limit?: number): Promise<any[]> {
    return this.client.get(`/domains/schemas/${encodeURIComponent(schema)}/audit${buildQuery({limit})}`);
  }

  tableAudit(schema: string, table: string, limit?: number): Promise<any[]> {
    return this.client.get(`${this.rows(schema, table)}/audit${buildQuery({limit})}`);
  }

  rowAudit(schema: string, table: string, id: string): Promise<any[]> {
    return this.client.get(`${this.rows(schema, table)}/${encodeURIComponent(id)}/audit`);
  }

  grants(entityId: string): Promise<any[]> {
    return this.client.get(`/domains/grants/${encodeURIComponent(entityId)}`);
  }

  grant(entityId: string, group: string, permission: string): Promise<any> {
    return this.client.post(`/domains/grants/${encodeURIComponent(entityId)}`, {group, permission});
  }

  revoke(entityId: string, group: string, permission?: string): Promise<any> {
    return this.client.del(`/domains/grants/${encodeURIComponent(entityId)}${buildQuery({group, permission})}`);
  }

  capabilities(schema: string, table: string): Promise<any> {
    return this.client.get(`${this.rows(schema, table)}/capabilities`);
  }

  /** JSON rows; spec = {filter, sort, columns, expand, limit, offset} (10k row cap). */
  query(schema: string, table: string, spec: Record<string, any> = {}): Promise<any[]> {
    return this.client.post(`${this.rows(schema, table)}/query`, spec);
  }

  /** The same query as a d42 DataFrame blob (10M row cap). */
  queryD42(schema: string, table: string, spec: Record<string, any> = {}): Promise<Buffer> {
    return this.client.postForBytes(`${this.rows(schema, table)}/query`, {...spec, format: 'd42'});
  }

  aggregate(schema: string, table: string, spec: Record<string, any>): Promise<any[]> {
    return this.client.post(`${this.rows(schema, table)}/aggregate`, spec);
  }

  async count(schema: string, table: string, filter?: any): Promise<number> {
    const spec: Record<string, any> = {measures: [{fn: 'count'}]};
    if (filter) spec.filter = filter;
    const rows = await this.aggregate(schema, table, spec);
    return Number(rows?.[0]?.count ?? 0);
  }

  /** One row, or null when it does not exist or is not visible (the server answers 404). */
  async getRow(schema: string, table: string, id: string): Promise<any> {
    try {
      return await this.client.get(`${this.rows(schema, table)}/${encodeURIComponent(id)}`);
    } catch (err: any) {
      if (err?.apiError?.errorCode === 404) return null;
      throw err;
    }
  }

  /** [rows] is one row object or an array; resolves to per-row `{id, created}` reports. */
  async insert(schema: string, table: string, rows: any, errorOnDuplicate: boolean = false): Promise<any[]> {
    const q = buildQuery({errorOnDuplicate: errorOnDuplicate ? 'true' : undefined});
    const res = await this.client.post(`${this.rows(schema, table)}${q}`, rows);
    return Array.isArray(res) ? res : [res];
  }

  update(schema: string, table: string, id: string, values: any, version?: number): Promise<any> {
    const body: Record<string, any> = {values};
    if (version !== undefined) body.version = version;
    return this.client.request('PATCH', `${this.rows(schema, table)}/${encodeURIComponent(id)}`, body);
  }

  deleteRow(schema: string, table: string, id: string): Promise<any> {
    return this.client.del(`${this.rows(schema, table)}/${encodeURIComponent(id)}`);
  }

  /** Soft-deletes up to [limit] (≤1000) matching rows in one transaction; `{deleted, hasMore}`. */
  deleteWhere(schema: string, table: string, filter: any, limit?: number): Promise<{deleted: number; hasMore: boolean}> {
    const body: Record<string, any> = {filter};
    if (limit !== undefined) body.limit = limit;
    return this.client.post(`${this.rows(schema, table)}/delete`, body);
  }

  /**
   * Bulk upload: [bytes] is the file content and [contentType] selects how the server
   * reads it — `text/csv`, `application/octet-stream` (d42), or `application/json` (a row
   * array, bare or under `rows`). Resolves to the batch report `{inserted, updated,
   * skipped, errorCount, rows}`; a report-carrying failure rejects with the report in
   * `apiError.body`.
   */
  batch(schema: string, table: string, bytes: Buffer, contentType: string, options: DomainBatchOptions = {}): Promise<any> {
    const q = buildQuery({
      mode: options.mode ?? 'insert',
      allOrNothing: options.allOrNothing === false ? 'false' : 'true',
      errorOnDuplicate: options.errorOnDuplicate ? 'true' : undefined,
    });
    return this.client.putBytes(`${this.rows(schema, table)}/batch${q}`, bytes, contentType);
  }

  /** Ordered ops (`{op, table, ref?, values?, id?, expectedVersion?}`) applied atomically. */
  transaction(schema: string, ops: any[]): Promise<any[]> {
    return this.client.post(`/domains/${encodeURIComponent(schema)}/transaction`, ops);
  }
}

export class NodeDapi {
  constructor(public client: NodeApiClient) {}

  get users(): NodeUsersDataSource { return new NodeUsersDataSource(this.client); }
  get groups(): NodeGroupsDataSource { return new NodeGroupsDataSource(this.client); }
  get functions(): NodeFuncsDataSource { return new NodeFuncsDataSource(this.client); }
  get connections(): NodeConnectionsDataSource { return new NodeConnectionsDataSource(this.client); }
  get queries(): InternalDataSource { return this.internal('/connectors/queries'); }
  get scripts(): InternalDataSource { return this.internal('/scripts'); }
  get packages(): NodePackagesDataSource { return new NodePackagesDataSource(this.client); }
  get reports(): InternalDataSource { return this.internal('/reports'); }
  get files(): NodeFilesDataSource { return new NodeFilesDataSource(this.client); }
  get shares(): NodeSharesDataSource { return new NodeSharesDataSource(this.client); }
  get tables(): NodeTablesDataSource { return new NodeTablesDataSource(this.client); }
  get domains(): NodeDomainsDataSource { return new NodeDomainsDataSource(this.client); }

  internal(route: string): InternalDataSource { return new InternalDataSource(this.client, route); }

  async serverInfo(): Promise<{version: string; commit?: string}> {
    const raw = await this.client.get('/info/server');
    const info = typeof raw === 'string' ? JSON.parse(raw) : raw;
    return {version: info?.Version ?? '', commit: info?.Commit};
  }

  /** Any API endpoint; [path] is API-relative (`/users/current`), a leading `/api` is accepted. */
  async raw(method: string, path: string, body?: any): Promise<any> {
    return this.client.request(method.toUpperCase(), apiPath(path), body);
  }

  async batch(request: BatchRequest): Promise<BatchResponse> {
    return this.client.post('/public/v1/batch', request);
  }

  /**
   * The shape of an entity type: its registry record (`/entities/types`) plus the top-level
   * fields of one existing entity of that type, since the server publishes no JSON schema.
   * [nameOrAlias] is a `grok s` entity (`connections`) or a type name (`DataConnection`).
   */
  async describe(nameOrAlias: string): Promise<{type: any; fields: DescribedField[]; sample: any}> {
    const alias = DESCRIBE_ALIASES[nameOrAlias.toLowerCase()];
    const typeName = alias?.type ?? nameOrAlias;
    const types: any[] = await this.client.get(`/entities/types${buildQuery({showSystem: 'true', text: `name="${typeName}"`, limit: 5})}`);
    const type = types.find((t) => (t?.name ?? '').toLowerCase() === typeName.toLowerCase()) ?? null;
    const sampleFrom = alias?.sample ?? DESCRIBE_SAMPLES[typeName];
    if (!type && !sampleFrom)
      throw new Error(`Unknown entity type '${nameOrAlias}'. Try one of: ${Object.keys(DESCRIBE_ALIASES).join(', ')}, or a type name such as Project`);
    const sample = sampleFrom ? (await sampleFrom(this))[0] ?? null : null;
    const fields = sample ? Object.entries(sample).map(([field, v]) => ({field, type: describeType(v), example: describeExample(v)})) : [];
    return {type, fields, sample};
  }
}

export interface DescribedField { field: string; type: string; example: string }

type SampleFetch = (dapi: NodeDapi) => Promise<any[]>;

const DESCRIBE_SAMPLES: Record<string, SampleFetch> = {
  User: (d) => d.users.by(1).list(),
  UserGroup: (d) => d.groups.by(1).list(),
  DataConnection: (d) => d.connections.by(1).list(),
  DataQuery: (d) => d.queries.list({limit: 1}),
  Script: (d) => d.scripts.list({limit: 1}),
  Func: (d) => d.functions.by(1).list(),
  Package: (d) => d.packages.by(1).list(),
  UserReport: (d) => d.reports.list({limit: 1}),
  TableInfo: (d) => d.tables.by(1).list(),
  Project: (d) => d.internal('/projects').list({limit: 1}),
  FileInfo: (d) => d.internal('/files').list({limit: 1}),
};

const DESCRIBE_ALIASES: Record<string, {type: string; sample?: SampleFetch}> = {
  users: {type: 'User'}, groups: {type: 'UserGroup'}, connections: {type: 'DataConnection'},
  queries: {type: 'DataQuery'}, scripts: {type: 'Script'}, functions: {type: 'Func'},
  packages: {type: 'Package'}, reports: {type: 'UserReport'}, tables: {type: 'TableInfo'},
  projects: {type: 'Project'}, files: {type: 'FileInfo'},
};

function describeType(v: any): string {
  if (v === null || v === undefined) return 'null';
  if (Array.isArray(v)) return `array(${v.length})`;
  if (typeof v === 'object') return v['#type'] ?? (v.id ? 'ref' : 'object');
  return typeof v;
}

function describeExample(v: any): string {
  const s = v === null || v === undefined ? '' : typeof v === 'object' ? JSON.stringify(v) : String(v);
  return s.length > 60 ? `${s.slice(0, 57)}...` : s;
}
