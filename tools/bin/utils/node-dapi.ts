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

export class NodeApiClient {
  constructor(public baseUrl: string, public token: string) {}

  static async login(baseUrl: string, devKey: string): Promise<NodeApiClient> {
    const res = await fetch(`${baseUrl}/users/login/dev/${devKey}`, {method: 'POST'});
    const json = await res.json() as any;
    if (!json.token)
      throw new Error('Login failed. Check your developer key.');
    return new NodeApiClient(baseUrl, json.token);
  }

  async request(method: string, path: string, body?: any, headers?: Record<string, string>): Promise<any> {
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

    const res = await fetch(url, opts);

    if (!res.ok) {
      // Read as text first to avoid "Body has already been read" when JSON.parse fails
      const rawText = await res.text();
      let errBody: any;
      try { errBody = JSON.parse(rawText); } catch { errBody = {error: rawText || `HTTP ${res.status}`}; }
      const err: NodeApiError = {
        error: errBody?.message ?? errBody?.error ?? `HTTP ${res.status}`,
        source: errBody?.source ?? 'Server',
        errorCode: errBody?.errorCode ?? res.status,
        stackTrace: errBody?.stackTrace,
      };
      throw Object.assign(new Error(err.error), {apiError: err});
    }

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
}

function buildQuery(params: Record<string, any>): string {
  const entries = Object.entries(params).filter(([, v]) => v !== undefined && v !== null && v !== '');
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

export class NodeFilesDataSource {
  constructor(private client: NodeApiClient) {}

  private splitPath(filePath: string): {connector: string; path: string} {
    const colonIdx = filePath.indexOf(':');
    if (colonIdx === -1)
      return {connector: filePath.replace(':', '.'), path: ''};
    const connector = filePath.slice(0, colonIdx).replace(':', '.');
    const path = filePath.slice(colonIdx + 1).replace(/^\//, '');
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
    return this.client.get(`/public/v1/files/${connector}/${path}`);
  }

  async delete(filePath: string): Promise<void> {
    const {connector, path} = this.splitPath(filePath);
    await this.client.del(`/public/v1/files/${connector}/${path}`);
  }
}

export class NodeDapi {
  constructor(private client: NodeApiClient) {}

  get users(): NodeUsersDataSource { return new NodeUsersDataSource(this.client); }
  get groups(): NodeGroupsDataSource { return new NodeGroupsDataSource(this.client); }
  get functions(): NodeFuncsDataSource { return new NodeFuncsDataSource(this.client, 'functions'); }
  get connections(): NodeConnectionsDataSource { return new NodeConnectionsDataSource(this.client); }
  get queries(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'queries'); }
  get scripts(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'scripts'); }
  get packages(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'packages'); }
  get reports(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'reports'); }
  get files(): NodeFilesDataSource { return new NodeFilesDataSource(this.client); }
  get shares(): NodeSharesDataSource { return new NodeSharesDataSource(this.client); }

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
