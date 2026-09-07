/// Docs: [Grok Dapi](/docs/plans/grok-dapi/)
/// Intentional subset mirror of `public/tools/bin/utils/node-dapi.ts` — the CLI keeps the
/// full client (entity data sources, batch, migrate); keep the shared parts in sync.
/// `raw` diverges on purpose: here it is a thin `client.request` (paths relative to the
/// client's base url), while the CLI's builds the url from the server root for `grok s raw`.
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
    const res = await fetch(`${this.baseUrl}${path}`, {
      method: 'POST',
      headers: {
        'Authorization': this.token,
        'Content-Type': contentType,
      },
      body: bytes as any,
    });
    if (!res.ok)
      await throwHttpError(res);
    const ct = res.headers.get('content-type') ?? '';
    if (ct.includes('application/json'))
      return res.json();
    return res.text();
  }

  /** GET raw bytes — d42 table data, file content, model blobs. */
  async getBytes(path: string): Promise<Buffer> {
    const res = await fetch(`${this.baseUrl}${path}`, {headers: {'Authorization': this.token}});
    if (!res.ok)
      await throwHttpError(res);
    return Buffer.from(await res.arrayBuffer());
  }
}

// Read as text first to avoid "Body has already been read" when JSON.parse fails
async function throwHttpError(res: Response): Promise<never> {
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

function tryParseJson(s: string): any {
  try { return JSON.parse(s); } catch { return null; }
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
}

export class NodeDapi {
  constructor(public client: NodeApiClient) {}

  get users(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'users'); }
  get groups(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'groups'); }
  get functions(): NodeFuncsDataSource { return new NodeFuncsDataSource(this.client, 'functions'); }
  get connections(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'connections'); }
  get queries(): InternalDataSource { return this.internal('/connectors/queries'); }
  get scripts(): InternalDataSource { return this.internal('/scripts'); }
  get packages(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'packages'); }
  get reports(): InternalDataSource { return this.internal('/reports'); }
  get files(): NodeFilesDataSource { return new NodeFilesDataSource(this.client); }

  internal(route: string): InternalDataSource { return new InternalDataSource(this.client, route); }

  async serverInfo(): Promise<{version: string; commit?: string}> {
    const raw = await this.client.get('/info/server');
    const info = typeof raw === 'string' ? JSON.parse(raw) : raw;
    return {version: info?.Version ?? '', commit: info?.Commit};
  }

  async raw(method: string, path: string, body?: any): Promise<any> {
    return this.client.request(method.toUpperCase(), path, body);
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
