
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
      let errBody: any;
      try { errBody = await res.json(); } catch { errBody = {error: await res.text()}; }
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

export class NodeFuncsDataSource extends NodeHttpDataSource {
  async run(name: string, params?: Record<string, any>): Promise<any> {
    const normalizedName = name.replace(':', '.');
    return this.client.post(`/public/v1/functions/${encodeURIComponent(normalizedName)}/call`, params ?? {});
  }
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

  get users(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'users'); }
  get groups(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'groups'); }
  get functions(): NodeFuncsDataSource { return new NodeFuncsDataSource(this.client, 'functions'); }
  get connections(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'connections'); }
  get queries(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'queries'); }
  get scripts(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'scripts'); }
  get packages(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'packages'); }
  get reports(): NodeHttpDataSource { return new NodeHttpDataSource(this.client, 'reports'); }
  get files(): NodeFilesDataSource { return new NodeFilesDataSource(this.client); }

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

  async describe(entityType: string): Promise<any> {
    try {
      return await this.client.get(`/public/v1/entity-types/${encodeURIComponent(entityType)}`);
    } catch {
      return await this.client.get(`/entities/types${buildQuery({name: entityType})}`);
    }
  }
}
