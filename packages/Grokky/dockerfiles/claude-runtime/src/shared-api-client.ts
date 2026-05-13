// Duplicated in mcp-server/src/shared-api-client.ts — keep in sync.
// Can't be shared because each Docker container builds with its own context.
import {AsyncLocalStorage} from 'node:async_hooks';

interface RequestContext {
  apiKey: string;
  apiUrl: string;
}

const ctxStore = new AsyncLocalStorage<RequestContext>();

function _ctx(): RequestContext {
  const ctx = ctxStore.getStore();
  if (!ctx)
    throw new Error('No request context — missing x-user-api-key / x-datagrok-api-url headers');
  return ctx;
}

export function runWithContext<T>(ctx: RequestContext, fn: () => T): T {
  return ctxStore.run(ctx, fn);
}

export interface Connection {
  id: string;
  name: string;
}

export function toGrokPath(name: string): string {
  return name.replace(/:/g, '.');
}

export async function request<T = unknown>(
  method: string, path: string, body?: unknown, headers?: Record<string, string>,
): Promise<T> {
  const ctx = _ctx();
  const url = `${ctx.apiUrl}/${path.replace(/^\//, '')}`;
  const res = await fetch(url, {
    method,
    headers: {
      'Content-Type': 'application/json',
      'Authorization': ctx.apiKey,
      ...headers,
    },
    body: body !== undefined ? JSON.stringify(body) : undefined,
  });

  if (!res.ok) {
    const text = await res.text();
    throw new Error(`HTTP ${res.status}: ${text}`);
  }
  if (res.headers.has('api-error')) {
    const text = await res.text();
    throw new Error(`API error: ${text}`);
  }

  const ct = res.headers.get('content-type') ?? '';
  if (ct.includes('application/json'))
    return res.json();
  return res.text() as T;
}

export async function requestBinary(
  method: string, path: string,
): Promise<Buffer> {
  const ctx = _ctx();
  const url = `${ctx.apiUrl}/${path.replace(/^\//, '')}`;
  const res = await fetch(url, {
    method,
    headers: {'Authorization': ctx.apiKey},
  });
  if (!res.ok) {
    const text = await res.text();
    throw new Error(`HTTP ${res.status}: ${text}`);
  }
  return Buffer.from(await res.arrayBuffer());
}

export async function listFiles(connector: string, path?: string): Promise<unknown[]> {
  const searchName = connector.includes(':') ? connector.split(':').pop()! : connector;
  const conns = await request<Connection[]>(
    'GET', `/public/v1/connections?text=${encodeURIComponent(searchName)}`,
  );
  if (!conns.length)
    throw new Error(`Connection "${connector}" not found`);
  const filePath = path ? `${path.replace(/^\//, '')}/` : '';
  return request('GET', `/connectors/connections/${conns[0].id}/files/${filePath}`);
}

export function downloadFile(connector: string, path: string): Promise<unknown> {
  return request('GET', `/public/v1/files/${toGrokPath(connector)}/${path}`, undefined,
    {'Accept': 'application/octet-stream'});
}

export function uploadFile(connector: string, path: string, content: string): Promise<unknown> {
  return request('POST', `/public/v1/files/${toGrokPath(connector)}/${path}`, content,
    {'Content-Type': 'application/octet-stream'});
}
