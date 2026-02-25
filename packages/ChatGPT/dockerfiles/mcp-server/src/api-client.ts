const BASE_URL = process.env.DATAGROK_API_URL ?? 'http://localhost:8082/api';
const API_KEY = process.env.DATAGROK_API_KEY ?? '';

interface Connection {
  id: string;
  name: string;
}

function toGrokPath(name: string): string {
  return name.replace(/:/g, '.');
}

async function request<T = unknown>(
  method: string, path: string, body?: unknown, headers?: Record<string, string>,
): Promise<T> {
  const url = `${BASE_URL}/${path.replace(/^\//, '')}`;
  const res = await fetch(url, {
    method,
    headers: {
      'Content-Type': 'application/json',
      'Authorization': API_KEY,
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

export function listFunctions(filter?: string): Promise<unknown[]> {
  const qs = filter ? `?text=${encodeURIComponent(filter)}` : '';
  return request('GET', `/public/v1/functions${qs}`);
}

export function getFunction(id: string): Promise<unknown> {
  return request('GET', `/public/v1/functions/${toGrokPath(id)}`);
}

export function callFunction(name: string, params?: Record<string, unknown>): Promise<unknown> {
  return request('POST', `/public/v1/functions/${toGrokPath(name)}/call`, params ?? {});
}

export function saveFunction(data: Record<string, unknown>): Promise<unknown> {
  return request('POST', '/public/v1/functions', data);
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

export function getCurrentUser(): Promise<unknown> {
  return request('GET', '/public/v1/users/current');
}
