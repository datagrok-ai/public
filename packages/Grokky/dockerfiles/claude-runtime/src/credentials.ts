import {createHash} from 'node:crypto';
import {PROVIDER_CREDENTIALS, PROVIDER_ENV, rewriteForDocker} from './constants';

export const CONTAINER_TOKEN_HEADER = 'x-datagrok-container-token';

/** Provider credentials an agent process starts with; `id` tells two sets apart without holding the values. */
export interface ProviderCredentials {
  values: Record<string, string>;
  id: string;
}

export function providerCredentials(values: Record<string, string>): ProviderCredentials {
  const id = createHash('sha256').update(JSON.stringify(Object.entries(values).sort())).digest('hex');
  return {values, id};
}

export const containerCredentials = providerCredentials(PROVIDER_ENV);

const API_URL = process.env.DATAGROK_API_URL && rewriteForDocker(process.env.DATAGROK_API_URL);

/** The provider credentials of the user whose proxied connection brought `token`, resolved by the platform from the
 * package's and this container's credentials by the user's groups; the container's own when it has none for them. */
export async function userCredentials(token: string | undefined): Promise<ProviderCredentials> {
  if (!token || !API_URL)
    return containerCredentials;
  const failed = (why: string) => new Error(`The assistant could not read its credentials: ${why}`);
  const res = await fetch(`${API_URL}/docker/credentials`, {
    headers: {'Authorization': `Bearer ${token}`},
    signal: AbortSignal.timeout(10_000),
  }).catch((e: Error) => {
    throw failed(e.message);
  });
  const status = Number(res.headers.get('api-error-code') ?? res.status);
  if (status === 404)
    return containerCredentials;
  if (!res.ok || res.headers.has('api-error'))
    throw failed(`HTTP ${status} from /docker/credentials`);
  const all: Record<string, unknown> = await res.json();
  // The platform lays the container's credentials over the package's key by key: the container's key goes out only
  // with the container's own provider settings, never to an endpoint a group's credentials name.
  if (containerCredentials.values.apiKey && all.apiKey === containerCredentials.values.apiKey)
    return containerCredentials;
  const values: Record<string, string> = {};
  for (const k of PROVIDER_CREDENTIALS) {
    const v = all[k];
    if (typeof v === 'string' && v)
      values[k] = v;
  }
  return Object.keys(values).length ? providerCredentials(values) : containerCredentials;
}
