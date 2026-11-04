import {createHash} from 'node:crypto';
import {PROVIDER_CREDENTIALS, rewriteForDocker} from './constants';

export const CONTAINER_TOKEN_HEADER = 'x-datagrok-container-token';

/** `id` tells two sets apart without exposing the values. */
export interface ProviderCredentials {
  values: Record<string, string>;
  id: string;
}

export function providerCredentials(values: Record<string, string>): ProviderCredentials {
  return {values, id: createHash('sha256').update(JSON.stringify(Object.entries(values).sort())).digest('hex')};
}

const API_URL = process.env.DATAGROK_API_URL && rewriteForDocker(process.env.DATAGROK_API_URL);
const NOT_SET =
  'No AI model credentials are set for you: an administrator sets them in Plugins | Grokky | Credentials...';
const failed = (why: string) => new Error(`The assistant could not read its credentials: ${why}`);

/** Grokky's package credentials for the user whose connection brought `token`, resolved by their groups. */
export async function userCredentials(token?: string): Promise<ProviderCredentials> {
  if (!API_URL)
    throw failed('DATAGROK_API_URL is not set');
  if (!token)
    throw failed('the connection did not come through Datagrok');
  const res = await fetch(`${API_URL}/docker/credentials`,
    {headers: {Authorization: `Bearer ${token}`}, signal: AbortSignal.timeout(10_000)})
    .catch((e: Error) => Promise.reject(failed(e.message)));
  const status = Number(res.headers.get('api-error-code') ?? res.status);
  if (!res.ok || res.headers.has('api-error')) {
    await res.body?.cancel();
    throw status === 404 ? new Error(NOT_SET) : failed(`HTTP ${status} from /docker/credentials`);
  }
  const all: Record<string, unknown> = await res.json().catch(() => Promise.reject(failed('the answer is not JSON')));
  const values: Record<string, string> = {};
  for (const k of PROVIDER_CREDENTIALS) {
    const v = all[k];
    if (typeof v === 'string' && v)
      values[k] = v;
  }
  if (!Object.keys(values).length)
    throw new Error(NOT_SET);
  return providerCredentials(values);
}
