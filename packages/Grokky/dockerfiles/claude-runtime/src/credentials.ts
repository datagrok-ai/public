import {createHash} from 'node:crypto';
import {PROVIDER_CREDENTIALS, rewriteForDocker} from './constants';

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

const API_URL = process.env.DATAGROK_API_URL && rewriteForDocker(process.env.DATAGROK_API_URL);

const NOT_SET = 'No AI model credentials are set for you: an administrator sets them in Plugins | Grokky | Credentials...';

/** The provider credentials of the user whose proxied connection brought `token`: Grokky's package credentials,
 * resolved by the platform by the user's groups. */
export async function userCredentials(token: string | undefined): Promise<ProviderCredentials> {
  const failed = (why: string) => new Error(`The assistant could not read its credentials: ${why}`);
  if (!API_URL)
    throw failed('DATAGROK_API_URL is not set');
  if (!token)
    throw failed('the connection did not come through Datagrok');
  const res = await fetch(`${API_URL}/docker/credentials`, {
    headers: {'Authorization': `Bearer ${token}`},
    signal: AbortSignal.timeout(10_000),
  }).catch((e: Error) => {
    throw failed(e.message);
  });
  const status = Number(res.headers.get('api-error-code') ?? res.status);
  if (status === 404)
    throw new Error(NOT_SET);
  if (!res.ok || res.headers.has('api-error'))
    throw failed(`HTTP ${status} from /docker/credentials`);
  const all: Record<string, unknown> = await res.json();
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
