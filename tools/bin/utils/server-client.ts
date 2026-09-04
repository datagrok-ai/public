import {NodeApiClient} from './node-dapi';
import {getDevKey} from './test-utils';

/**
 * `--admin` asks the server for an admin session, which lifts the permission filter for this run:
 * without it a stand-wide pull sees only what the key's own account can, and content in other
 * people's spaces is invisible. The server authorises it (`START_ADMIN_SESSION`) and signs the
 * flag into the token, and a dev key mints a fresh session per invocation, so it cannot outlive
 * the command or reach another session.
 */
export async function createClient(hostArg?: string, admin: boolean = false): Promise<NodeApiClient> {
  const {url, key} = getDevKey(hostArg ?? '');
  const client = await NodeApiClient.login(url, key);
  if (!admin)
    return client;
  const token = (await client.post('/users/sessions/current/admin'))?.token;
  if (!token)
    throw new Error(`${url} refused an admin session — the account behind this key cannot start one`);
  client.token = token;
  client.adminMode = true;
  return client;
}
