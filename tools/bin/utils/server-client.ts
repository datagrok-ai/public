import {NodeApiClient, NodeDapi} from './node-dapi';
import {getServerCredentials} from './keypair';
import {printError, OutputFormat} from './server-output';

/** Opens a session on a host alias or URL (the configured default when omitted). */
export type Connect = (host?: string) => Promise<NodeDapi>;

/** `--host` as minimist leaves it: absent, one value, or an array when repeated. */
export function hostList(host: any): string[] {
  if (host === undefined || host === null || host === true) return [];
  return (Array.isArray(host) ? host : [host]).map(String);
}

/** Commands that talk to one deployment refuse a repeated `--host`. */
export function singleHost(argv: any, command: string): string | undefined {
  const hosts = hostList(argv.host);
  if (hosts.length > 1)
    throw new Error(`'grok s ${command}' takes one --host`);
  return hosts[0];
}

/**
 * Runs [fn] against every `--host`. With several hosts a failing one is reported on stderr and
 * the rest still answer; the run then exits 1.
 */
export async function eachHost(argv: any, connect: Connect,
                               fn: (dapi: NodeDapi, host: string, multi: boolean) => Promise<void>): Promise<void> {
  const hosts = hostList(argv.host);
  if (hosts.length < 2)
    return fn(await connect(hosts[0]), hosts[0] ?? '', false);
  let answered = 0;
  for (const host of hosts) {
    try {
      await fn(await connect(host), host, true);
      answered++;
    }
    catch (err: any) {
      const apiError = err?.apiError ? {...err.apiError, error: `${host}: ${err.apiError.error}`} : undefined;
      printError({message: `${host}: ${err?.message ?? err}`, apiError});
      process.exitCode = 1;
    }
  }
  if (!answered)
    throw new Error('No host answered');
}

/**
 * Rows from every `--host`, each prefixed with the host when there are several: a `HOST` column
 * in a table, a `host` field in JSON.
 */
export async function forEachHost(argv: any, connect: Connect, fn: (dapi: NodeDapi, host: string) => Promise<any[]>,
                                  output: OutputFormat): Promise<any[]> {
  const column = output === 'json' ? 'host' : 'HOST';
  const rows: any[] = [];
  await eachHost(argv, connect, async (dapi, host, multi) => {
    for (const row of await fn(dapi, host))
      rows.push(multi ? {[column]: host, ...row} : row);
  });
  return rows;
}

/**
 * `--admin` asks the server for an admin session, which lifts the permission filter for this run:
 * without it a stand-wide pull sees only what the key's own account can, and content in other
 * people's spaces is invisible. The server authorises it (`START_ADMIN_SESSION`) and signs the
 * flag into the token, and a dev key mints a fresh session per invocation, so it cannot outlive
 * the command or reach another session.
 */
export async function createClient(hostArg?: string, admin: boolean = false): Promise<NodeApiClient> {
  // Resolved from the alias the caller named, not from its URL: two aliases can point at the
  // same server, and only this side knows which one was asked for.
  const cred = getServerCredentials(hostArg ?? '');
  const client = await NodeApiClient.login(cred.url, cred.key ?? '', cred.privateKey);
  if (!admin)
    return client;
  const token = (await client.post('/users/sessions/current/admin'))?.token;
  if (!token)
    throw new Error(`${cred.url} refused an admin session — the account behind this key cannot start one`);
  client.token = token;
  client.adminMode = true;
  return client;
}
