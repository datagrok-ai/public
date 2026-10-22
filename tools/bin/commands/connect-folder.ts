import fs from 'fs';
import os from 'os';
import path from 'path';
import * as color from '../utils/color-utils';
import {generateKeyPair, getServerCredentials, Jwk, ServerCredentials} from '../utils/keypair';
import {createClient} from '../utils/server-client';
import {NodeApiClient, NodeDapi} from '../utils/node-dapi';
import {Indexable} from '../utils/utils';
import {printError, printOutput, setOutputFormat, OutputFormat} from '../utils/server-output';
import * as dk from '../utils/connect-folder/docker';

/**
 * `grok connect-folder <folder>` - makes a local folder a Datagrok connection. The connection is created over REST
 * with the keypair session, and served by the proxy-storage agent in a Docker container that dials out
 * to the server; `list`, `start`, `stop`, `rm` and `logs` manage the containers of this machine.
 */
export async function connectFolder(args: ShareArgs): Promise<boolean> {
  const follow = flag(args, 'f') || args.follow === true;
  const output: OutputFormat = (args.output ?? args.o ?? 'table') as OutputFormat;
  setOutputFormat(output);
  const [, sub, ref, ...extra] = args._.map((a) => `${a}`);
  if (sub == null || extra.length > 0)
    return false;
  try {
    switch (sub) {
      case 'list':
        return ref == null && await list(args.host, output);
      case 'start':
      case 'stop':
      case 'logs':
        return ref != null && await manage(sub, ref, args.host, follow);
      case 'rm':
        return ref != null && await remove(ref, args.host);
      default:
        return ref == null && await create(sub, {host: args.host, image: args.image,
          readOnly: args['read-only'] === true, foreground: args.foreground === true, force: args.force === true});
    }
  } catch (e: any) {
    printError(e);
    return false;
  }
}

interface CreateOptions { host?: string, image?: string, readOnly: boolean, foreground: boolean, force: boolean }

async function create(folderArg: string, opts: CreateOptions): Promise<boolean> {
  const cred = shareCredentials(opts.host);
  await dk.dockerAvailable();
  let folder: string;
  let note: string | undefined;
  try {
    ({path: folder, note} = dk.hostFolder(folderArg));
  } catch (e: any) {
    throw e.code === 'ENOENT' ? new Error(`${folderArg} does not exist`) : e;
  }
  if (!fs.statSync(folder).isDirectory())
    throw new Error(`${folder} is not a folder`);
  checkFolderSafety(folder, opts.force, os.homedir(), process.platform,
    cred.privateKeySource && cred.privateKeySource !== 'GROK_PRIVATE_KEY' ? [cred.privateKeySource] : []);
  if (note)
    color.warn(note);

  const alias = cred.alias ?? cred.url;
  const client = await createClient(alias);
  const before: any[] = await client.get('/local_shares');
  const created = await color.step(`Connecting ${folder} to ${cred.url}`, () => client.post('/local_shares', {
    host: dk.canonicalHost(), hostLabel: os.hostname(), folder: dk.canonicalFolder(folder), folderLabel: folder,
    readOnly: opts.readOnly,
  })).catch((e: any) => {
    throw e?.apiError?.errorCode === 403 && `${e.apiError.error}`.includes('ConnectLocalFolders')
      ? new Error(`You don't have permission to connect folders on ${new URL(cred.url).host}. Ask an administrator.`) : e;
  });
  const conn = created.connection;
  const name = conn.friendlyName ?? conn.name;
  const spec: dk.RunSpec = {
    image: opts.image ?? created.agentImage, connectionId: conn.id, alias, apiRoot: cred.url, name, folder,
    hostLabel: os.hostname(), readOnly: opts.readOnly, envFile: '', foreground: opts.foreground, user: dk.linuxUser(),
  };

  const containers = await dk.ps({[`${dk.LABEL}.connection`]: conn.id, [`${dk.LABEL}.url`]: cred.url});
  const existing = containers.find((c) => c.name === dk.containerName(cred.url, conn.id));
  const keys = await shareKeys(client, conn.id);
  if (existing && !opts.foreground && dk.sameSettings(existing, spec) && keys.some((k) => k.id === existing.labels[`${dk.LABEL}.key`])) {
    // A share deleted on the server leaves its agent idle (4410), so the recreated one needs a fresh start.
    const wasDeleted = !before.some((s) => s.id === conn.id);
    if (existing.state === 'running' && !wasDeleted)
      color.info(`"${name}" is already connected from this machine (container ${existing.name})`);
    else
      await color.step(`Starting container ${existing.name}`, () => dk.dockerOk([wasDeleted ? 'restart' : 'start', existing.name]));
    return await reportOnline(client, conn.id, name, existing.name);
  }
  for (const c of containers)
    await color.step(`Replacing container ${c.name}`, () => dk.rm(c.name));

  // A key is registered only once the image is here: a Ctrl+C during a long pull then leaves nothing behind.
  await color.step(`Getting the folder agent image ${spec.image}`, () => dk.ensureImage(spec.image));
  const key = await issueShareKey(client, conn.id);
  const revokeMine = () => client.del(`/users/keys/${key.id}`).catch(() => {});
  spec.keyId = key.id;
  const container = dk.containerName(cred.url, conn.id, opts.foreground);
  if (opts.foreground)
    return await runForeground(client, spec, key.privateKey, container, name, revokeMine);
  try {
    await withKeyFile(key.privateKey, (envFile) => color.step('Starting the folder agent',
      () => dk.dockerOk(dk.runArgs({...spec, envFile}))), revokeMine);
  } catch (e) {
    await revokeMine();
    throw e;
  }
  if (await dk.isRunning(container))
    await revokeShareKeys(client, conn.id, key.id);
  if (process.platform === 'win32' || process.platform === 'darwin')
    color.warn('The folder comes back after a reboot only when Docker does: let Docker Desktop start when you sign in.');
  return await reportOnline(client, conn.id, name, container);
}

/** `docker create` reads the key file, which is deleted before the agent starts; `start -a` forwards Ctrl+C. */
async function runForeground(client: NodeApiClient, spec: dk.RunSpec, key: Jwk, container: string,
  name: string, revokeMine: () => Promise<void>): Promise<boolean> {
  try {
    await withKeyFile(key, (envFile) => dk.dockerOk(['create', ...dk.runArgs({...spec, envFile}).slice(1)]), revokeMine);
    color.info(`Serving "${name}" until Ctrl+C`);
    process.on('SIGINT', () => { dk.docker(['stop', container]).catch(() => {}); });
    const attached = dk.docker(['start', '-a', container], true);
    if (await dk.isRunning(container, 10))
      await revokeShareKeys(client, spec.connectionId, spec.keyId).catch(() => {});
    const r = await attached;
    return r.code === 0 || r.code === 130 || r.code === 143;
  } finally {
    await revokeMine();
  }
}

/**
 * [key], base64 in a 0600 `--env-file` that exists only while [action] runs; Ctrl+C meanwhile deletes it and
 * runs [onAbort] before exiting.
 */
export async function withKeyFile<T>(key: Jwk, action: (envFile: string) => Promise<T>,
  onAbort: () => Promise<void> = async () => {}): Promise<T> {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-folder-'));
  const file = path.join(dir, 'key.env');
  const abort = (signal: NodeJS.Signals) => {
    fs.rmSync(dir, {recursive: true, force: true});
    onAbort().finally(() => process.exit(signal === 'SIGINT' ? 130 : 143));
  };
  process.once('SIGINT', abort);
  process.once('SIGTERM', abort);
  try {
    fs.writeFileSync(file, `GROK_SHARE_KEY=${Buffer.from(JSON.stringify(key)).toString('base64')}\n`, {mode: 0o600});
    return await action(file);
  } finally {
    process.off('SIGINT', abort);
    process.off('SIGTERM', abort);
    fs.rmSync(dir, {recursive: true, force: true});
  }
}

/** The container never gets the user's own key: the server issues each share one that can only serve it. */
async function issueShareKey(client: NodeApiClient, connectionId: string): Promise<{id: string, privateKey: Jwk}> {
  const {publicKey, privateKey} = generateKeyPair();
  const key = await client.post(`/local_shares/${connectionId}/agent_key`, {publicKey: JSON.stringify(publicKey)});
  return {id: key.id, privateKey};
}

const folderScope = (connectionId: string) => `connected-folder:${connectionId}`;

async function shareKeys(client: NodeApiClient, connectionId: string): Promise<any[]> {
  return ((await client.get('/users/current/keys')) as any[] ?? []).filter((k) => k.scope === folderScope(connectionId));
}

async function revokeShareKeys(client: NodeApiClient, connectionId: string, keep?: string): Promise<void> {
  for (const k of await shareKeys(client, connectionId)) {
    if (k.id !== keep)
      await client.del(`/users/keys/${k.id}`);
  }
}

async function reportOnline(client: NodeApiClient, id: string, name: string, container: string): Promise<boolean> {
  const online = await color.step('Waiting for the agent to connect', async () => {
    for (let i = 0; i < 30; i++) {
      const s = ((await client.get('/local_shares')) as any[]).find((s) => s.id === id);
      if (s != null && s.isAvailable !== false)
        return true;
      await new Promise((r) => setTimeout(r, 1000));
    }
    return false;
  });
  if (!online) {
    color.warn(`"${name}" is connected, but its agent has not dialed in yet. See: grok connect-folder logs ${id.slice(0, 8)}`);
    return true;
  }
  color.success(`"${name}" is online (connection ${id}, container ${container})`);
  return true;
}

async function list(host: string | undefined, output: OutputFormat): Promise<boolean> {
  const containers = await dockerContainersOrEmpty();
  const own = getServerCredentials(host ?? '');
  const aliases = new Set<string>([own.alias ?? own.url]);
  if (host == null) {
    for (const c of containers)
      aliases.add(c.labels[`${dk.LABEL}.alias`]);
  }
  const rows: Indexable[] = [];
  const seen = new Set<string>();
  for (const alias of aliases) {
    let shares: any[] | null = null;
    let keys: any[] = [];
    let url: string;
    try {
      const cred = getServerCredentials(alias);
      url = cred.url;
      const client = await createClient(alias);
      shares = await client.get('/local_shares');
      keys = (await client.get('/users/current/keys')) ?? [];
    } catch (e: any) {
      color.warn(`${alias}: ${e.message ?? e}`);
      url = '';
    }
    const ofServer = containers.filter((c) => c.labels[`${dk.LABEL}.url`] === url || c.labels[`${dk.LABEL}.alias`] === alias);
    for (const s of shares ?? []) {
      const c = ofServer.find((c) => c.labels[`${dk.LABEL}.connection`] === s.id);
      seen.add(`${alias}|${s.id}`);
      // The server refuses an agent whose key is not limited to its folder: its container runs, but serves nothing.
      const outdated = c != null && !keys.some((k) => k.id === c.labels[`${dk.LABEL}.key`] && k.scope === folderScope(s.id));
      rows.push(row(s.friendlyName ?? s.name, s.parameters?.localShareFolder, s.parameters?.localShareReadOnly === 'true',
        c?.state, outdated ? 'outdated key: re-run connect-folder' : s.isAvailable === false ? 'offline' : 'online',
        s.parameters?.localShareHost, s.id, alias));
    }
    for (const c of ofServer) {
      const id = c.labels[`${dk.LABEL}.connection`];
      if (seen.has(`${alias}|${id}`))
        continue;
      seen.add(`${alias}|${id}`);
      rows.push(row(c.labels[`${dk.LABEL}.name`], c.labels[`${dk.LABEL}.folder`], c.labels[`${dk.LABEL}.readonly`] === 'true',
        c.state, shares == null ? 'unknown' : 'revoked', os.hostname(), id, alias));
    }
  }
  printOutput(rows, output);
  return true;
}

function row(name: string, folder: string, readOnly: boolean, container: string | undefined, server: string,
  host: string, id: string, alias: string): Indexable {
  return {name, folder, access: readOnly ? 'read-only' : 'read-write', container: container ?? '-', server, host, id, alias};
}

async function dockerContainersOrEmpty(): Promise<dk.ShareContainer[]> {
  try {
    return await dk.ps();
  } catch (e: any) {
    color.warn(`${e.message ?? e} Only the server side is listed.`);
    return [];
  }
}

async function manage(command: 'start' | 'stop' | 'logs', ref: string, host: string | undefined, follow: boolean): Promise<boolean> {
  await dk.dockerAvailable();
  const c = findContainer(inScope(await dk.ps(), host), ref);
  if (command === 'logs')
    return (await dk.logs(c.name, follow)).code === 0;
  const name = c.labels[`${dk.LABEL}.name`];
  await color.step(`${command === 'start' ? 'Starting' : 'Stopping'} "${name}"`,
    () => command === 'start' ? dk.start(c.name) : dk.stop(c.name));
  if (command === 'stop')
    return true;
  let client: NodeApiClient;
  try {
    client = await createClient(c.labels[`${dk.LABEL}.alias`]);
  } catch (e: any) {
    color.warn(`"${name}" is started, but ${c.labels[`${dk.LABEL}.url`]} cannot be reached to check it: ${e.message ?? e}`);
    return true;
  }
  return await reportOnline(client, c.labels[`${dk.LABEL}.connection`], name, c.name);
}

/** The containers of the server [host] names (all of them without one), by alias or by URL. */
function inScope(containers: dk.ShareContainer[], host: string | undefined): dk.ShareContainer[] {
  if (host == null)
    return containers;
  let url: string | undefined;
  try {
    url = getServerCredentials(host).url;
  } catch { /* an alias that is gone can still name its containers */ }
  return containers.filter((c) => c.labels[`${dk.LABEL}.alias`] === host || c.labels[`${dk.LABEL}.url`] === url);
}

/**
 * Removes the share's containers on this machine, then its connection and keys on the server. The containers go
 * even when the server cannot be reached.
 */
async function remove(ref: string, host: string | undefined): Promise<boolean> {
  const containers = await dockerContainersOrEmpty();
  const found = matchContainers(inScope(containers, host), ref);
  const ids = [...new Set(found.map((c) => c.labels[`${dk.LABEL}.connection`]))];
  if (new Set(found.map((c) => `${c.labels[`${dk.LABEL}.alias`]}|${c.labels[`${dk.LABEL}.connection`]}`)).size > 1)
    throw ambiguous(ref, found.map((c) => `${c.labels[`${dk.LABEL}.name`]} (${c.labels[`${dk.LABEL}.alias`]})`));
  const alias = host ?? found[0]?.labels[`${dk.LABEL}.alias`] ?? '';
  let id: string | undefined = ids[0];
  let client: NodeApiClient | undefined;
  let serverError: any;
  try {
    client = await createClient(alias);
    if (id == null) {
      const shares = (await client.get('/local_shares') as any[]).filter((s) => matches(ref,
        s.id, s.friendlyName ?? s.name, s.parameters?.localShareFolder));
      if (shares.length > 1)
        throw ambiguous(ref, shares.map((s) => s.friendlyName ?? s.name));
      id = shares[0]?.id;
    }
  } catch (e: any) {
    serverError = new Error(`${alias || 'the default server'}: ${e?.message ?? e}`);
  }
  if (id == null)
    throw serverError ?? new Error(`No connected folder matches "${ref}". See: grok connect-folder list`);

  for (const c of inScope(containers, alias).filter((c) => c.labels[`${dk.LABEL}.connection`] === id))
    await color.step(`Removing container ${c.name}`, () => dk.rm(c.name));
  if (client && !serverError) {
    let doing = 'revoke the folder key';
    try {
      await color.step('Revoking the folder key', () => revokeShareKeys(client!, id!));
      doing = 'delete the connection';
      await color.step(`Deleting connection ${id}`, () => new NodeDapi(client!).connections.delete({id: id!}));
    } catch (e: any) {
      serverError = new Error(`could not ${doing}: ${e?.message ?? e}`);
    }
  }
  if (serverError) {
    color.error(`Removed on this machine, but the connection is still on the server: ${serverError.message ?? serverError}. ` +
      `Run grok connect-folder rm ${id.slice(0, 8)} again, or delete the connection in Datagrok`);
    return false;
  }
  return true;
}

function findContainer(containers: dk.ShareContainer[], ref: string): dk.ShareContainer {
  const found = matchContainers(containers, ref);
  if (found.length === 0)
    throw new Error(`No folder container on this machine matches "${ref}". See: grok connect-folder list`);
  if (found.length > 1)
    throw ambiguous(ref, found.map((c) => `${c.labels[`${dk.LABEL}.name`]} (${c.labels[`${dk.LABEL}.alias`]})`));
  return found[0];
}

function matchContainers(containers: dk.ShareContainer[], ref: string): dk.ShareContainer[] {
  return containers.filter((c) => matches(ref, c.labels[`${dk.LABEL}.connection`], c.labels[`${dk.LABEL}.name`],
    c.labels[`${dk.LABEL}.folder`]));
}

/** A share's `<ref>`: its name, a connection id prefix of 8 or more characters, or its folder. */
export function matches(ref: string, id: string, name: string | undefined, folder: string | undefined): boolean {
  if (name != null && name.toLowerCase() === ref.toLowerCase())
    return true;
  if (ref.length >= 8 && id?.toLowerCase().startsWith(ref.toLowerCase()))
    return true;
  if (folder == null)
    return false;
  let abs = path.resolve(ref);
  try {
    abs = fs.realpathSync.native(abs);
  } catch { /* a folder that is gone can still name its share */ }
  return dk.canonicalFolder(abs) === dk.canonicalFolder(folder);
}

function ambiguous(ref: string, names: string[]): Error {
  return new Error(`"${ref}" matches several connected folders: ${names.join(', ')}. Use the connection id, or --host`);
}

/** The server alias with the keypair `grok connect-folder` signs in with (the agent gets a key of its own). */
function shareCredentials(host: string | undefined): ServerCredentials {
  const cred = getServerCredentials(host ?? '');
  if (cred.privateKey == null)
    throw new Error(`grok connect-folder signs in with a keypair, and ${cred.alias ?? cred.url} has none: run grok login ${cred.url}`);
  return cred;
}

/**
 * Refuses a folder that is, holds or is inside `~/.grok` or one of [keyFiles] (private keys); a filesystem or
 * drive root, or the home folder, needs [force].
 */
export function checkFolderSafety(folder: string, force: boolean, home: string = os.homedir(),
  platform: string = process.platform, keyFiles: string[] = []): void {
  const f = dk.canonicalFolder(folder, platform);
  const within = (child: string, parent: string) => child === parent || child.startsWith(parent.endsWith('/') ? parent : `${parent}/`);
  for (const secret of [path.join(home, '.grok'), ...keyFiles]) {
    const s = dk.canonicalFolder(secret, platform);
    if (within(f, s) || within(s, f))
      throw new Error(`${folder} is or holds ${secret}, where your private keys are; it cannot be connected`);
  }
  if (!force && (f === '/' || /^[a-z]:\/?$/i.test(f) || f === dk.canonicalFolder(home, platform)))
    throw new Error(`${folder} is a root or your home folder; pass --force to connect it anyway`);
}

/** minimist reads `-f <connection>` as a value of the flag, so the value goes back to the positionals. */
function flag(args: ShareArgs, name: string): boolean {
  const v = args[name];
  if (typeof v === 'string' && v !== '') {
    args._.push(v);
    return true;
  }
  return v === true;
}

interface ShareArgs extends Indexable {
  _: string[],
  host?: string,
  image?: string,
  output?: string,
  o?: string,
}
