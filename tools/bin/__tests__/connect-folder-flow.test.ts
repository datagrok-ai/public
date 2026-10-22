import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import {EventEmitter} from 'events';
import fs from 'fs';
import os from 'os';
import path from 'path';
import * as kp from '../utils/keypair';
import * as dk from '../utils/connect-folder/docker';
import {connectFolder as share, withKeyFile} from '../commands/connect-folder';

const A = 'http://stand-a/api';
const B = 'http://stand-b/api';
const ID = '7ff0714e-1d2c-5a3b-9e8f-0123456789ab';
const USER_KEY = kp.generateKeyPair().privateKey;

interface FakeContainer { id: string, name: string, state: string, image: string, labels: Record<string, string>, rm: boolean }

/** A `docker` CLI over an in-memory container list; `run`/`create` also record what the key file held at that moment. */
const docker = {
  containers: [] as FakeContainer[],
  calls: [] as string[][],
  failRun: false,
  /** `run` starts a container that exits at once. */
  runExits: false,
  images: new Set<string>(),
  keyFiles: [] as {path: string, existed: boolean, content: string, mode: number}[],
};

/** Docker and server calls in the order they happened. */
const events: string[] = [];

function handle(args: string[]): {code: number, out?: string, err?: string, delay?: number, after?: () => void} {
  const byName = (n: string) => docker.containers.find((c) => c.name === n || c.id === n);
  switch (args[0]) {
    case 'version': return {code: 0, out: '27.0.0'};
    case 'info': return {code: 0, out: 'linux\n'};
    case 'ps': {
      const filters = args.filter((_, i) => args[i - 1] === '--filter').map((f) => f.slice('label='.length));
      return {code: 0, out: docker.containers.filter((c) => filters.every((f) => {
        const eq = f.indexOf('=');
        return c.labels[f.slice(0, eq)] === f.slice(eq + 1);
      })).map((c) => c.id).join('\n')};
    }
    case 'image':
      return docker.images.has(args[4]) ? {code: 0, out: 'sha256:1'} : {code: 1, err: 'No such image'};
    case 'pull':
      docker.images.add(args[1]);
      return {code: 0};
    case 'inspect':
      if (args[1] === '--format') {
        const c = byName(args[3]);
        return c ? {code: 0, out: `${c.state === 'running'}\n`} : {code: 1, err: 'No such object'};
      }
      return {code: 0, out: JSON.stringify(args.slice(1).map(byName).map((c) => ({Id: c!.id, Name: `/${c!.name}`,
        State: {Status: c!.state}, Config: {Image: c!.image, Labels: c!.labels}, Mounts: []})))};
    case 'run':
    case 'create': {
      const valued = ['--name', '--label', '--mount', '--env-file', '--cap-drop', '--security-opt', '--memory', '-e',
        '--user', '--add-host', '--restart'];
      const labels: Record<string, string> = {};
      let name = '';
      let image = '';
      for (let i = 1; i < args.length && !image; i++) {
        if (valued.includes(args[i])) {
          const v = args[++i];
          if (args[i - 1] === '--name')
            name = v;
          else if (args[i - 1] === '--label')
            labels[v.slice(0, v.indexOf('='))] = v.slice(v.indexOf('=') + 1);
          else if (args[i - 1] === '--env-file') {
            const existed = fs.existsSync(v);
            docker.keyFiles.push({path: v, existed, content: existed ? fs.readFileSync(v, 'utf-8') : '',
              mode: existed ? fs.statSync(v).mode & 0o777 : 0});
          }
        } else if (!args[i].startsWith('-'))
          image = args[i];
      }
      if (docker.failRun)
        return {code: 125, err: 'docker: Error response from daemon: pull access denied'};
      if (byName(name))
        return {code: 125, err: `Conflict. The container name "/${name}" is already in use`};
      docker.containers.push({id: `c${docker.containers.length + 1}${Date.now()}`, name, image, labels,
        state: args[0] === 'run' && !docker.runExits ? 'running' : args[0] === 'run' ? 'exited' : 'created',
        rm: args.includes('--rm')});
      return {code: 0, out: 'id'};
    }
    case 'start':
    case 'restart': {
      const c = byName(args[args.length - 1]);
      if (!c)
        return {code: 1, err: 'No such container'};
      c.state = 'running';
      return args.includes('-a') ? {code: 0, delay: 50, after: () => {
        c.state = 'exited';
        if (c.rm)
          docker.containers = docker.containers.filter((x) => x !== c);
      }} : {code: 0};
    }
    case 'stop': {
      const c = byName(args[1]);
      if (c)
        c.state = 'exited';
      return {code: c ? 0 : 1};
    }
    case 'rm':
      docker.containers = docker.containers.filter((c) => c.name !== args[2]);
      return {code: 0};
    default:
      return {code: 0};
  }
}

dk.proc.spawn = ((_cmd: string, args: string[]) => {
  docker.calls.push(args);
  events.push(`docker ${args[0]}`);
  const child: any = new EventEmitter();
  child.stdout = new EventEmitter();
  child.stderr = new EventEmitter();
  setImmediate(() => {
    const r = handle(args);
    if (r.out)
      child.stdout.emit('data', r.out);
    if (r.err)
      child.stderr.emit('data', r.err);
    setTimeout(() => {
      r.after?.();
      child.emit('close', r.code);
    }, r.delay ?? 0);
  });
  return child;
}) as any;

/** Two stands that know the same user; a share's id depends on the login, host and folder only. */
interface Stand {
  shares: Map<string, any>, keys: any[], down: boolean, agentImage: string,
  fail?: (method: string, path: string) => boolean,
  /** `POST /local_shares` answers 403 with this message. */
  forbidden?: string,
  /** The agent shows up online only after this many more share listings. */
  offlinePolls?: number,
}
const stands: Record<string, Stand> = {};

const json = (body: any, status: number = 200) =>
  new Response(JSON.stringify(body), {status, headers: {'content-type': 'application/json'}});

async function serve(url: string, init: RequestInit = {}): Promise<Response> {
  const root = [A, B].find((r) => url.startsWith(r))!;
  const s = stands[root];
  if (s.down)
    throw new Error('connect ECONNREFUSED');
  const p = url.slice(root.length).split('?')[0];
  const method = init.method ?? 'GET';
  events.push(`${method} ${p}`);
  if (s.fail?.(method, p))
    return json({error: 'Internal error'}, 500);
  const body = init.body ? JSON.parse(init.body as string) : null;
  if (s.forbidden && p === '/local_shares' && method === 'POST') {
    return new Response(JSON.stringify({'#type': 'ApiError', message: s.forbidden}),
      {status: 200, headers: {'content-type': 'application/json', 'api-error-code': '403'}});
  }
  if (p === '/users/login/key/challenge')
    return json({nonce: 'n'});
  if (p === '/users/login/key')
    return json({isSuccess: true, token: 'token'});
  if (p === '/local_shares' && method === 'GET') {
    const online = !(s.offlinePolls && s.offlinePolls-- > 0);
    return json([...s.shares.values()].map((c) => ({...c, isAvailable: online})));
  }
  if (p === '/local_shares' && method === 'POST') {
    const conn = {id: ID, friendlyName: 'laptop data', isAvailable: true,
      parameters: {localShareHost: 'laptop', localShareFolder: body.folderLabel, localShareReadOnly: `${body.readOnly}`}};
    s.shares.set(ID, conn);
    return json({connection: conn, agentImage: s.agentImage});
  }
  if (p === '/users/current/keys' && method === 'GET')
    return json(s.keys);
  if (p === '/users/current/keys' && method === 'POST') {
    const key = {id: `key-${s.keys.length + 1}-${Math.random()}`, name: body.name, publicKey: body.publicKey};
    s.keys.push(key);
    return json(key);
  }
  const agentKey = p.match(/^\/local_shares\/([^/]+)\/agent_key$/);
  if (agentKey && method === 'POST') {
    const key = {id: `key-${s.keys.length + 1}-${Math.random()}`, publicKey: body.publicKey,
      name: `grok connect-folder laptop ${s.shares.get(agentKey[1])?.parameters?.localShareFolder} [${agentKey[1]}]`,
      scope: `connected-folder:${agentKey[1]}`};
    s.keys.push(key);
    return json(key);
  }
  if (p.startsWith('/users/keys/') && method === 'DELETE') {
    s.keys = s.keys.filter((k) => k.id !== p.slice('/users/keys/'.length));
    return json('ok');
  }
  if (p.startsWith('/public/v1/connections/') && method === 'DELETE') {
    const id = decodeURIComponent(p.slice('/public/v1/connections/'.length));
    if (!s.shares.delete(id))
      return json({error: 'Connection not found'}, 404);
    return json('ok');
  }
  return json({error: `unexpected ${method} ${p}`}, 500);
}

let folder: string;
const args = (...rest: any[]) => {
  const flags = typeof rest[rest.length - 1] === 'object' ? rest.pop() : {};
  return {_: ['connect-folder', ...rest], ...flags};
};
const containersOf = (root: string) => docker.containers.filter((c) => c.labels['ai.datagrok.folder.url'] === root);
const ran = () => docker.calls.filter((c) => c[0] === 'run' || c[0] === 'create');

beforeEach(() => {
  folder = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-folder-test-'));
  docker.containers = [];
  docker.calls = [];
  docker.failRun = false;
  docker.runExits = false;
  docker.images = new Set(['datagrok/proxy_storage:1.0', 'datagrok/proxy_storage:2.0']);
  events.length = 0;
  docker.keyFiles = [];
  for (const r of [A, B])
    stands[r] = {shares: new Map(), keys: [], down: false, agentImage: 'datagrok/proxy_storage:1.0'};
  process.env.GROK_PRIVATE_KEY = JSON.stringify(USER_KEY);
  vi.stubGlobal('fetch', vi.fn(serve));
  vi.spyOn(console, 'log').mockImplementation(() => {});
  vi.spyOn(process.stdout, 'write').mockImplementation(() => true);
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
  delete process.env.GROK_PRIVATE_KEY;
  fs.rmSync(folder, {recursive: true, force: true});
});

describe('grok connect-folder <folder>', () => {
  it('runs the agent with a key of its own, in a 0600 file deleted right after docker run', async () => {
    expect(await share(args(folder, {host: A}))).toBe(true);
    expect(ran()).toHaveLength(1);
    const [c] = containersOf(A);
    expect(c.name).toBe(dk.containerName(A, ID));
    expect(c.state).toBe('running');

    const keyFile = docker.keyFiles[0];
    expect(keyFile.existed).toBe(true);
    expect(fs.existsSync(keyFile.path)).toBe(false);
    if (process.platform !== 'win32')
      expect(keyFile.mode).toBe(0o600);
    const b64 = keyFile.content.trim().replace(/^GROK_SHARE_KEY=/, '');
    const shareKey = JSON.parse(Buffer.from(b64, 'base64').toString());
    expect(shareKey.d).toBeTruthy();
    expect(shareKey.d).not.toBe(USER_KEY.d);
    for (const a of ran()[0])
      expect(a.includes(b64) || a.includes(shareKey.d) || a.includes(USER_KEY.d as string)).toBe(false);

    const keys = stands[A].keys;
    expect(keys).toHaveLength(1);
    expect(keys[0].name).toMatch(new RegExp(`^grok connect-folder .* \\[${ID}\\]$`));
    expect(kp.fingerprint(JSON.parse(keys[0].publicKey))).toBe(kp.fingerprint(kp.publicPart(shareKey)));
    expect(c.labels['ai.datagrok.folder.key']).toBe(keys[0].id);
  });

  it('resumes a running container, starts a stopped one, and restarts one whose share was deleted', async () => {
    await share(args(folder, {host: A}));
    const key = stands[A].keys[0].id;
    docker.calls = [];
    expect(await share(args(folder, {host: A}))).toBe(true);
    expect(docker.calls.map((c) => c[0])).not.toContain('run');
    expect(docker.calls.map((c) => c[0])).not.toContain('start');

    containersOf(A)[0].state = 'exited';
    docker.calls = [];
    await share(args(folder, {host: A}));
    expect(docker.calls.filter((c) => c[0] === 'start')).toEqual([['start', dk.containerName(A, ID)]]);

    stands[A].shares.clear();
    docker.calls = [];
    await share(args(folder, {host: A}));
    expect(docker.calls.filter((c) => c[0] === 'restart')).toEqual([['restart', dk.containerName(A, ID)]]);
    expect(ran()).toHaveLength(0);
    expect(stands[A].keys.map((k) => k.id)).toEqual([key]);
  });

  it('recreates the container, with a new key, when the mode, the image, the folder or the key changed', async () => {
    const changes: [string, () => any][] = [
      ['mode', () => args(folder, {host: A, 'read-only': true})],
      ['server image', () => { stands[A].agentImage = 'datagrok/proxy_storage:2.0'; return args(folder, {host: A, 'read-only': true}); }],
      ['folder', () => {
        containersOf(A)[0].labels['ai.datagrok.folder.folder'] = '/mnt/c/other';
        return args(folder, {host: A, 'read-only': true});
      }],
      ['revoked key', () => { stands[A].keys = []; return args(folder, {host: A, 'read-only': true}); }],
    ];
    await share(args(folder, {host: A}));
    for (const [what, next] of changes) {
      const before = containersOf(A)[0];
      docker.calls = [];
      expect(await share(next()), what).toBe(true);
      expect(docker.calls.filter((c) => c[0] === 'rm'), what).toEqual([['rm', '-f', before.name]]);
      expect(ran(), what).toHaveLength(1);
      const now = containersOf(A);
      expect(now, what).toHaveLength(1);
      expect(stands[A].keys.map((k) => k.id), what).toEqual([now[0].labels['ai.datagrok.folder.key']]);
    }
    expect(containersOf(A)[0].image).toBe('datagrok/proxy_storage:2.0');
    expect(containersOf(A)[0].labels['ai.datagrok.folder.readonly']).toBe('true');
  });

  it('says so when the user may not connect folders, and leaves no key and no container', async () => {
    stands[A].forbidden = 'Permission ConnectLocalFolders is required';
    const stderr = vi.spyOn(process.stderr, 'write').mockImplementation(() => true);
    expect(await share(args(folder, {host: A}))).toBe(false);
    expect(stderr.mock.calls.map((c) => `${c[0]}`).join(''))
      .toContain('You don\'t have permission to connect folders on stand-a. Ask an administrator.');
    expect(stands[A].keys).toHaveLength(0);
    expect(ran()).toHaveLength(0);
  });

  it('passes any other 403 through as the server says it', async () => {
    stands[A].forbidden = 'The share belongs to another user';
    const stderr = vi.spyOn(process.stderr, 'write').mockImplementation(() => true);
    expect(await share(args(folder, {host: A}))).toBe(false);
    const printed = stderr.mock.calls.map((c) => `${c[0]}`).join('');
    expect(printed).toContain('The share belongs to another user');
    expect(printed).not.toContain('permission to connect folders');
    expect(stands[A].keys).toHaveLength(0);
    expect(ran()).toHaveLength(0);
  });

  it('keeps one id on two servers apart', async () => {
    await share(args(folder, {host: A}));
    await share(args(folder, {host: B}));
    expect(containersOf(A)).toHaveLength(1);
    expect(containersOf(B)).toHaveLength(1);
    expect(containersOf(A)[0].name).not.toBe(containersOf(B)[0].name);
    expect(docker.calls.filter((c) => c[0] === 'ps').pop()).toContain(`label=ai.datagrok.folder.url=${B}`);

    docker.calls = [];
    await share(args(folder, {host: A}));
    expect(ran()).toHaveLength(0);
  });

  it('runs a foreground share as create + start -a, and revokes its key when it ends', async () => {
    expect(await share(args(folder, {host: A, foreground: true}))).toBe(true);
    const name = dk.containerName(A, ID, true);
    expect(docker.calls.filter((c) => ['create', 'start', 'run'].includes(c[0])).map((c) => c.slice(0, 3)))
      .toEqual([['create', '--name', name], ['start', '-a', name]]);
    expect(fs.existsSync(docker.keyFiles[0].path)).toBe(false);
    expect(docker.containers).toHaveLength(0);
    expect(stands[A].keys).toHaveLength(0);
  });
});

describe('grok connect-folder start', () => {
  it('waits for the agent to come online, as a new connection does', async () => {
    await share(args(folder, {host: A}));
    await share(args('stop', ID.slice(0, 8)));
    stands[A].offlinePolls = 1;
    events.length = 0;
    expect(await share(args('start', ID.slice(0, 8)))).toBe(true);
    const started = events.indexOf('docker start');
    const polls = events.slice(started).filter((e) => e === 'GET /local_shares');
    expect(started).toBeGreaterThan(-1);
    expect(polls).toHaveLength(2);
    expect((console.log as any).mock.calls.map((c: any[]) => c.join(' ')).join(' ')).toContain('is online');
  });
});

describe('withKeyFile', () => {
  it('deletes the file when the action fails', async () => {
    let file = '';
    await expect(withKeyFile(USER_KEY, async (f) => { file = f; throw new Error('boom'); })).rejects.toThrow('boom');
    expect(file).not.toBe('');
    expect(fs.existsSync(file)).toBe(false);
  });
});

describe('grok connect-folder rm', () => {
  it('removes the container, the connection and the share key', async () => {
    await share(args(folder, {host: A}));
    expect(await share(args('rm', ID.slice(0, 8)))).toBe(true);
    expect(docker.containers).toHaveLength(0);
    expect(stands[A].shares.size).toBe(0);
    expect(stands[A].keys).toHaveLength(0);
  });

  it('asks which server when the same share runs for two, and takes --host', async () => {
    await share(args(folder, {host: A}));
    await share(args(folder, {host: B}));
    expect(await share(args('rm', 'laptop data'))).toBe(false);
    expect(docker.containers).toHaveLength(2);
    expect(await share(args('rm', 'laptop data', {host: B}))).toBe(true);
    expect(containersOf(A)).toHaveLength(1);
    expect(containersOf(B)).toHaveLength(0);
    expect(stands[A].shares.size).toBe(1);
  });

  it('removes the container even when the server cannot be reached, and fails', async () => {
    await share(args(folder, {host: A}));
    stands[A].down = true;
    expect(await share(args('rm', ID.slice(0, 8)))).toBe(false);
    expect(docker.containers).toHaveLength(0);
    expect(stands[A].shares.size).toBe(1);
  });

  it('matches containers by connection id and alias, whatever URL they were created with', async () => {
    await share(args(folder, {host: A}));
    containersOf(A)[0].labels['ai.datagrok.folder.alias'] = B;
    containersOf(A)[0].labels['ai.datagrok.folder.url'] = 'http://172.31.192.1/api';
    expect(await share(args('rm', ID, {host: B}))).toBe(false);
    expect(docker.containers).toHaveLength(0);
  });
});

describe('keys of an interrupted or failed connection', () => {
  it('gets the image before it registers a key', async () => {
    docker.images.clear();
    await share(args(folder, {host: A}));
    const pulled = events.indexOf('docker pull');
    expect(pulled).toBeGreaterThan(-1);
    expect(pulled).toBeLessThan(events.indexOf(`POST /local_shares/${ID}/agent_key`));
  });

  it('deletes the key file and revokes only its own key when docker run fails', async () => {
    await share(args(folder, {host: A}));
    const old = stands[A].keys[0].id;
    docker.failRun = true;
    docker.keyFiles = [];
    expect(await share(args(folder, {host: A, 'read-only': true}))).toBe(false);
    expect(docker.keyFiles[0].existed).toBe(true);
    expect(fs.existsSync(docker.keyFiles[0].path)).toBe(false);
    expect(stands[A].keys.map((k) => k.id)).toEqual([old]);
  });

  it('keeps the older keys while its container is not running', async () => {
    await share(args(folder, {host: A}));
    const old = stands[A].keys[0].id;
    docker.runExits = true;
    await share(args(folder, {host: A, 'read-only': true}));
    expect(stands[A].keys).toHaveLength(2);
    expect(stands[A].keys.map((k) => k.id)).toContain(old);
  });

  it('revokes the older keys of a foreground share once it runs, and its own when it ends', async () => {
    await share(args(folder, {host: A}));
    let during: string[] = [];
    const watch = setInterval(() => {
      if (docker.containers.some((c) => c.name.endsWith('-fg') && c.state === 'running'))
        during = stands[A].keys.map((k) => k.id);
    }, 5);
    await share(args(folder, {host: A, foreground: true}));
    clearInterval(watch);
    expect(during).toHaveLength(1);
    expect(stands[A].keys).toHaveLength(0);
  });

  it('deletes the key file and runs the abort hook on Ctrl+C', async () => {
    const exit = vi.spyOn(process, 'exit').mockImplementation((() => undefined) as any);
    let aborted = false;
    let file = '';
    await withKeyFile(USER_KEY, async (f) => {
      file = f;
      process.emit('SIGINT', 'SIGINT');
      await new Promise((r) => setTimeout(r, 10));
    }, async () => { aborted = true; });
    expect(fs.existsSync(file)).toBe(false);
    expect(aborted).toBe(true);
    expect(exit).toHaveBeenCalledWith(130);
    expect(process.listenerCount('SIGTERM')).toBe(0);
  });
});

describe('grok connect-folder rm, step by step', () => {
  const printed = () => (console.log as any).mock.calls.map((c: any[]) => c.join(' ')).join('\n');

  it('revokes the share keys before it deletes the connection', async () => {
    await share(args(folder, {host: A}));
    events.length = 0;
    await share(args('rm', ID.slice(0, 8)));
    const revoked = events.findIndex((e) => e.startsWith('DELETE /users/keys/'));
    expect(revoked).toBeGreaterThan(-1);
    expect(revoked).toBeLessThan(events.findIndex((e) => e.startsWith('DELETE /public/v1/connections/')));
  });

  it('says which step failed', async () => {
    await share(args(folder, {host: A}));
    stands[A].fail = (m, p) => m === 'DELETE' && p.startsWith('/users/keys/');
    expect(await share(args('rm', ID.slice(0, 8)))).toBe(false);
    expect(printed()).toContain('could not revoke the folder key');
    expect(stands[A].shares.size).toBe(1);

    stands[A].fail = (m, p) => m === 'DELETE' && p.startsWith('/public/v1/connections/');
    expect(await share(args('rm', ID.slice(0, 8), {host: A}))).toBe(false);
    expect(printed()).toContain('could not delete the connection');
    expect(stands[A].keys).toHaveLength(0);
  });
});

describe('the folder-limited agent key', () => {
  it('is issued by the share, never as a key of the user', async () => {
    expect(await share(args(folder, {host: A}))).toBe(true);
    expect(events).toContain(`POST /local_shares/${ID}/agent_key`);
    expect(events).not.toContain('POST /users/current/keys');
    expect(stands[A].keys.map((k) => k.scope)).toEqual([`connected-folder:${ID}`]);
  });

  it('keeps a running agent whose key is limited to the folder, with no new key', async () => {
    await share(args(folder, {host: A}));
    const keys = stands[A].keys.map((k) => k.id);
    docker.calls = [];
    events.length = 0;
    expect(await share(args(folder, {host: A}))).toBe(true);
    expect(ran()).toHaveLength(0);
    expect(events).not.toContain(`POST /local_shares/${ID}/agent_key`);
    expect(stands[A].keys.map((k) => k.id)).toEqual(keys);
  });

  it('list flags a container whose key is not limited to the folder', async () => {
    await share(args(folder, {host: A}));
    const listed = async () => {
      (console.log as any).mockClear();
      expect(await share(args('list', {host: A, output: 'json'}))).toBe(true);
      return (console.log as any).mock.calls.map((c: any[]) => c.join(' ')).join('\n');
    };
    expect(await listed()).not.toContain('outdated key');
    delete stands[A].keys[0].scope;
    expect(await listed()).toContain('outdated key: re-run connect-folder');
  });
});
