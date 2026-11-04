import {describe, it, expect, beforeEach} from 'vitest';
import {EventEmitter} from 'events';
import * as dk from '../utils/connect-folder/docker';
import {checkFolderSafety, matches} from '../commands/connect-folder';

const spawned = {calls: [] as string[][], next: [] as ((child: any) => void)[]};

dk.proc.spawn = ((_cmd: string, args: string[]) => {
  spawned.calls.push(args);
  const child: any = new EventEmitter();
  child.stdout = new EventEmitter();
  child.stderr = new EventEmitter();
  const behave = spawned.next.shift() ?? ((c: any) => c.emit('close', 0));
  setImmediate(() => behave(child));
  return child;
}) as any;

beforeEach(() => { spawned.calls = []; spawned.next = []; });

const ID = '7ff0714e-1d2c-5a3b-9e8f-0123456789ab';

const spec = (over: Partial<dk.RunSpec> = {}): dk.RunSpec => ({
  image: 'datagrok/proxy_storage:local', connectionId: ID, alias: 'dev', apiRoot: 'https://dev.datagrok.ai/api',
  name: 'laptop data', folder: '/home/a/data', hostLabel: 'laptop', readOnly: false, envFile: '/tmp/k/key.env',
  platform: 'linux', user: '1000:1000', keyId: 'key-1', ...over,
});

const after = (args: string[], name: string) => args[args.indexOf(name) + 1];

describe('runArgs', () => {
  it('runs a detached read-write share with one bind mount, the key file and the hardening flags', () => {
    expect(dk.runArgs(spec())).toEqual([
      'run', '--name', dk.containerName('https://dev.datagrok.ai/api', ID), '-d', '--restart', 'unless-stopped',
      '--label', 'ai.datagrok.folder=1',
      '--label', `ai.datagrok.folder.connection=${ID}`,
      '--label', 'ai.datagrok.folder.alias=dev',
      '--label', 'ai.datagrok.folder.url=https://dev.datagrok.ai/api',
      '--label', 'ai.datagrok.folder.name=laptop data',
      '--label', 'ai.datagrok.folder.folder=/home/a/data',
      '--label', 'ai.datagrok.folder.readonly=false',
      '--label', 'ai.datagrok.folder.user=1000:1000',
      '--label', 'ai.datagrok.folder.key=key-1',
      '--mount', 'type=bind,source=/home/a/data,target=/share',
      '--env-file', '/tmp/k/key.env',
      '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges', '--memory', '256m', '-e', 'DEBUG_MODE=false',
      '--user', '1000:1000',
      'datagrok/proxy_storage:local', '--tunnel', '--server', 'https://dev.datagrok.ai/api', '--root', '/share',
      '--connection', ID, '--host-label', 'laptop', '--folder-label', '/home/a/data',
    ]);
  });

  it('mounts read-only and tells the agent so', () => {
    const args = dk.runArgs(spec({readOnly: true}));
    expect(after(args, '--mount')).toBe('type=bind,source=/home/a/data,target=/share,readonly');
    expect(args[args.length - 1]).toBe('--read-only');
    expect(args).toContain('ai.datagrok.folder.readonly=true');
  });

  it('passes a Windows path as is, without --user', () => {
    const args = dk.runArgs(spec({platform: 'win32', folder: 'C:\\Users\\a\\my data', user: undefined}));
    expect(after(args, '--mount')).toBe('type=bind,source=C:\\Users\\a\\my data,target=/share');
    expect(after(args, '--folder-label')).toBe('C:\\Users\\a\\my data');
    expect(args).not.toContain('--user');
    expect(args).not.toContain('--add-host');
  });

  it('never passes --user off Linux, even when given', () => {
    expect(dk.runArgs(spec({platform: 'darwin'}))).not.toContain('--user');
  });

  it('quotes a mount source holding a comma', () => {
    expect(after(dk.runArgs(spec({folder: '/home/a/x,y'})), '--mount')).toBe('type=bind,"source=/home/a/x,y",target=/share');
  });

  it('passes a WSL path as is', () => {
    expect(after(dk.runArgs(spec({folder: '/mnt/c/Users/a'})), '--mount')).toBe('type=bind,source=/mnt/c/Users/a,target=/share');
  });

  it('points the agent at host.docker.internal for a localhost server, and maps it on Linux', () => {
    const args = dk.runArgs(spec({apiRoot: 'http://localhost:55287'}));
    expect(after(args, '--server')).toBe('http://host.docker.internal:55287');
    expect(after(args, '--add-host')).toBe('host.docker.internal:host-gateway');
    expect(args).toContain('ai.datagrok.folder.url=http://localhost:55287');
    expect(dk.runArgs(spec({apiRoot: 'http://127.0.0.1:8080/api', platform: 'win32'}))).not.toContain('--add-host');
  });

  it('runs a foreground share removed on exit, without a restart policy', () => {
    const args = dk.runArgs(spec({foreground: true}));
    expect(after(args, '--name')).toBe(`${dk.containerName('https://dev.datagrok.ai/api', ID)}-fg`);
    expect(args).toContain('--rm');
    expect(args).not.toContain('--restart');
    expect(args).not.toContain('-d');
  });
});

describe('containerUrl / containerName', () => {
  it('leaves a real host alone', () => {
    expect(dk.containerUrl('https://dev.datagrok.ai/api')).toEqual({url: 'https://dev.datagrok.ai/api', rewritten: false});
    expect(dk.containerUrl('http://[::1]:8082/api').url).toBe('http://host.docker.internal:8082/api');
  });

  it('derives a valid name from the server and the id, so one id on two servers gets two', () => {
    expect(dk.containerName('https://a/api', ID)).toMatch(/^grok-folder-[0-9a-f]{12}$/);
    expect(dk.containerName('https://a/api', ID)).toBe(dk.containerName('https://a/api', ID));
    expect(dk.containerName('https://a/api', ID)).not.toBe(dk.containerName('https://b/api', ID));
  });
});

describe('hostFolder', () => {
  const same = (p: string) => p;

  it('refuses UNC paths', () => {
    expect(() => dk.hostFolder('\\\\server\\share\\data', 'win32', same)).toThrow('network path');
    expect(() => dk.hostFolder('//server/share', 'win32', same)).toThrow('network path');
  });

  it('sends WSL paths to WSL', () => {
    expect(() => dk.hostFolder('\\\\wsl$\\Ubuntu\\home\\a', 'win32', same)).toThrow('run grok connect-folder inside WSL');
    expect(() => dk.hostFolder('\\\\wsl.localhost\\Ubuntu\\home', 'win32', same)).toThrow('run grok connect-folder inside WSL');
  });

  it('refuses a folder that resolves to a network path', () => {
    expect(() => dk.hostFolder('Z:\\data', 'win32', () => '\\\\nas\\data')).toThrow('network path');
  });

  it('resolves a Windows path and notes the slow WSL view of a Windows drive', () => {
    expect(dk.hostFolder('C:\\Users\\a\\data\\', 'win32', same)).toEqual({path: 'C:\\Users\\a\\data', note: undefined});
    expect(dk.hostFolder('/mnt/c/Users/a', 'linux', same).note).toContain('slow');
    expect(dk.hostFolder('/home/a', 'linux', same).note).toBeUndefined();
  });
});

describe('canonicalFolder / canonicalHost', () => {
  it('matches the server-side id inputs', () => {
    expect(dk.canonicalFolder('C:\\Users\\A\\Data\\', 'win32')).toBe('c:/users/a/data');
    expect(dk.canonicalFolder('C:\\', 'win32')).toBe('c:/');
    expect(dk.canonicalFolder('/mnt/c/Users/A/Data', 'linux')).toBe('c:/users/a/data');
    expect(dk.canonicalFolder('/home/A/Data/', 'linux')).toBe('/home/A/Data');
    expect(dk.canonicalFolder('/Users/A/Data', 'darwin')).toBe('/users/a/data');
    expect(dk.canonicalFolder('/', 'linux')).toBe('/');
  });

  it('drops the domain and the case of the host', () => {
    expect(dk.canonicalHost('Mac.local')).toBe('mac');
    expect(dk.canonicalHost('LAPTOP')).toBe('laptop');
  });
});

describe('docker', () => {
  it('says Docker is not installed when the binary is missing', async () => {
    spawned.next.push((c) => c.emit('error', Object.assign(new Error('spawn docker ENOENT'), {code: 'ENOENT'})));
    await expect(dk.dockerAvailable()).rejects.toThrow('which is not installed');
  });

  it('says Docker is not running when the daemon does not answer', async () => {
    spawned.next.push((c) => { c.stderr.emit('data', 'Cannot connect to the Docker daemon'); c.emit('close', 1); });
    await expect(dk.dockerAvailable()).rejects.toThrow('installed but not running');
  });

  it('refuses Windows containers', async () => {
    spawned.next.push((c) => { c.stdout.emit('data', '27.1.1'); c.emit('close', 0); });
    spawned.next.push((c) => { c.stdout.emit('data', 'windows\n'); c.emit('close', 0); });
    await expect(dk.dockerAvailable()).rejects.toThrow('needs Linux containers');
  });

  it('lists share containers by label, with exact labels from inspect', async () => {
    spawned.next.push((c) => { c.stdout.emit('data', 'abc\n'); c.emit('close', 0); });
    spawned.next.push((c) => {
      c.stdout.emit('data', JSON.stringify([{Id: 'abc', Name: '/grok-folder-7ff0714e1d2c', State: {Status: 'running'},
        Image: 'sha256:1', Config: {Image: 'img', Labels: {'ai.datagrok.folder.folder': '/home/a/x,y'}}, Mounts: []}]));
      c.emit('close', 0);
    });
    const found = await dk.ps({'ai.datagrok.folder.connection': ID});
    expect(spawned.calls[0]).toEqual(['ps', '-a', '-q', '--filter', 'label=ai.datagrok.folder=1',
      '--filter', `label=ai.datagrok.folder.connection=${ID}`]);
    expect(spawned.calls[1]).toEqual(['inspect', 'abc']);
    expect(found).toEqual([{id: 'abc', name: 'grok-folder-7ff0714e1d2c', state: 'running', image: 'img',
      imageId: 'sha256:1', labels: {'ai.datagrok.folder.folder': '/home/a/x,y'}, mounts: []}]);
  });

  it('pulls a tag every time, and falls back to the local copy when the pull fails', async () => {
    const id = (c: any) => { c.stdout.emit('data', 'sha256:2\n'); c.emit('close', 0); };
    spawned.next.push((c) => c.emit('close', 0), id);
    expect(await dk.ensureImage('datagrok/proxy_storage:bleeding-edge')).toEqual({id: 'sha256:2'});
    expect(spawned.calls.map((c) => c[0])).toEqual(['pull', 'image']);

    spawned.calls = [];
    spawned.next.push((c) => { c.stderr.emit('data', 'dial tcp: lookup registry-1.docker.io\n'); c.emit('close', 1); }, id);
    expect(await dk.ensureImage('datagrok/proxy_storage:local'))
      .toEqual({id: 'sha256:2', pullError: 'dial tcp: lookup registry-1.docker.io'});

    spawned.next.push((c) => { c.stderr.emit('data', 'pull access denied\n'); c.emit('close', 1); },
      (c) => c.emit('close', 1));
    await expect(dk.ensureImage('datagrok/proxy_storage:gone')).rejects.toThrow('pull access denied');
  });

  it('pulls a digest only when Docker does not have it', async () => {
    const digest = 'datagrok/proxy_storage@sha256:abc';
    spawned.next.push((c) => { c.stdout.emit('data', 'sha256:3\n'); c.emit('close', 0); });
    expect(await dk.ensureImage(digest)).toEqual({id: 'sha256:3'});
    expect(spawned.calls.map((c) => c[0])).toEqual(['image']);

    spawned.calls = [];
    spawned.next.push((c) => c.emit('close', 1), (c) => c.emit('close', 0),
      (c) => { c.stdout.emit('data', 'sha256:3\n'); c.emit('close', 0); });
    expect(await dk.ensureImage(digest)).toEqual({id: 'sha256:3'});
    expect(spawned.calls.map((c) => c[0])).toEqual(['image', 'pull', 'image']);
  });

  it('removes with -f and surfaces the CLI error', async () => {
    spawned.next.push((c) => { c.stderr.emit('data', 'Error: No such container: x\n'); c.emit('close', 1); });
    await expect(dk.rm('x')).rejects.toThrow('No such container');
    expect(spawned.calls[0]).toEqual(['rm', '-f', 'x']);
  });
});

describe('checkFolderSafety', () => {
  const home = '/home/a';

  it('refuses ~/.grok, folders inside it, and folders holding it, even with --force', () => {
    for (const f of ['/home/a/.grok', '/home/a/.grok/keys', '/home/a', '/'])
      expect(() => checkFolderSafety(f, true, home, 'linux')).toThrow('private keys');
  });

  it('wants --force for another root', () => {
    expect(() => checkFolderSafety('D:\\', false, 'C:\\Users\\a', 'win32')).toThrow('--force');
    expect(() => checkFolderSafety('D:\\', true, 'C:\\Users\\a', 'win32')).not.toThrow();
  });

  it('compares Windows paths without case', () => {
    expect(() => checkFolderSafety('c:\\users\\A\\.GROK', true, 'C:\\Users\\a', 'win32')).toThrow('private keys');
    expect(() => checkFolderSafety('C:\\Users\\a\\data', false, 'C:\\Users\\a', 'win32')).not.toThrow();
  });

  it('does not take a sibling with a common prefix for ~/.grok', () => {
    expect(() => checkFolderSafety('/home/a/.grokker', false, home, 'linux')).not.toThrow();
  });
});

describe('matches', () => {
  it('takes a name, an id prefix of 8 or more, or the folder', () => {
    expect(matches('Laptop Data', ID, 'laptop data', '/nope')).toBe(true);
    expect(matches('7ff0714e', ID, 'x', undefined)).toBe(true);
    expect(matches('7ff07', ID, 'x', undefined)).toBe(false);
    expect(matches(process.cwd(), ID, 'x', process.cwd())).toBe(true);
  });
});
