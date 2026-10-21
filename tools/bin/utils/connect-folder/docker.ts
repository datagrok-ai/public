import childProcess from 'child_process';
import crypto from 'crypto';
import fs from 'fs';
import os from 'os';
import path from 'path';

/** The `docker` CLI for `grok connect-folder`: always an argument vector, never a shell, so no quoting or MSYS path mangling. */

export const LABEL = 'ai.datagrok.folder';
export const SHARE_ROOT = '/share';
const DOCKER_HOST = 'host.docker.internal';

export class DockerError extends Error {
  constructor(message: string, public readonly stderr: string = '') {
    super(message);
    this.name = 'DockerError';
  }
}

/** Where `docker` is started; the tests replace it, since a module mock cannot reach the transpiled require. */
export const proc = {spawn: childProcess.spawn};

export interface DockerResult { code: number, stdout: string, stderr: string }

export interface ShareContainer {
  id: string,
  name: string,
  /** `running`, `exited`, `restarting`, `created`, ... */
  state: string,
  image: string,
  labels: Record<string, string>,
  mounts: {Source: string, Destination: string, RW: boolean}[],
}

/** Runs `docker <args>`; [inherit] attaches the terminal (logs -f, a foreground share). */
export function docker(args: string[], inherit: boolean = false): Promise<DockerResult> {
  return new Promise((resolve, reject) => {
    const child = proc.spawn('docker', args, {stdio: inherit ? 'inherit' : ['ignore', 'pipe', 'pipe'], windowsHide: true});
    let stdout = '';
    let stderr = '';
    child.stdout?.on('data', (d) => stdout += d);
    child.stderr?.on('data', (d) => stderr += d);
    child.on('error', (e: any) => reject(e?.code === 'ENOENT'
      ? new DockerError('grok connect-folder runs the folder agent in Docker, which is not installed. ' +
        'Install Docker Desktop (Windows/macOS) or Docker Engine (Linux).')
      : e));
    child.on('close', (code) => resolve({code: code ?? 1, stdout, stderr}));
  });
}

/** [docker], failing with the CLI's own message on a non-zero exit. */
export async function dockerOk(args: string[]): Promise<string> {
  const r = await docker(args);
  if (r.code !== 0)
    throw new DockerError(r.stderr.trim() || `docker ${args[0]} failed with exit code ${r.code}`, r.stderr);
  return r.stdout;
}

export async function dockerAvailable(): Promise<void> {
  const version = await docker(['version', '--format', '{{.Server.Version}}']);
  if (version.code !== 0)
    throw new DockerError('Docker is installed but not running: start Docker Desktop / the docker service.', version.stderr);
  const osType = (await dockerOk(['info', '--format', '{{.OSType}}'])).trim();
  if (osType !== 'linux') {
    throw new DockerError(`Docker runs ${osType} containers; the folder agent needs Linux containers. ` +
      'Switch Docker Desktop to Linux containers.');
  }
}

/**
 * The absolute path Docker mounts for [folder]. Refuses network and `\\wsl$` paths, which Docker Desktop
 * cannot bind; [note] is set for a Windows drive seen from WSL, which works but is slow.
 */
export function hostFolder(folder: string, platform: string = process.platform,
  resolve: (p: string) => string = (p) => fs.realpathSync.native(p)): {path: string, note?: string} {
  const raw = folder.trim();
  if (/^[\\/]{2}(wsl\$|wsl\.localhost)([\\/]|$)/i.test(raw))
    throw new Error(`${raw} is inside WSL: run grok connect-folder inside WSL for WSL paths`);
  if (/^[\\/]{2}[^\\/]/.test(raw))
    throw new Error(`${raw} is a network path; Docker cannot mount it. Connect a local folder instead`);
  const p = platform === 'win32' ? path.win32 : path.posix;
  const resolved = resolve(p.resolve(raw));
  if (/^[\\/]{2}/.test(resolved))
    throw new Error(`${raw} resolves to the network path ${resolved}; Docker cannot mount it`);
  const note = platform === 'linux' && /^\/mnt\/[a-z](\/|$)/i.test(resolved)
    ? 'Windows drives seen from WSL (/mnt/...) are slow to read through Docker; consider a folder in the WSL filesystem'
    : undefined;
  return {path: resolved, note};
}

/**
 * The folder the server derives the share id from: `/` separators, no trailing slash, lower case where paths are
 * case-insensitive, and a Windows drive seen from WSL (`/mnt/c/x`) the same as on Windows (`c:/x`).
 */
export function canonicalFolder(abs: string, platform: string = process.platform): string {
  let f = abs.replace(/\\/g, '/');
  const wslDrive = platform === 'linux' ? /^\/mnt\/([a-z])(\/.*)?$/i.exec(f) : null;
  if (wslDrive)
    f = `${wslDrive[1]}:${wslDrive[2] ?? '/'}`;
  while (f.length > 1 && f.endsWith('/') && !f.endsWith(':/'))
    f = f.slice(0, -1);
  return platform === 'win32' || platform === 'darwin' || wslDrive ? f.toLowerCase() : f;
}

/** `Mac.local` and `Mac` are one machine: DHCP and macOS flip between them. */
export function canonicalHost(hostname: string = os.hostname()): string {
  return hostname.trim().toLowerCase().split('.')[0];
}

/** `localhost` in a container is the container itself; Docker reaches the host as `host.docker.internal`. */
export function containerUrl(url: string): {url: string, rewritten: boolean} {
  const u = new URL(url);
  if (!['localhost', '127.0.0.1', '[::1]'].includes(u.hostname))
    return {url, rewritten: false};
  u.hostname = DOCKER_HOST;
  return {url: u.href.replace(/\/$/, ''), rewritten: true};
}

/** The same login shares a folder under the same id on every server, so the name hashes the server in. */
export function containerName(apiRoot: string, connectionId: string, foreground: boolean = false): string {
  const hash = crypto.createHash('sha256').update(`${apiRoot}|${connectionId}`).digest('hex').slice(0, 12);
  return `grok-folder-${hash}${foreground ? '-fg' : ''}`;
}

/** A `--mount` field: the value is CSV, so a path with a comma or a quote must be quoted. */
function mountField(field: string): string {
  return /[",]/.test(field) ? `"${field.replace(/"/g, '""')}"` : field;
}

export interface RunSpec {
  image: string,
  connectionId: string,
  alias: string,
  /** The API root as the user configured it (labels, list). */
  apiRoot: string,
  name: string,
  /** The absolute host folder, as Docker mounts it. */
  folder: string,
  hostLabel: string,
  readOnly: boolean,
  /** The 0600 file holding `GROK_SHARE_KEY`. */
  envFile: string,
  foreground?: boolean,
  platform?: string,
  /** `uid:gid` on Linux, so files the agent writes belong to the user. */
  user?: string,
  /** The server id of the share's own key, whose private half is in [envFile]. */
  keyId?: string,
}

/** Whether [c] runs [spec] as is: the mount, the mode, the image and the user are fixed at `docker run`. */
export function sameSettings(c: ShareContainer, spec: RunSpec): boolean {
  const l = c.labels;
  return c.image === spec.image && l[`${LABEL}.folder`] === spec.folder && l[`${LABEL}.readonly`] === `${spec.readOnly}` &&
    (l[`${LABEL}.user`] ?? '') === (runUser(spec) ?? '');
}

function runUser(spec: RunSpec): string | undefined {
  return (spec.platform ?? process.platform) === 'linux' ? spec.user : undefined;
}

export function runArgs(spec: RunSpec): string[] {
  const platform = spec.platform ?? process.platform;
  const server = containerUrl(spec.apiRoot);
  const labels: Record<string, string> = {
    [LABEL]: '1',
    [`${LABEL}.connection`]: spec.connectionId,
    [`${LABEL}.alias`]: spec.alias,
    [`${LABEL}.url`]: spec.apiRoot,
    [`${LABEL}.name`]: spec.name,
    [`${LABEL}.folder`]: spec.folder,
    [`${LABEL}.readonly`]: `${spec.readOnly}`,
    [`${LABEL}.user`]: runUser(spec) ?? '',
    [`${LABEL}.key`]: spec.keyId ?? '',
  };
  const args = ['run', '--name', containerName(spec.apiRoot, spec.connectionId, spec.foreground)];
  args.push(...(spec.foreground ? ['--rm'] : ['-d', '--restart', 'unless-stopped']));
  for (const [k, v] of Object.entries(labels))
    args.push('--label', `${k}=${v}`);
  args.push('--mount', ['type=bind', mountField(`source=${spec.folder}`), `target=${SHARE_ROOT}`,
    ...(spec.readOnly ? ['readonly'] : [])].join(','));
  args.push('--env-file', spec.envFile);
  args.push('--cap-drop', 'ALL', '--security-opt', 'no-new-privileges', '--memory', '256m', '-e', 'DEBUG_MODE=false');
  if (runUser(spec))
    args.push('--user', runUser(spec)!);
  // Docker Desktop resolves host.docker.internal by itself; Docker Engine on Linux (and WSL) needs it mapped to the host.
  if (platform === 'linux' && server.rewritten)
    args.push('--add-host', `${DOCKER_HOST}:host-gateway`);
  args.push(spec.image, '--tunnel', '--server', server.url, '--root', SHARE_ROOT, '--connection', spec.connectionId,
    '--host-label', spec.hostLabel, '--folder-label', spec.folder);
  if (spec.readOnly)
    args.push('--read-only');
  return args;
}

export function linuxUser(platform: string = process.platform): string | undefined {
  return platform === 'linux' && process.getuid && process.getgid ? `${process.getuid()}:${process.getgid()}` : undefined;
}

/** Share containers with all of [labels] (`label key → value`), stopped ones included. */
export async function ps(labels: Record<string, string> = {}): Promise<ShareContainer[]> {
  const filters = [`label=${LABEL}=1`, ...Object.entries(labels).map(([k, v]) => `label=${k}=${v}`)];
  const ids = (await dockerOk(['ps', '-a', '-q', ...filters.flatMap((f) => ['--filter', f])])).split(/\s+/).filter((s) => s);
  if (!ids.length)
    return [];
  // `docker ps --format` joins labels with commas, which a folder name can contain; inspect keeps them apart.
  return (JSON.parse(await dockerOk(['inspect', ...ids])) as any[]).map((c) => ({
    id: c.Id,
    name: String(c.Name ?? '').replace(/^\//, ''),
    state: c.State?.Status ?? '',
    image: c.Config?.Image ?? '',
    labels: c.Config?.Labels ?? {},
    mounts: c.Mounts ?? [],
  }));
}

/** Pulls [image] unless Docker already has it (a local build has no registry to pull from). */
export async function ensureImage(image: string): Promise<void> {
  if ((await docker(['image', 'inspect', '--format', '{{.Id}}', image])).code !== 0)
    await dockerOk(['pull', image]);
}

/** Whether [container] is running, waiting up to [seconds] for it to get there. */
export async function isRunning(container: string, seconds: number = 0): Promise<boolean> {
  for (let i = 0; ; i++) {
    const r = await docker(['inspect', '--format', '{{.State.Running}}', container]);
    if (r.code === 0 && r.stdout.trim() === 'true')
      return true;
    if (i >= seconds * 4)
      return false;
    await new Promise((resolve) => setTimeout(resolve, 250));
  }
}

export const start = (container: string) => dockerOk(['start', container]);
export const stop = (container: string) => dockerOk(['stop', container]);
export const rm = (container: string) => dockerOk(['rm', '-f', container]);
export const logs = (container: string, follow: boolean = false) =>
  docker(['logs', ...(follow ? ['-f'] : []), container], true);
