import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {rewriteForDocker} from '../constants';

const USERS_DIR = '/users';
const VERIFY_TTL_MS = 10 * 60 * 1000;
const REJECT_TTL_MS = 30 * 1000;
const ID_RE = /^[a-zA-Z0-9-]{1,64}$/;

export class UserDirectory {
  private apiUrl = process.env.DATAGROK_API_URL && rewriteForDocker(process.env.DATAGROK_API_URL);
  private readonly known = new Set<string>();
  private readonly creating = new Map<string, Promise<string>>();
  private readonly verified = new Map<string, {id: string; at: number}>();
  private readonly verifying = new Map<string, Promise<string>>();
  private readonly rejected = new Map<string, {msg: string; at: number}>();
  private readonly ignoredUrls = new Set<string>();

  async resolveUserId(apiKey: string, apiUrl?: string): Promise<string> {
    this.pinApiUrl(apiUrl);
    if (!this.apiUrl)
      throw new Error('user verification failed: no Datagrok API URL configured');
    const hit = this.verified.get(apiKey);
    if (hit && Date.now() - hit.at < VERIFY_TTL_MS)
      return hit.id;
    const bad = this.rejected.get(apiKey);
    if (bad && Date.now() - bad.at < REJECT_TTL_MS)
      throw new Error(bad.msg);
    const existing = this.verifying.get(apiKey);
    if (existing)
      return existing;
    const p = this.verify(apiKey).finally(() => this.verifying.delete(apiKey));
    this.verifying.set(apiKey, p);
    return p;
  }

  private pinApiUrl(url?: string): void {
    if (!url || this.apiUrl === url)
      return;
    if (this.apiUrl) {
      if (!this.ignoredUrls.has(url)) {
        this.ignoredUrls.add(url);
        console.warn(`user-dir: ignoring API URL change to ${url} (pinned to ${this.apiUrl})`);
      }
    } else {
      this.apiUrl = url;
      console.log(`user-dir: verifying identities against ${url}`);
    }
  }

  private async verify(apiKey: string): Promise<string> {
    const res = await fetch(`${this.apiUrl}/users/current`, {
      headers: {'Authorization': apiKey},
      signal: AbortSignal.timeout(10_000),
    });
    if (!res.ok || res.headers.has('api-error'))
      this.reject(apiKey, `user verification failed: HTTP ${res.status} from /users/current`);
    const id = ((await res.json()) as {id?: unknown}).id;
    if (typeof id !== 'string' || !ID_RE.test(id))
      this.reject(apiKey, 'user verification failed: /users/current returned an invalid id');
    UserDirectory.sweep(this.verified, VERIFY_TTL_MS);
    this.verified.set(apiKey, {id, at: Date.now()});
    return id;
  }

  // Only server rejections are cached — a transient network error must not lock a user out.
  private reject(apiKey: string, msg: string): never {
    UserDirectory.sweep(this.rejected, REJECT_TTL_MS);
    this.rejected.set(apiKey, {msg, at: Date.now()});
    throw new Error(msg);
  }

  private static sweep<V extends {at: number}>(map: Map<string, V>, ttl: number): void {
    const now = Date.now();
    for (const [k, v] of map) {
      if (now - v.at >= ttl)
        map.delete(k);
    }
  }

  dirFromId(userId: string): string {
    return path.join(USERS_DIR, userId);
  }

  // Undefined only while no API URL is known (e.g. direct local dev) — the session then runs
  // with no user dir and the access guard blocks all of /users. A rejected token throws.
  async ensureDir(apiKey: string, apiUrl?: string): Promise<string | undefined> {
    if (!this.apiUrl && !apiUrl)
      return undefined;
    const dir = this.dirFromId(await this.resolveUserId(apiKey, apiUrl));
    if (this.known.has(dir))
      return dir;
    const existing = this.creating.get(dir);
    if (existing)
      return existing;
    const p = (async () => {
      await fs.mkdir(path.join(dir, 'agents'), {recursive: true});
      await fs.mkdir(path.join(dir, 'workspace'), {recursive: true});
      this.known.add(dir);
      console.log(`user-dir: ensured layout at ${dir}`);
      return dir;
    })().finally(() => this.creating.delete(dir));
    this.creating.set(dir, p);
    return p;
  }
}

export const userDirs = new UserDirectory();
