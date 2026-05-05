import * as fs from 'node:fs/promises';
import * as path from 'node:path';

const USERS_DIR = '/users';
const known = new Set<string>();
const inFlight = new Map<string, Promise<string>>();

export function userIdFromKey(apiKey: string): string {
  try {
    const token = apiKey.startsWith('Bearer ') ? apiKey.slice(7) : apiKey;
    const payload = JSON.parse(Buffer.from(token.split('.')[1], 'base64url').toString());
    const id = payload.sub ?? payload.usr?.id;
    if (!id) {
      console.warn('user-dir: JWT has no sub/usr.id, falling back to "default"');
      return 'default';
    }
    return id;
  } catch {
    return apiKey.length > 64 ? 'default' : apiKey;
  }
}

export function getUserDir(apiKey: string): string {
  return path.join(USERS_DIR, userIdFromKey(apiKey));
}

export function userDirFromId(userId: string): string {
  return path.join(USERS_DIR, userId);
}

export function ensureUserDir(apiKey: string): Promise<string> {
  const dir = getUserDir(apiKey);
  if (known.has(dir)) return Promise.resolve(dir);
  const existing = inFlight.get(dir);
  if (existing) return existing;
  const p = (async () => {
    await fs.mkdir(path.join(dir, 'agents'), {recursive: true});
    await fs.mkdir(path.join(dir, 'workspace'), {recursive: true});
    known.add(dir);
    console.log(`user-dir: ensured layout at ${dir}`);
    return dir;
  })().finally(() => inFlight.delete(dir));
  inFlight.set(dir, p);
  return p;
}
