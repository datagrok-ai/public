import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {request, downloadFile} from './shared-api-client';
import type {ConnectionInfo, RemoteFileEntry} from './sync-utils';
import {connectorPath, collectFiles} from './sync-utils';

// Per-user file metadata cache: relative path → updatedOn timestamp.
const userFileCache = new Map<string, Map<string, string>>();

export async function resolveHomeConnection(): Promise<{connId: string; connectorPath: string} | null> {
  console.log('user-files: resolving Home connection...');
  const conns = await request<ConnectionInfo[]>('GET', '/public/v1/connections?text=My%20files');
  const home = conns.find((c) => c.name === 'Home');
  if (!home) {
    console.warn('user-files: Home connection not found among', conns.map((c) => c.name));
    return null;
  }
  console.log(`user-files: resolved Home connection id=${home.id}`);
  return {connId: home.id, connectorPath: connectorPath(home)};
}

// Compares remote file timestamps against a local cache and only
// downloads files that are new or have changed. Removes local files that
// no longer exist on the remote side.
export async function syncHomeFiles(
  userDir: string, connId: string, connPath: string,
): Promise<string[]> {
  const userId = path.basename(userDir);
  if (!userFileCache.has(userId))
    userFileCache.set(userId, new Map());
  const cached = userFileCache.get(userId)!;

  const remoteFiles = await collectFiles(connId, 'agents');
  console.log(`user-files: remote has ${remoteFiles.length} file(s) in agents/`);

  const remoteSet = new Set<string>();
  const downloaded: string[] = [];

  for (const remote of remoteFiles) {
    remoteSet.add(remote.path);
    const cachedTs = cached.get(remote.path);
    console.log(`user-files: checking ${remote.path} — remote ts=${remote.updatedOn || '<empty>'}, cached ts=${cachedTs ?? 'none'}, match=${cachedTs === remote.updatedOn}`);
    if (cachedTs === remote.updatedOn)
      continue;

    const label = cachedTs ? 'updated' : 'new';
    console.log(`user-files: downloading ${label} file ${remote.path}`);
    try {
      const content = await downloadFile(connPath, remote.path) as string;
      const dest = path.join(userDir, remote.path);
      await fs.mkdir(path.dirname(dest), {recursive: true});
      await fs.writeFile(dest, content);
      cached.set(remote.path, remote.updatedOn);
      downloaded.push(remote.path.replace(/^agents\//, ''));
    } catch (e: any) {
      console.warn(`user-files: failed to download ${remote.path}:`, e.message);
    }
  }

  // Remove local files that were deleted on the remote side
  const agentsDir = path.join(userDir, 'agents');
  try {
    const localEntries = await fs.readdir(agentsDir, {recursive: true, withFileTypes: true});
    for (const entry of localEntries) {
      if (!entry.isFile())
        continue;
      const rel = path.relative(agentsDir, path.join(entry.parentPath, entry.name));
      const remotePath = `agents/${rel}`;
      // Only remove files that were synced from My files (not package files).
      if (!remoteSet.has(remotePath) && cached.has(remotePath)) {
        console.log(`user-files: removing deleted file ${remotePath}`);
        await fs.rm(path.join(entry.parentPath, entry.name), {force: true});
        cached.delete(remotePath);
      }
    }
  } catch { /* agents dir may not exist yet */ }

  return downloaded;
}
