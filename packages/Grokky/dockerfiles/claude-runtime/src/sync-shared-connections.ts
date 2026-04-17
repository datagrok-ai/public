import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {request, downloadFile} from './shared-api-client';
import type {ConnectionInfo} from './sync-utils';
import {connectorPath, collectFiles} from './sync-utils';

const SHARED_AGENTS_CONN_NAME = 'MyFilesAgents';

// Per-user shared files cache: "connId/filePath" → updatedOn timestamp.
const sharedCache = new Map<string, Map<string, string>>();

// Discovers shared "Home/agents" connections and syncs their files.
// Skips the user's own Home connection (already synced separately).
// Files are stored as agents/shared-{connId}-{flatPath} to avoid collisions.
export async function syncSharedConnections(userDir: string, homeConnId: string | null): Promise<void> {
  console.log('shared-connections: searching for shared agents folders...');
  const connections = await request<ConnectionInfo[]>(
    'GET', '/public/v1/connections?text=agents',
  );
  console.log('shared-connections: API returned', (connections ?? []).map((c) => ({id: c.id, name: c.name})));
  const shared = (connections ?? []).filter(
    (c) => c.name === SHARED_AGENTS_CONN_NAME && c.id !== homeConnId,
  );
  if (!shared.length) {
    console.log('shared-connections: no shared agents folders found');
    return;
  }
  console.log(`shared-connections: found ${shared.length} shared agents folder(s)`);

  const userId = path.basename(userDir);
  if (!sharedCache.has(userId))
    sharedCache.set(userId, new Map());
  const cached = sharedCache.get(userId)!;

  const agentsDir = path.join(userDir, 'agents');
  const activeKeys = new Set<string>();

  for (const conn of shared) {
    const prefix = `shared-${conn.id}-`;
    console.log(`shared-connections: syncing ${conn.id}`);

    // The shared connection root IS the agents folder, so list from root
    const remoteFiles = await collectFiles(conn.id, '');
    if (!remoteFiles.length) {
      console.log(`shared-connections: no files in ${conn.id}`);
      continue;
    }

    const cPath = connectorPath(conn);
    for (const remote of remoteFiles) {
      const cacheKey = `${conn.id}/${remote.path}`;
      activeKeys.add(cacheKey);
      if (cached.get(cacheKey) === remote.updatedOn)
        continue;

      const label = cached.has(cacheKey) ? 'updated' : 'new';
      console.log(`shared-connections: downloading ${label} file ${cacheKey}`);
      try {
        const content = await downloadFile(cPath, remote.path) as string;
        const fileName = remote.path.replace(/\//g, '-');
        const dest = path.join(agentsDir, `${prefix}${fileName}`);
        await fs.writeFile(dest, content);
        cached.set(cacheKey, remote.updatedOn);
      } catch (e: any) {
        console.warn(`shared-connections: failed to download ${cacheKey}:`, e.message);
      }
    }
  }

  // Clean up files from connections that are no longer shared
  for (const [key] of cached) {
    if (!activeKeys.has(key)) {
      const slashIdx = key.indexOf('/');
      const connId = key.substring(0, slashIdx);
      const filePath = key.substring(slashIdx + 1).replace(/\//g, '-');
      const localFile = path.join(agentsDir, `shared-${connId}-${filePath}`);
      try {
        await fs.rm(localFile, {force: true});
        console.log(`shared-connections: removed stale file ${localFile}`);
      } catch { /* ignore */ }
      cached.delete(key);
    }
  }
}
