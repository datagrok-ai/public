import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {runWithContext, request} from './shared-api-client';
import {resolveHomeConnection, syncHomeFiles} from './sync-home-files';
import {syncPackages} from './sync-packages';
import {syncSharedConnections} from './sync-shared-connections';
import {loadPackageKnowledge} from './package-knowledge-tool';
import {getInstalledPackages} from './installed-packages';
import {ensureUserDir, getUserDir} from './user-dir';
import {buildStagedWorkspace} from './staged-workspace';

export {ensureUserDir, getUserDir, userIdFromKey} from './user-dir';

export const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

// ── Shared types ───────────────────────────────────────────────────────

export interface EntityTag {
  id?: string;
  tag: string;
  type?: string | null;
}

export interface ConnectionInfo {
  id: string;
  name: string;
  namespace?: string;
  entityTags?: EntityTag[];
}

export interface RemoteFileEntry {
  path: string;       // relative path from the root dir, e.g. "subdir/file.md"
  updatedOn: string;  // ISO timestamp from the server
}

// Response shape from the Datagrok connectors file-listing API (FileInfo in Dart).
interface FileEntry {
  name?: string;
  friendlyName?: string;
  isFile?: boolean;
  updatedOn?: string;
}

// ── Shared utilities ───────────────────────────────────────────────────

export function connectorPath(conn: ConnectionInfo): string {
  const ns = conn.namespace ?? 'System:';
  return ns.endsWith(':') ? `${ns}${conn.name}` : `${ns}:${conn.name}`;
}

export async function collectFiles(connId: string, dirPath: string): Promise<RemoteFileEntry[]> {
  let entries: FileEntry[];
  try {
    const filePath = dirPath ? `${dirPath.replace(/^\//, '')}/` : '';
    entries = await request<FileEntry[]>('GET', `/connectors/connections/${connId}/files/${filePath}`);
  } catch (e: any) {
    console.warn(`user-files: failed to list ${dirPath}:`, e.message);
    return [];
  }
  
  if (!Array.isArray(entries))
    return [];

  const files: RemoteFileEntry[] = [];
  for (const entry of entries) {
    const name = entry.name ?? entry.friendlyName ?? '';
    if (!name)
      continue;
    const rel = dirPath ? `${dirPath}/${name}` : name;
    if (!entry.isFile) {
      console.log(`user-files: descending into directory ${rel}`);
      files.push(...await collectFiles(connId, rel));
    } else {
      files.push({path: rel, updatedOn: entry.updatedOn ?? ''});
    }
  }
  return files;
}

// ── List agent files on disk ───────────────────────────────────────────

async function listAgentFiles(userDir: string): Promise<string[]> {
  const agentsDir = path.join(userDir, 'agents');
  try {
    const entries = await fs.readdir(agentsDir, {recursive: true, withFileTypes: true});
    const files: string[] = [];
    for (const entry of entries) {
      if (!entry.isFile())
        continue;
      files.push(path.relative(agentsDir, path.join(entry.parentPath, entry.name)));
    }
    console.log(`user-files: listed ${files.length} agent file(s) on disk`);
    return files;
  } catch {
    return [];
  }
}

// ── Sync orchestration ─────────────────────────────────────────────────

const syncInFlight = new Map<string, Promise<{dir: string; files: string[]}>>();
const initialSyncDone = new Set<string>();

export type SyncScope = 'all' | 'user-files' | 'packages' | 'shared';

export async function syncUserFiles(
  apiUrl: string, apiKey: string, scope: SyncScope = 'all', packageName?: string,
): Promise<{dir: string; files: string[]}> {
  const dir = getUserDir(apiKey);

  // Fast path: already synced and scope is 'all' (user message) → just return current state
  if (initialSyncDone.has(apiKey) && scope === 'all') {
    console.log('user-files: already synced, returning cached state');
    const files = await listAgentFiles(dir);
    return {dir, files};
  }

  // If a sync is already running for this key, wait for it
  const existing = syncInFlight.get(apiKey);
  if (existing && scope === 'all') {
    console.log('user-files: sync already in flight, awaiting...');
    return existing;
  }

  const promise = doSync(apiUrl, apiKey, scope, packageName);
  syncInFlight.set(apiKey, promise);
  try {
    return await promise;
  }
  finally {
    if (syncInFlight.get(apiKey) === promise)
      syncInFlight.delete(apiKey);
  }
}

async function doSync(
  apiUrl: string, apiKey: string, scope: SyncScope, packageName?: string,
): Promise<{dir: string; files: string[]}> {
  const dir = getUserDir(apiKey);
  const userId = path.basename(dir);
  console.log(`user-files: starting sync for user ${userId} (scope=${scope})`);

  await ensureUserDir(apiKey);

  let homeConnId: string | null = null;
  if (scope === 'all' || scope === 'user-files') {
    await runWithContext({apiKey, apiUrl}, async () => {
      const home = await resolveHomeConnection();
      if (!home) {
        console.log('user-files: skipping user file sync (no Home connection)');
        return;
      }
      homeConnId = home.connId;
      const downloaded = await syncHomeFiles(dir, home.connId, home.connectorPath);
      if (downloaded.length)
        console.log(`user-files: synced ${downloaded.length} user file(s)`);
      else
        console.log('user-files: all user files up-to-date');
    });
  }

  if (scope === 'all' || scope === 'shared') {
    await runWithContext({apiKey, apiUrl}, async () => {
      try {
        await syncSharedConnections(dir, homeConnId);
      } catch (e: any) {
        console.warn('shared-connections: sync failed:', e.message);
      }
    });
  }

  if (scope === 'all' || scope === 'packages') {
    await runWithContext({apiKey, apiUrl}, async () => {
      try {
        await syncPackages(dir, packageName);
      } catch (e: any) {
        console.warn('package-agents: sync failed:', e.message);
      }
    });

    const installed = getInstalledPackages(userId);
    if (installed) {
      try {
        await buildStagedWorkspace(dir, installed.keys());
      } catch (e: any) {
        console.warn('staged-workspace: build failed:', e.message);
      }
    }
  }

  if (scope === 'all')
    initialSyncDone.add(apiKey);

  const files = await listAgentFiles(dir);
  console.log(`user-files: sync complete — ${files.length} total agent file(s)`);
  return {dir, files};
}

// ── Package index generation ──────────────────────────────────────────

export async function generatePackageIndex(userId?: string): Promise<string | null> {
  const map = await loadPackageKnowledge();
  if (map.size === 0)
    return null;

  const installed = userId ? getInstalledPackages(userId) : undefined;
  const visible = installed
    ? [...map.values()].filter((p) => installed.has(p.packageName))
    : [...map.values()];
  if (!visible.length)
    return null;

  visible.sort((a, b) => a.packageName.localeCompare(b.packageName));

  let md = '| Package | Description | Keywords |\n';
  md += '|---------|-------------|----------|\n';
  for (const pkg of visible)
    md += `| ${pkg.packageName} | ${pkg.description} | ${pkg.keywords.join(', ')} |\n`;

  return md;
}
