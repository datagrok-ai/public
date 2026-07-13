import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {runWithContext, request} from '../shared-api-client';
import {resolveHomeConnection, syncHomeFiles} from './home-files';
import {syncPackages} from './packages';
import {syncSharedConnections} from './shared-connections';
import {getInstalledPackages} from '../user/installed-packages';
import {ensureUserDir, getUserDir} from '../user/user-dir';
import {buildStagedWorkspace} from '../user/staged-workspace';

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
// apiKey → epoch ms of the last full ('all') sync. Drives the on-demand TTL refresh.
const lastFullSync = new Map<string, number>();
const REFRESH_TTL_MS = 15 * 60 * 1000;

export type SyncScope = 'all' | 'user-files' | 'packages' | 'shared';

export async function syncUserFiles(
  apiUrl: string, apiKey: string, scope: SyncScope = 'all', packageName?: string,
): Promise<{dir: string; files: string[]}> {
  const dir = getUserDir(apiKey);

  // scope 'all' = a Claude query is starting (chat, panel, or search — all route
  // through handleMessage). After the first full sync, serve disk immediately and
  // refresh remote sources in the background once the TTL has elapsed, so a turn
  // never pays sync latency. Freshness lags by at most one turn — fine for
  // low-churn shared skills / package agents.
  if (scope === 'all' && lastFullSync.has(apiKey)) {
    const age = Date.now() - lastFullSync.get(apiKey)!;
    if (age > REFRESH_TTL_MS && !syncInFlight.has(apiKey)) {
      console.log('user-files: refresh TTL elapsed — refreshing packages + shared in background');
      const refresh = (async () => {
        // This user's own agent-folder edits are caught by the package.ts file
        // listeners; the recurring refresh only needs remote-actor changes.
        await doSync(apiUrl, apiKey, 'shared');
        await doSync(apiUrl, apiKey, 'packages');
        lastFullSync.set(apiKey, Date.now());
        return {dir, files: await listAgentFiles(dir)};
      })();
      syncInFlight.set(apiKey, refresh);
      refresh
        .catch((e: any) => console.warn('user-files: background refresh failed:', e.message))
        .finally(() => {
          if (syncInFlight.get(apiKey) === refresh)
            syncInFlight.delete(apiKey);
        });
    }
    return {dir, files: await listAgentFiles(dir)};
  }

  // First 'all' sync (blocking, so files exist for this turn) or an explicit
  // partial scope from an event ('user-files' / 'packages' / 'shared').
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
    lastFullSync.set(apiKey, Date.now());

  const files = await listAgentFiles(dir);
  console.log(`user-files: sync complete — ${files.length} total agent file(s)`);
  return {dir, files};
}
