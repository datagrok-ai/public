import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import * as crypto from 'node:crypto';
import AdmZip from 'adm-zip';
import {runWithContext, request, requestBinary, downloadFile} from './shared-api-client';

const USERS_DIR = '/users';
const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

// ── Sync state ──────────────────────────────────────────────────────────
// Guards against concurrent syncs for the same apiKey. When a sync is in
// progress the promise is stored here; subsequent callers await it instead
// of starting a second sync.
const syncInFlight = new Map<string, Promise<{dir: string; files: string[]}>>();

// Tracks which users have completed their initial sync. After the first
// sync, subsequent messages just return the cached dir/file list. New
// syncs only happen when triggered by events (file edit, package load)
// or manual sync — those specify a scope to re-sync only the relevant files.
const initialSyncDone = new Set<string>();

// Per-user file metadata cache: relative path → updatedOn timestamp.
// Used for incremental sync — only files whose remote timestamp differs
// from the cached one are re-downloaded.
const userFileCache = new Map<string, Map<string, string>>();

// Per-user package metadata cache: packageName → updatedOn timestamp.
// A package ZIP is only re-downloaded when its updatedOn changes.
const packageTimestamps = new Map<string, Map<string, string>>();

// Per-user shared skills cache: "connName/filePath" → updatedOn timestamp.
const sharedSkillsCache = new Map<string, Map<string, string>>();

// Connections we've already verified have the ai-skills tag.

// ── User ID extraction ──────────────────────────────────────────────────
// Datagrok passes either a JWT (Bearer <token>) or a plain API key.
// We extract a stable user ID from the JWT payload so that each user gets
// their own directory under /users/{userId}/.
function userIdFromKey(apiKey: string): string {
  try {
    const token = apiKey.startsWith('Bearer ') ? apiKey.slice(7) : apiKey;
    const payload = JSON.parse(Buffer.from(token.split('.')[1], 'base64url').toString());
    const id = payload.sub ?? payload.usr?.id;
    if (!id) {
      console.warn('user-files: JWT has no sub/usr.id, falling back to "default"');
      return 'default';
    }
    return id;
  }
  catch {
    // Not a JWT — likely a plain API key. Use it directly if short enough,
    // otherwise fall back to "default".
    console.log('user-files: apiKey is not a JWT, using fallback user ID');
    return apiKey.length > 64 ? 'default' : apiKey;
  }
}

function getUserDir(apiKey: string): string {
  return path.join(USERS_DIR, userIdFromKey(apiKey));
}

// ── Directory setup ─────────────────────────────────────────────────────
// Creates the per-user directory structure:
//   /users/{userId}/
//   /users/{userId}/agents/           ← synced knowledge files
//   /users/{userId}/workspace -> /workspace  ← symlink to the Datagrok repo
async function ensureUserDir(apiKey: string): Promise<string> {
  const dir = getUserDir(apiKey);
  console.log(`user-files: ensuring user dir ${dir}`);
  await fs.mkdir(path.join(dir, 'agents'), {recursive: true});

  // Create or validate the workspace symlink
  const link = path.join(dir, 'workspace');
  try {
    const target = await fs.readlink(link);
    if (target !== WORKSPACE) {
      console.warn(`user-files: symlink ${link} points to ${target} instead of ${WORKSPACE}, recreating`);
      await fs.rm(link);
      await fs.symlink(WORKSPACE, link);
    }
  }
  catch (e: any) {
    if (e.code === 'ENOENT') {
      await fs.symlink(WORKSPACE, link);
      console.log(`user-files: created workspace symlink ${link} -> ${WORKSPACE}`);
    }
    else
      throw e;
  }
  return dir;
}

// ── Path safety ─────────────────────────────────────────────────────────
// Prevents path traversal attacks. A file path from the remote API could
// contain "../" sequences that would escape the target directory.
function safePath(base: string, relative: string): string {
  const dest = path.join(base, relative);
  const resolved = path.resolve(dest);
  const resolvedBase = path.resolve(base);
  if (!resolved.startsWith(resolvedBase + '/') && resolved !== resolvedBase)
    throw new Error(`Path traversal blocked: ${relative}`);
  return dest;
}

// ── My files connection resolution ──────────────────────────────────────
// Every Datagrok user has a "Home" connection (displayed as "My files" in
// the UI). We find it via the connections API and build the connector path
// that the file download API expects (e.g. "System:Home").

interface EntityTag {
  id?: string;
  tag: string;
  type?: string | null;
}

interface ConnectionInfo {
  id: string;
  name: string;
  namespace?: string;
  entityTags?: EntityTag[];
}

async function resolveHomeConnection(): Promise<{connId: string; connectorPath: string} | null> {
  console.log('user-files: resolving Home connection...');
  const conns = await request<ConnectionInfo[]>('GET', '/public/v1/connections?text=My%20files');
  const home = conns.find((c) => c.name === 'Home');
  if (!home) {
    console.warn('user-files: Home connection not found among', conns.map((c) => c.name));
    return null;
  }
  const ns = (home.namespace ?? 'System:').replace(/:$/, '');
  const connectorPath = `${ns}:${home.name}`;
  console.log(`user-files: resolved Home connection id=${home.id}, path=${connectorPath}`);
  return {connId: home.id, connectorPath};
}

// ── Remote file listing ─────────────────────────────────────────────────
// Recursively lists files under a directory in a Datagrok file connection.
// Returns file metadata (relative path + updatedOn timestamp) for
// incremental sync decisions.

interface RemoteFileEntry {
  path: string;       // relative path from the root dir, e.g. "subdir/file.md"
  updatedOn: string;  // ISO timestamp from the server
}

async function collectFiles(connId: string, dirPath: string): Promise<RemoteFileEntry[]> {
  let entries: any[];
  try {
    const filePath = dirPath ? `${dirPath.replace(/^\//, '')}/` : '';
    entries = await request<any[]>('GET', `/connectors/connections/${connId}/files/${filePath}`);
  }
  catch (e: any) {
    console.warn(`user-files: failed to list ${dirPath}:`, e.message);
    return [];
  }
  if (!Array.isArray(entries))
    return [];

  const files: RemoteFileEntry[] = [];
  for (const entry of entries) {
    const name = entry.name ?? entry.friendlyName ?? entry.fileName ?? '';
    if (!name)
      continue;
    const rel = dirPath ? `${dirPath}/${name}` : name;
    if (!entry.isFile) {
      console.log(`user-files: descending into directory ${rel}`);
      files.push(...await collectFiles(connId, rel));
    }
    else {
      files.push({path: rel, updatedOn: entry.updatedOn ?? ''});
    }
  }
  return files;
}

// ── Incremental user file sync ──────────────────────────────────────────
// Compares remote file timestamps against a local cache and only
// downloads files that are new or have changed. Removes local files that
// no longer exist on the remote side.
async function syncUserFilesIncremental(
  userDir: string, connId: string, connectorPath: string,
): Promise<string[]> {
  const userId = path.basename(userDir);
  if (!userFileCache.has(userId))
    userFileCache.set(userId, new Map());
  const cached = userFileCache.get(userId)!;

  // Fetch the full remote file list with timestamps
  const remoteFiles = await collectFiles(connId, 'agents');
  console.log(`user-files: remote has ${remoteFiles.length} file(s) in agents/`);

  const remoteSet = new Set<string>();
  const downloaded: string[] = [];

  // Download new or updated files
  for (const remote of remoteFiles) {
    remoteSet.add(remote.path);
    const cachedTs = cached.get(remote.path);
    console.log(`user-files: checking ${remote.path} — remote ts=${remote.updatedOn || '<empty>'}, cached ts=${cachedTs ?? 'none'}, match=${cachedTs === remote.updatedOn}`);
    if (cachedTs === remote.updatedOn)
      continue;

    const label = cachedTs ? 'updated' : 'new';
    console.log(`user-files: downloading ${label} file ${remote.path}`);
    try {
      const content = await downloadFile(connectorPath, remote.path) as string;
      const dest = safePath(userDir, remote.path);
      await fs.mkdir(path.dirname(dest), {recursive: true});
      await fs.writeFile(dest, content);
      cached.set(remote.path, remote.updatedOn);
      downloaded.push(remote.path.replace(/^agents\//, ''));
    }
    catch (e: any) {
      console.warn(`user-files: failed to download ${remote.path}:`, e.message);
    }
  }

  // Remove local files that were deleted on the remote side
  const agentsDir = path.join(userDir, 'agents');
  try {
    const localEntries = await fs.readdir(agentsDir, {recursive: true});
    for (const entry of localEntries) {
      const entryStr = String(entry);
      const fullPath = path.join(agentsDir, entryStr);
      const stat = await fs.stat(fullPath);
      if (!stat.isFile())
        continue;
      // Local files from user sync have paths like "agents/filename"
      const remotePath = `agents/${entryStr}`;
      // Only remove files that were synced from My files (not package files).
      // Package files have the pattern "{PackageName}-{filename}".
      if (!remoteSet.has(remotePath) && cached.has(remotePath)) {
        console.log(`user-files: removing deleted file ${remotePath}`);
        await fs.rm(fullPath, {force: true});
        cached.delete(remotePath);
      }
    }
  }
  catch { /* agents dir may not exist yet */ }

  return downloaded;
}

const AI_SKILLS_TAG = 'ai-skills';

// ── Shared skills sync ─────────────────────────────────────────────────
// Discovers connections tagged with "ai-skills" (server-filtered) and
// syncs their agents/ folders. Skips the user's own Home connection
// (already synced separately). Files are stored as
// agents/{connName}-{filename} to avoid collisions between sources.
function connectorPath(conn: ConnectionInfo): string {
  const ns = (conn.namespace ?? 'System:').replace(/:$/, '');
  return `${ns}:${conn.name}`;
}

async function syncSharedSkills(userDir: string, homeConnId: string | null): Promise<void> {
  console.log('shared-skills: discovering connections with ai-skills tag...');
  const connections = await request<ConnectionInfo[]>(
    'GET', `/public/v1/connections?tags=${AI_SKILLS_TAG}`,
  );
  if (!Array.isArray(connections) || !connections.length) {
    console.log('shared-skills: no tagged connections found');
    return;
  }
  console.log(`shared-skills: found ${connections.length} tagged connection(s)`);

  const userId = path.basename(userDir);
  if (!sharedSkillsCache.has(userId))
    sharedSkillsCache.set(userId, new Map());
  const cached = sharedSkillsCache.get(userId)!;

  const agentsDir = path.join(userDir, 'agents');
  const activeKeys = new Set<string>();

  for (const conn of connections) {
    if (conn.id === homeConnId) {
      console.log(`shared-skills: skipping own Home connection ${conn.id}`);
      continue;
    }

    const prefix = `shared-${conn.name}-`;
    console.log(`shared-skills: syncing connection "${conn.name}" (${conn.id})`);

    // List files at root — the shared connection itself is the skills folder,
    // unlike Home where skills live under agents/.
    const remoteFiles = await collectFiles(conn.id, '');
    if (!remoteFiles.length) {
      console.log(`shared-skills: no files in "${conn.name}"`);
      continue;
    }
    console.log(`shared-skills: found ${remoteFiles.length} file(s) in "${conn.name}"`);

    const cPath = connectorPath(conn);
    for (const remote of remoteFiles) {
      const cacheKey = `${conn.name}/${remote.path}`;
      activeKeys.add(cacheKey);
      if (cached.get(cacheKey) === remote.updatedOn)
        continue;

      const label = cached.has(cacheKey) ? 'updated' : 'new';
      console.log(`shared-skills: downloading ${label} file ${cacheKey}`);
      try {
        const content = await downloadFile(cPath, remote.path) as string;
        const fileName = remote.path.replace(/\//g, '-');
        const dest = safePath(agentsDir, `${prefix}${fileName}`);
        await fs.writeFile(dest, content);
        cached.set(cacheKey, remote.updatedOn);
      }
      catch (e: any) {
        console.warn(`shared-skills: failed to download ${cacheKey}:`, e.message);
      }
    }
  }

  // Clean up files from connections that are no longer shared or tagged
  for (const [key] of cached) {
    if (!activeKeys.has(key)) {
      const [connName, ...rest] = key.split('/');
      const fileName = rest.join('/').replace(/^agents\//, '').replace(/\//g, '-');
      const localFile = path.join(agentsDir, `shared-${connName}-${fileName}`);
      try {
        await fs.rm(localFile, {force: true});
        console.log(`shared-skills: removed stale file ${localFile}`);
      }
      catch { /* ignore */ }
      cached.delete(key);
    }
  }
}

// ── Package agent file sync ─────────────────────────────────────────────
// Downloads published Datagrok packages as ZIPs, extracts files from their
// agents/ folders, and stores them as /users/{userId}/agents/{PkgName}-{file}.
// Only re-downloads packages whose updatedOn timestamp has changed.

interface PackageInfo {
  id: string;
  name: string;
  updatedOn?: string;
}

async function syncPackageAgentFiles(userDir: string): Promise<void> {
  console.log('package-agents: fetching published packages...');
  const packages = await request<PackageInfo[]>('GET', '/packages/published/current');
  if (!Array.isArray(packages) || !packages.length) {
    console.log('package-agents: no published packages found');
    return;
  }
  console.log(`package-agents: found ${packages.length} published package(s)`);

  const userId = path.basename(userDir);
  if (!packageTimestamps.has(userId))
    packageTimestamps.set(userId, new Map());
  const cached = packageTimestamps.get(userId)!;

  const agentsDir = path.join(userDir, 'agents');
  const currentPackages = new Set<string>();

  for (const pkg of packages) {
    currentPackages.add(pkg.name);
    const ts = pkg.updatedOn ?? '';
    if (cached.get(pkg.name) === ts) {
      console.log(`package-agents: ${pkg.name} unchanged (ts=${ts}), skipping`);
      continue;
    }

    console.log(`package-agents: syncing ${pkg.name} (remote ts=${ts}, cached ts=${cached.get(pkg.name) ?? 'none'})`);
    try {
      const buf = await requestBinary('GET', `/packages/published/${pkg.id}/zip`);
      if (!buf || buf.length === 0) {
        console.warn(`package-agents: empty ZIP for ${pkg.name}, skipping`);
        continue;
      }

      const zip = new AdmZip(buf);
      const entries = zip.getEntries();
      const agentEntries = entries.filter((e) =>
        !e.isDirectory && e.entryName.startsWith('agents/'));

      if (!agentEntries.length) {
        console.log(`package-agents: ${pkg.name} has no agents/ files`);
        cached.set(pkg.name, ts);
        continue;
      }

      // Clean old files for this package before extracting new ones
      const prefix = `${pkg.name}-`;
      try {
        const existing = await fs.readdir(agentsDir);
        const removed = existing.filter((f) => f.startsWith(prefix));
        for (const f of removed)
          await fs.rm(path.join(agentsDir, f), {force: true});
        if (removed.length)
          console.log(`package-agents: removed ${removed.length} old file(s) for ${pkg.name}`);
      }
      catch { /* dir may not exist */ }

      // Extract agent files, flattening subdirectory paths with "-"
      for (const entry of agentEntries) {
        const rel = entry.entryName.replace(/^agents\//, '').replace(/\//g, '-');
        const dest = safePath(agentsDir, `${pkg.name}-${rel}`);
        await fs.writeFile(dest, entry.getData());
      }
      cached.set(pkg.name, ts);
      console.log(`package-agents: extracted ${agentEntries.length} file(s) from ${pkg.name}`);
    }
    catch (e: any) {
      console.warn(`package-agents: failed to sync ${pkg.name}:`, e.message);
    }
  }

  // Clean up files from packages that are no longer published
  for (const [name] of cached) {
    if (!currentPackages.has(name)) {
      const prefix = `${name}-`;
      try {
        const existing = await fs.readdir(agentsDir);
        const removed = existing.filter((f) => f.startsWith(prefix));
        for (const f of removed)
          await fs.rm(path.join(agentsDir, f), {force: true});
        if (removed.length)
          console.log(`package-agents: cleaned up ${removed.length} file(s) for removed package ${name}`);
      }
      catch { /* ignore */ }
      cached.delete(name);
    }
  }
}

// ── List agent files on disk ────────────────────────────────────────────
async function listAgentFiles(userDir: string): Promise<string[]> {
  const agentsDir = path.join(userDir, 'agents');
  try {
    const entries = await fs.readdir(agentsDir, {recursive: true});
    const files: string[] = [];
    for (const entry of entries) {
      const full = path.join(agentsDir, String(entry));
      const stat = await fs.stat(full);
      if (stat.isFile())
        files.push(String(entry));
    }
    console.log(`user-files: listed ${files.length} agent file(s) on disk`);
    return files;
  }
  catch {
    return [];
  }
}

// Sync scopes — determines which category of files to re-sync.
// 'all' syncs everything (used for initial sync and user messages).
// Scoped syncs only run the relevant step, keeping other caches intact.
export type SyncScope = 'all' | 'user-files' | 'packages' | 'shared';

// ── Main sync entry point ───────────────────────────────────────────────
// Called on every user message from handleMessage() with scope='all'.
// Event-driven syncs from the browser specify a narrower scope so only
// the relevant files are re-checked against the server.
// All syncs are incremental — timestamps are always compared, never cleared.
export async function syncUserFiles(
  apiUrl: string, apiKey: string, scope: SyncScope = 'all',
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

  const promise = doSync(apiUrl, apiKey, scope);
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
  apiUrl: string, apiKey: string, scope: SyncScope,
): Promise<{dir: string; files: string[]}> {
  const dir = getUserDir(apiKey);
  const userId = path.basename(dir);
  console.log(`user-files: starting sync for user ${userId} (scope=${scope})`);

  await ensureUserDir(apiKey);

  // ── Sync user files from My files/agents/ ──
  let homeConnId: string | null = null;
  if (scope === 'all' || scope === 'user-files') {
    await runWithContext({apiKey, apiUrl}, async () => {
      const home = await resolveHomeConnection();
      if (!home) {
        console.log('user-files: skipping user file sync (no Home connection)');
        return;
      }
      homeConnId = home.connId;
      const downloaded = await syncUserFilesIncremental(dir, home.connId, home.connectorPath);
      if (downloaded.length)
        console.log(`user-files: synced ${downloaded.length} user file(s)`);
      else
        console.log('user-files: all user files up-to-date');
    });
  }

  // ── Sync agent files from published packages ──
  if (scope === 'all' || scope === 'packages') {
    await runWithContext({apiKey, apiUrl}, async () => {
      try {
        await syncPackageAgentFiles(dir);
      }
      catch (e: any) {
        console.warn('package-agents: sync failed:', e.message);
      }
    });
  }

  // ── Sync skills from shared connections ──
  if (scope === 'all' || scope === 'shared') {
    // Need homeConnId to skip own connection — resolve if not already done
    if (!homeConnId) {
      await runWithContext({apiKey, apiUrl}, async () => {
        const home = await resolveHomeConnection();
        if (home)
          homeConnId = home.connId;
      });
    }
    await runWithContext({apiKey, apiUrl}, async () => {
      try {
        await syncSharedSkills(dir, homeConnId);
      }
      catch (e: any) {
        console.warn('shared-skills: sync failed:', e.message);
      }
    });
  }

  initialSyncDone.add(apiKey);

  const files = await listAgentFiles(dir);
  console.log(`user-files: sync complete — ${files.length} total agent file(s)`);
  return {dir, files};
}
