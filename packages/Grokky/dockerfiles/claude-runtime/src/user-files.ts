import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import AdmZip from 'adm-zip';
import {runWithContext, request, requestBinary, downloadFile} from './shared-api-client';

const USERS_DIR = '/users';
const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

const syncedUsers = new Set<string>();
const packageTimestamps = new Map<string, Map<string, string>>();

function userIdFromKey(apiKey: string): string {
  try {
    const token = apiKey.startsWith('Bearer ') ? apiKey.slice(7) : apiKey;
    const payload = JSON.parse(Buffer.from(token.split('.')[1], 'base64url').toString());
    return payload.sub ?? payload.usr?.id ?? 'default';
  }
  catch {
    return apiKey.length > 64 ? 'default' : apiKey;
  }
}

function getUserDir(apiKey: string): string {
  return path.join(USERS_DIR, userIdFromKey(apiKey));
}

async function ensureUserDir(apiKey: string): Promise<string> {
  const dir = getUserDir(apiKey);
  console.log(`user-files: ensuring user dir ${dir}`);
  await fs.mkdir(path.join(dir, 'agents'), {recursive: true});
  const link = path.join(dir, 'workspace');
  try {
    await fs.symlink(WORKSPACE, link);
  }
  catch (e: any) {
    if (e.code !== 'EEXIST')
      throw e;
  }
  return dir;
}

interface ConnectionInfo {
  id: string;
  name: string;
  namespace?: string;
}

async function resolveHomeConnection(): Promise<{connId: string; connectorPath: string} | null> {
  const conns = await request<ConnectionInfo[]>('GET', '/public/v1/connections?text=My%20files');
  const home = conns.find((c) => c.name === 'Home');
  if (!home)
    return null;
  const ns = (home.namespace ?? 'System:').replace(/:$/, '');
  return {connId: home.id, connectorPath: `${ns}:${home.name}`};
}

async function collectFiles(connId: string, dirPath: string): Promise<string[]> {
  let entries: any[];
  try {
    const filePath = dirPath ? `${dirPath.replace(/^\//, '')}/` : '';
    entries = await request<any[]>('GET', `/connectors/connections/${connId}/files/${filePath}`);
  }
  catch {
    return [];
  }
  if (!Array.isArray(entries))
    return [];

  const files: string[] = [];
  for (const entry of entries) {
    const name = entry.name ?? entry.fileName ?? '';
    if (!name)
      continue;
    const rel = dirPath ? `${dirPath}/${name}` : name;
    if (!entry.isFile)
      files.push(...await collectFiles(connId, rel));
    else
      files.push(rel);
  }
  return files;
}

interface PackageInfo {
  id: string;
  name: string;
  updatedOn?: string;
}

let standardPackageNames: Set<string> | null = null;

async function getStandardPackageNames(): Promise<Set<string>> {
  if (standardPackageNames)
    return standardPackageNames;
  try {
    const entries = await fs.readdir(path.join(WORKSPACE, 'packages'));
    standardPackageNames = new Set(entries);
  }
  catch {
    standardPackageNames = new Set();
  }
  return standardPackageNames;
}

async function syncPackageAgentFiles(userDir: string): Promise<void> {
  const packages = await request<PackageInfo[]>('GET', '/packages/published/current');
  if (!Array.isArray(packages) || !packages.length)
    return;

  const standard = await getStandardPackageNames();
  const custom = packages.filter((p) => !standard.has(p.name));
  if (!custom.length)
    return;

  const userId = path.basename(userDir);
  if (!packageTimestamps.has(userId))
    packageTimestamps.set(userId, new Map());
  const cached = packageTimestamps.get(userId)!;

  const agentsDir = path.join(userDir, 'agents');
  const currentPackages = new Set<string>();

  for (const pkg of custom) {
    currentPackages.add(pkg.name);
    const ts = pkg.updatedOn ?? '';
    if (cached.get(pkg.name) === ts)
      continue;

    console.log(`package-agents: syncing ${pkg.name} (updated: ${ts})`);
    try {
      const buf = await requestBinary('GET', `/packages/published/${pkg.id}/zip`);
      const zip = new AdmZip(buf);
      const entries = zip.getEntries();
      const agentEntries = entries.filter((e) =>
        !e.isDirectory && e.entryName.startsWith('agents/'));

      // Clean old files for this package
      const prefix = `${pkg.name}-`;
      try {
        const existing = await fs.readdir(agentsDir);
        for (const f of existing)
          if (f.startsWith(prefix))
            await fs.rm(path.join(agentsDir, f), {force: true});
      }
      catch { /* dir may not exist */ }

      for (const entry of agentEntries) {
        const rel = entry.entryName.replace(/^agents\//, '').replace(/\//g, '-');
        const dest = path.join(agentsDir, `${pkg.name}-${rel}`);
        await fs.writeFile(dest, entry.getData());
      }
      cached.set(pkg.name, ts);
      if (agentEntries.length)
        console.log(`package-agents: extracted ${agentEntries.length} file(s) from ${pkg.name}`);
    }
    catch (e: any) {
      console.warn(`package-agents: failed to sync ${pkg.name}:`, e.message);
    }
  }

  // Clean up files from removed packages
  for (const [name] of cached) {
    if (!currentPackages.has(name)) {
      const prefix = `${name}-`;
      try {
        const existing = await fs.readdir(agentsDir);
        for (const f of existing)
          if (f.startsWith(prefix))
            await fs.rm(path.join(agentsDir, f), {force: true});
      }
      catch { /* ignore */ }
      cached.delete(name);
      console.log(`package-agents: cleaned up files for removed package ${name}`);
    }
  }
}

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
    return files;
  }
  catch {
    return [];
  }
}

export async function syncUserFiles(
  apiUrl: string, apiKey: string, force?: boolean,
): Promise<{dir: string; files: string[]}> {
  const dir = getUserDir(apiKey);

  if (!syncedUsers.has(apiKey) || force) {
    console.log('user-files: starting sync...');
    await ensureUserDir(apiKey);

    const synced = await runWithContext({apiKey, apiUrl}, async () => {
      const home = await resolveHomeConnection();
      if (!home) {
        console.log('user-files: Home connection not found');
        return [];
      }
      console.log(`user-files: resolved connector path: ${home.connectorPath}`);

      const files = await collectFiles(home.connId, 'agents');
      console.log(`user-files: found ${files.length} file(s) in agents/`);
      if (!files.length)
        return [];

      const agentsDir = path.join(dir, 'agents');
      try {
        const entries = await fs.readdir(agentsDir);
        for (const entry of entries)
          await fs.rm(path.join(agentsDir, entry), {recursive: true, force: true});
      }
      catch { /* dir may be empty */ }

      const downloaded: string[] = [];
      for (const filePath of files) {
        try {
          const content = await downloadFile(home.connectorPath, filePath) as string;
          const dest = path.join(dir, filePath);
          await fs.mkdir(path.dirname(dest), {recursive: true});
          await fs.writeFile(dest, content);
          downloaded.push(filePath.replace(/^agents\//, ''));
        }
        catch (e: any) {
          console.warn(`Failed to download ${filePath}:`, e.message);
        }
      }
      return downloaded;
    });

    syncedUsers.add(apiKey);
    if (synced.length)
      console.log(`Synced ${synced.length} file(s) for user`);

    // Clear package timestamp cache on force sync so packages get re-downloaded
    if (force) {
      const userId = path.basename(dir);
      packageTimestamps.delete(userId);
    }
  }

  await runWithContext({apiKey, apiUrl}, async () => {
    try {
      await syncPackageAgentFiles(dir);
    }
    catch (e: any) {
      console.warn('package-agents: sync failed:', e.message);
    }
  });

  const files = await listAgentFiles(dir);
  return {dir, files};
}
