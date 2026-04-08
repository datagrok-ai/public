import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {runWithContext, request, downloadFile} from './shared-api-client';

const USERS_DIR = '/users';
const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

const syncedUsers = new Set<string>();

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
  }

  const files = await listAgentFiles(dir);
  return {dir, files};
}
