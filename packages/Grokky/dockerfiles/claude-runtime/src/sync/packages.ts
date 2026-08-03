import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {request, requestBinary, listPackageFolder} from '../shared-api-client';
import {WORKSPACE} from '../constants';
import {setInstalledPackages} from '../user/installed-packages';

export interface PackageInfo {
  id: string;
  name: string;
  version: string;
  buildHash: string;
  buildNumber: number;
  updatedOn?: string;
}

interface CacheEntry {
  buildHash: string;
  files: string[];
}

// Per-user, per-process cache: userId → (packageName → {buildHash, synced files}).
// In-memory only — /users is container-ephemeral, so persisting it would buy nothing.
const cache = new Map<string, Map<string, CacheEntry>>();

let workspacePackages: Set<string> | null = null;

async function getWorkspacePackages(): Promise<Set<string>> {
  if (workspacePackages)
    return workspacePackages;
  try {
    const entries = await fs.readdir(path.join(WORKSPACE, 'packages'), {withFileTypes: true});
    workspacePackages = new Set(entries.filter((e) => e.isDirectory()).map((e) => e.name.toLowerCase()));
    console.log(`package-agents: found ${workspacePackages.size} workspace package(s)`);
  } catch {
    workspacePackages = new Set();
  }
  return workspacePackages;
}

// Syncs one package's agents/ folder. Only fetches when the build changed.
async function syncSinglePackage(userDir: string, pkg: PackageInfo, userCache: Map<string, CacheEntry>): Promise<void> {
  const prev = userCache.get(pkg.name);
  if (prev?.buildHash === pkg.buildHash) {
    console.log(`package-agents: ${pkg.name} up-to-date (build=${pkg.buildHash}), skipping`);
    return;
  }

  console.log(`package-agents: syncing ${pkg.name} (build=${pkg.buildHash})`);
  const base = `/packages/published/files/${pkg.name}/${pkg.version}/${pkg.buildHash}/${pkg.buildNumber}`;
  const children = await listPackageFolder(`${base}/agents/`);
  const agentsDir = path.join(userDir, 'agents');

  // Drop the previous build's files before writing the current ones.
  if (prev?.files.length) {
    for (const f of prev.files)
      await fs.rm(path.join(agentsDir, f), {force: true});
    console.log(`package-agents: removed ${prev.files.length} old file(s) for ${pkg.name}`);
  }

  if (!children?.length) {
    userCache.set(pkg.name, {buildHash: pkg.buildHash, files: []});
    console.log(`package-agents: ${pkg.name} has no agents/ files`);
    return;
  }

  await fs.mkdir(agentsDir, {recursive: true});
  const files: string[] = [];
  for (const child of children) {
    if (child.endsWith('/')) {
      console.warn(`package-agents: ${pkg.name} skipping nested agents folder ${child}`);
      continue;
    }
    const buf = await requestBinary('GET', `${base}/agents/${child}`);
    const dest = `${pkg.name}-${child}`;
    await fs.writeFile(path.join(agentsDir, dest), buf);
    files.push(dest);
  }

  userCache.set(pkg.name, {buildHash: pkg.buildHash, files});
  console.log(`package-agents: extracted ${files.length} agent file(s) from ${pkg.name}`);
}

export async function syncPackages(userDir: string, packageName?: string): Promise<void> {
  const wsPackages = await getWorkspacePackages();
  if (packageName && wsPackages.has(packageName.toLowerCase())) {
    console.log(`package-agents: skipping workspace package ${packageName}`);
    return;
  }

  const packages = await request<PackageInfo[]>('GET', '/packages/published/current');
  if (!Array.isArray(packages) || !packages.length) {
    console.log('package-agents: no published packages found');
    return;
  }

  const userId = path.basename(userDir);
  setInstalledPackages(userId, packages);

  let userCache = cache.get(userId);
  if (!userCache) {
    userCache = new Map();
    cache.set(userId, userCache);
  }

  const agentsDir = path.join(userDir, 'agents');
  const installed = new Set(packages.map((p) => p.name));

  // Drop files + cache entries for packages no longer installed.
  let orphans = 0;
  for (const name of [...userCache.keys()]) {
    if (!installed.has(name)) {
      for (const f of userCache.get(name)!.files)
        await fs.rm(path.join(agentsDir, f), {force: true});
      userCache.delete(name);
      orphans++;
    }
  }
  if (orphans)
    console.log(`package-agents: removed ${orphans} orphan package(s)`);

  if (packageName) {
    const pkg = packages.find((p) => p.name === packageName);
    if (!pkg) {
      console.log(`package-agents: package "${packageName}" not found among published`);
      return;
    }
    try {
      await syncSinglePackage(userDir, pkg, userCache);
    } catch (e: any) {
      console.warn(`package-agents: failed to sync ${pkg.name}: ${e.message}`);
    }
    return;
  }

  console.log(`package-agents: syncing ${packages.length} published package(s), skipping ${wsPackages.size} workspace package(s)`);
  for (const pkg of packages) {
    if (wsPackages.has(pkg.name.toLowerCase())) {
      console.log(`package-agents: skipping workspace package ${pkg.name}`);
      continue;
    }
    try {
      await syncSinglePackage(userDir, pkg, userCache);
    } catch (e: any) {
      console.warn(`package-agents: failed to sync ${pkg.name}: ${e.message}`);
    }
  }
}
