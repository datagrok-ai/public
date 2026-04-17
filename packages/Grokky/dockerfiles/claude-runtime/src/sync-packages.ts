import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import AdmZip from 'adm-zip';
import {request, requestBinary} from './shared-api-client';
import {WORKSPACE} from './sync-utils';

interface PackageInfo {
  id: string;
  name: string;
  updatedOn?: string;
}

// Per-user package cache: packageName → updatedOn timestamp.
const packageCache = new Map<string, Map<string, string>>();

let workspacePackages: Set<string> | null = null;

async function getWorkspacePackages(): Promise<Set<string>> {
  if (workspacePackages)
    return workspacePackages;
  try {
    const entries = await fs.readdir(path.join(WORKSPACE, 'packages'), {withFileTypes: true});
    workspacePackages = new Set(entries.filter((e) => e.isDirectory()).map((e) => e.name));
    console.log(`package-agents: found ${workspacePackages.size} workspace package(s)`);
  } catch {
    workspacePackages = new Set();
  }
  return workspacePackages;
}

async function syncSinglePackage(userDir: string, pkg: PackageInfo, cached: Map<string, string>): Promise<void> {
  if (pkg.updatedOn && cached.get(pkg.name) === pkg.updatedOn) {
    console.log(`package-agents: ${pkg.name} up-to-date, skipping`);
    return;
  }
  console.log(`package-agents: syncing ${pkg.name}`);
  const buf = await requestBinary('GET', `/packages/published/${pkg.id}/zip`);
  if (!buf || buf.length === 0) {
    console.warn(`package-agents: empty ZIP for ${pkg.name}, skipping`);
    return;
  }

  const zip = new AdmZip(buf);
  const entries = zip.getEntries();
  const agentEntries = entries.filter((e) => !e.isDirectory && e.entryName.startsWith('agents/'));
  if (!agentEntries.length) {
    console.log(`package-agents: ${pkg.name} has no agents/ files`);
    return;
  }

  const agentsDir = path.join(userDir, 'agents');

  // Clean old files for this package before extracting new ones
  const prefix = `${pkg.name}-`;
  try {
    const existing = await fs.readdir(agentsDir);
    const removed = existing.filter((f) => f.startsWith(prefix));
    for (const f of removed)
      await fs.rm(path.join(agentsDir, f), {force: true});
    if (removed.length)
      console.log(`package-agents: removed ${removed.length} old file(s) for ${pkg.name}`);
  } catch { /* dir may not exist */ }

  // Extract agent files, flattening subdirectory paths with "-"
  for (const entry of agentEntries) {
    const rel = entry.entryName.replace(/^agents\//, '').replace(/\//g, '-');
    const dest = path.join(agentsDir, `${pkg.name}-${rel}`);
    await fs.writeFile(dest, entry.getData());
  }
  console.log(`package-agents: extracted ${agentEntries.length} file(s) from ${pkg.name}`);
  if (pkg.updatedOn)
    cached.set(pkg.name, pkg.updatedOn);
}

export async function syncPackages(userDir: string, packageName?: string): Promise<void> {
  const packages = await request<PackageInfo[]>('GET', '/packages/published/current');
  if (!Array.isArray(packages) || !packages.length) {
    console.log('package-agents: no published packages found');
    return;
  }

  const userId = path.basename(userDir);
  if (!packageCache.has(userId))
    packageCache.set(userId, new Map());
  const cached = packageCache.get(userId)!;

  const wsPackages = await getWorkspacePackages();

  if (packageName) {
    if (wsPackages.has(packageName)) {
      console.log(`package-agents: skipping workspace package ${packageName}`);
      return;
    }
    const pkg = packages.find((p) => p.name === packageName);
    if (!pkg) {
      console.log(`package-agents: package "${packageName}" not found among published`);
      return;
    }
    await syncSinglePackage(userDir, pkg, cached);
    return;
  }

  // Initial sync: fetch all, skip workspace packages
  console.log(`package-agents: syncing ${packages.length} published package(s), skipping ${wsPackages.size} workspace package(s)`);
  for (const pkg of packages) {
    if (wsPackages.has(pkg.name)) {
      console.log(`package-agents: skipping workspace package ${pkg.name}`);
      continue;
    }
    try {
      await syncSinglePackage(userDir, pkg, cached);
    } catch (e: any) {
      console.warn(`package-agents: failed to sync ${pkg.name}:`, e.message);
    }
  }
}
