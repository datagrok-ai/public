import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {WORKSPACE} from './sync-utils';

const PACKAGES_DIR = 'packages';
const lastBuilt = new Map<string, string>();

export async function buildStagedWorkspace(
  userDir: string, installed: Iterable<string>,
): Promise<void> {
  const userId = path.basename(userDir);
  const sorted = [...installed].sort();
  const fingerprint = sorted.join(',');
  if (lastBuilt.get(userId) === fingerprint) {
    console.log(`staged-workspace: unchanged for ${userId}, skip`);
    return;
  }

  const stage = path.join(userDir, 'workspace');
  await fs.mkdir(stage, {recursive: true});

  const topLevel = (await fs.readdir(WORKSPACE)).filter((e) => e !== '.git');
  for (const entry of topLevel) {
    if (entry === PACKAGES_DIR)
      continue;
    const link = path.join(stage, entry);
    try {
      await fs.symlink(path.join(WORKSPACE, entry), link);
    } catch (e: any) {
      if (e.code !== 'EEXIST')
        throw e;
    }
  }

  const pkgsStage = path.join(stage, PACKAGES_DIR);
  await fs.mkdir(pkgsStage, {recursive: true});

  const target = new Set(sorted);
  let current: Set<string>;
  try {
    current = new Set(await fs.readdir(pkgsStage));
  } catch {
    current = new Set();
  }

  const toAdd = sorted.filter((n) => !current.has(n));
  const toRemove = [...current].filter((n) => !target.has(n));

  for (const name of toRemove) {
    try {
      await fs.unlink(path.join(pkgsStage, name));
    } catch (e: any) {
      console.warn(`staged-workspace: unlink ${name} failed:`, e.message);
    }
  }
  for (const name of toAdd) {
    const src = path.join(WORKSPACE, PACKAGES_DIR, name);
    try {
      await fs.access(src);
      await fs.symlink(src, path.join(pkgsStage, name));
    } catch {
      // package not in canonical clone (private/custom) — skip
    }
  }

  lastBuilt.set(userId, fingerprint);
  console.log(
    `staged-workspace: ${userId} +${toAdd.length} -${toRemove.length} ` +
    `(total ${sorted.length})`,
  );
}
