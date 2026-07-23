import * as fs from 'node:fs/promises';
import * as path from 'node:path';
import {WORKSPACE} from '../constants';

const PACKAGES_DIR = 'packages';
// Repo dev-tooling that must not reach the model's cwd: the root CLAUDE.md instructs a
// `.kg/scripts/qq.py` knowledge-graph query first (self-installs a venv, ~30 s) — a measured
// help-turn timeout mechanism if the SDK auto-loads it. The model's instructions come from the
// system prompt, not the repo.
const EXCLUDED_ENTRIES = new Set(['CLAUDE.md', '.claude', '.kg']);
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

  const topLevel = (await fs.readdir(WORKSPACE)).filter((e) => e !== '.git' && !EXCLUDED_ENTRIES.has(e));
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
  // Drop previously-created links for now-excluded entries (existing user dirs).
  for (const entry of EXCLUDED_ENTRIES) {
    try {
      await fs.unlink(path.join(stage, entry));
    } catch { /* not present — fine */ }
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
