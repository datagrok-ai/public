/// The change set `tests-for --changed` answers about (change-tests/plan.md § Ops): what git says differs from a
/// base in the monorepo and in the `public/` submodule, working tree and untracked files included, the public
/// paths under the prefix the graph gives them.
import * as path from 'path';
import {spawnSync} from 'child_process';

export interface ChangedFile {
  path: string;
  repo: 'core' | 'public';
}

export interface ChangeSet {
  files: ChangedFile[];
  /** The base each repository was diffed against. */
  base: Record<string, string>;
  notes: string[];
}

/** The merge base with master, or with origin/master, or the parent commit when neither branch exists. */
const BASES = ['master', 'origin/master'];

function git(cwd: string, ...args: string[]): string | undefined {
  const r = spawnSync('git', ['-C', cwd, '-c', 'core.quotepath=false', ...args], {encoding: 'utf8', maxBuffer: 32 * 1024 * 1024});
  return r.status === 0 ? r.stdout.trim() : undefined;
}

function mergeBase(cwd: string, name: string, notes: string[]): string {
  for (const branch of BASES) {
    const base = git(cwd, 'merge-base', branch, 'HEAD');
    if (base) return base;
  }
  notes.push(`${name}: neither master nor origin/master exists; the base is HEAD~1`);
  return 'HEAD~1';
}

/**
 * The files changed against [ref] (default: the merge base with master) in the monorepo at [repoRoot] and in its
 * `public/` submodule: `git diff --name-only <base>` (committed and working-tree changes alike) plus the untracked
 * files. A monorepo ref names nothing inside the submodule, so the public base is the gitlink that ref recorded
 * (`<ref>:public`) when it resolves, and the submodule's own merge base otherwise.
 */
export function changeSet(repoRoot: string, ref?: string): ChangeSet {
  const notes: string[] = [];
  const base: Record<string, string> = {};
  const files = new Map<string, ChangedFile>();
  const roots = new Set<string>();
  for (const repo of ['core', 'public'] as const) {
    const cwd = repo === 'core' ? repoRoot : path.join(repoRoot, 'public');
    const top = git(cwd, 'rev-parse', '--show-toplevel');
    if (top === undefined || roots.has(top)) continue;
    roots.add(top);
    const from = repo === 'core' ? ref ?? mergeBase(cwd, repo, notes)
      : (ref === undefined ? undefined : git(repoRoot, 'rev-parse', `${ref}:public`)) ?? mergeBase(cwd, repo, notes);
    base[repo] = from;
    const diff = git(cwd, 'diff', '--name-only', from);
    if (diff === undefined) {
      notes.push(`${repo}: git diff ${from} failed; its changes are not in the answer`);
      continue;
    }
    const prefix = repo === 'public' ? 'public/' : '';
    for (const line of [...diff.split('\n'), ...(git(cwd, 'ls-files', '--others', '--exclude-standard') ?? '').split('\n')]) {
      const file = line.trim();
      if (file && file !== 'public') files.set(prefix + file, {path: prefix + file, repo});
    }
  }
  return {files: [...files.values()].sort((a, b) => a.path < b.path ? -1 : 1), base, notes};
}
