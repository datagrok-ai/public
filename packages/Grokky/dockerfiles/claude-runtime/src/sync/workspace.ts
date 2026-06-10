import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import {WORKSPACE} from '../constants';
import {invalidatePackageKnowledgeCache} from '../package-knowledge-tool';

const exec = promisify(execFile);
const SYNC_INTERVAL_MS = 30 * 60 * 1000;

let activeCount = 0;
let syncing = false;
let syncDone: Promise<void> = Promise.resolve();
let resolveIdle: (() => void) | null = null;
let allQueriesIdle: Promise<void> = Promise.resolve();

export function awaitWorkspaceSync(): Promise<void> {
  return syncDone;
}

export function markQueryStart(): void {
  if (activeCount === 0)
    allQueriesIdle = new Promise((r) => (resolveIdle = r));
  activeCount++;
}

export function markQueryEnd(): void {
  activeCount--;
  if (activeCount === 0 && resolveIdle) {
    resolveIdle();
    resolveIdle = null;
  }
}

async function syncWorkspace(): Promise<void> {
  if (syncing)
    return;
  syncing = true;
  let resolve!: () => void;
  syncDone = new Promise((r) => (resolve = r));
  try {
    await allQueriesIdle;
    const oldSha = (await exec('git', ['-C', WORKSPACE, 'rev-parse', 'HEAD'])).stdout.trim();
    // Sync the branch the workspace was cloned from (the Dockerfile's `git clone -b`),
    // not a hardcoded master. A pre-merge image cloned from a feature branch that
    // carries files absent on master (e.g. .kg) would otherwise be reset to master
    // and lose them. Falls back to master if HEAD is detached.
    const branchRef = (await exec('git', ['-C', WORKSPACE, 'rev-parse', '--abbrev-ref', 'HEAD'])).stdout.trim();
    const branch = branchRef && branchRef !== 'HEAD' ? branchRef : 'master';
    await exec('git', ['-C', WORKSPACE, 'fetch', '--depth=1', 'origin', branch]);
    await exec('git', ['-C', WORKSPACE, 'reset', '--hard', 'FETCH_HEAD']);
    const newSha = (await exec('git', ['-C', WORKSPACE, 'rev-parse', 'HEAD'])).stdout.trim();

    if (oldSha === newSha) {
      console.log('workspace: already up to date');
      return;
    }

    const {stdout} = await exec('git', ['-C', WORKSPACE, 'diff', '--name-only', oldSha, newSha]);
    const changed = stdout.split('\n');
    if (changed.some((f) => /^packages\/[^/]+\/agents\/package-knowledge\.yaml$/.test(f))) {
      invalidatePackageKnowledgeCache();
      console.log('workspace: knowledge cache invalidated');
    }
    console.log(`workspace: synced ${oldSha.slice(0, 7)} → ${newSha.slice(0, 7)}`);
  } catch (e: any) {
    console.warn('workspace: sync failed:', e.message);
  } finally {
    resolve();
    syncing = false;
  }
}

export function startWorkspaceSync(): void {
  setInterval(syncWorkspace, SYNC_INTERVAL_MS);
}
