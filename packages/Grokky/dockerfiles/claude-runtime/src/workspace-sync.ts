import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import {WORKSPACE} from './sync-utils';
import {invalidatePackageKnowledgeCache} from './package-knowledge-tool';

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
    await exec('git', ['-C', WORKSPACE, 'fetch', '--depth=1', 'origin', 'master']);
    await exec('git', ['-C', WORKSPACE, 'reset', '--hard', 'FETCH_HEAD']);
    invalidatePackageKnowledgeCache();
    console.log('workspace: synced to latest origin/master');
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
