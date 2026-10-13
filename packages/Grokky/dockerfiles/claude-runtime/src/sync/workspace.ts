import {execFile} from 'node:child_process';
import * as path from 'node:path';
import {promisify} from 'node:util';
import {WORKSPACE} from '../constants';
import {buildHelpIndex} from '../help-index';

const exec = promisify(execFile);
const SYNC_INTERVAL_MS = 30 * 60 * 1000;
const HELP_WORKSPACE = path.join(WORKSPACE, 'help');

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

async function syncRepo(dir: string): Promise<boolean> {
  const oldSha = (await exec('git', ['-C', dir, 'rev-parse', 'HEAD'])).stdout.trim();
  await exec('git', ['-C', dir, 'fetch', '--depth=1', 'origin', 'master']);
  await exec('git', ['-C', dir, 'reset', '--hard', 'FETCH_HEAD']);
  const newSha = (await exec('git', ['-C', dir, 'rev-parse', 'HEAD'])).stdout.trim();
  if (oldSha === newSha)
    return false;
  console.log(`workspace: synced ${path.basename(dir)} ${oldSha.slice(0, 7)} → ${newSha.slice(0, 7)}`);
  return true;
}

async function syncWorkspace(): Promise<void> {
  if (syncing)
    return;
  syncing = true;
  let resolve!: () => void;
  syncDone = new Promise((r) => (resolve = r));
  try {
    await allQueriesIdle;
    const publicChanged = await syncRepo(WORKSPACE);
    const helpChanged = await syncRepo(HELP_WORKSPACE);
    if (!publicChanged && !helpChanged) {
      console.log('workspace: already up to date');
      return;
    }

    try {
      buildHelpIndex(WORKSPACE);
    } catch (e: any) {
      console.warn('help-index: rebuild failed:', e.message);
    }
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
