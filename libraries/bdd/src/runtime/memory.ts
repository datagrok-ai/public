/* What a worker's browser holds when a feature ends, written as one JSON line per feature to the
   file BDD_MEMORY_LOG names: the renderer, GPU and browser processes by the operating system's count
   (working set and private bytes on Windows, the resident set elsewhere), and the page's JS heap and
   DOM counts by the DevTools protocol. BDD_MEMORY_GC=1 adds the JS heap and the renderer after a
   forced garbage collection, which tells what the page keeps from what it has not collected yet. */
import {execFile} from 'node:child_process';
import {appendFile} from 'node:fs/promises';
import {promisify} from 'node:util';
import type {Page} from '@playwright/test';

const run = promisify(execFile);
const RUN_OPTIONS = {timeout: 10000, windowsHide: true};
const MB = 1024 * 1024;

type Process = {type: string; pid: number};
type Resident = {mb: number; privateMb?: number};

/** Resident memory by process id: PowerShell on Windows, ps elsewhere. */
async function resident(pids: number[]): Promise<Map<number, Resident>> {
  const out = new Map<number, Resident>();
  if (pids.length === 0)
    return out;
  if (process.platform === 'win32') {
    const command = `Get-Process -Id ${pids.join(',')} -ErrorAction SilentlyContinue | ` +
      'ForEach-Object { "$($_.Id) $($_.WorkingSet64) $($_.PrivateMemorySize64)" }';
    const {stdout} = await run('powershell', ['-NoProfile', '-Command', command], RUN_OPTIONS);
    for (const line of stdout.trim().split(/\r?\n/).filter(Boolean)) {
      const [pid, ws, priv] = line.trim().split(/\s+/).map(Number);
      out.set(pid, {mb: Math.round(ws / MB), privateMb: Math.round(priv / MB)});
    }
  }
  else {
    const {stdout} = await run('ps', ['-o', 'pid=,rss=', '-p', pids.join(',')], RUN_OPTIONS);
    for (const line of stdout.trim().split(/\n/).filter(Boolean)) {
      const [pid, kb] = line.trim().split(/\s+/).map(Number);
      out.set(pid, {mb: Math.round(kb / 1024)});
    }
  }
  return out;
}

async function browserProcesses(page: Page): Promise<Process[]> {
  const cdp = await page.context().browser()!.newBrowserCDPSession();
  try {
    const {processInfo} = await cdp.send('SystemInfo.getProcessInfo') as {processInfo: {type: string; id: number}[]};
    return processInfo.map((p) => ({type: p.type, pid: p.id}));
  }
  finally {
    await cdp.detach().catch(() => undefined);
  }
}

type ByType = Record<string, {pids: number[]; mb: number; privateMb?: number}>;

/** Memory by process type: the sum over the processes of a type, with their ids. */
async function processMemory(page: Page): Promise<ByType> {
  const processes = await browserProcesses(page);
  const memory = await resident(processes.map((p) => p.pid));
  const byType: ByType = {};
  for (const p of processes) {
    const m = memory.get(p.pid);
    const entry = byType[p.type] ??= {pids: [], mb: 0};
    entry.pids.push(p.pid);
    entry.mb += m?.mb ?? 0;
    if (m?.privateMb !== undefined)
      entry.privateMb = (entry.privateMb ?? 0) + m.privateMb;
  }
  return byType;
}

/** The browser's renderer processes in MB (private bytes on Windows, the resident set elsewhere): the
 * page's heap and DOM, the dedicated workers it started, which run inside it, and the few renderers
 * of the worker's one browser besides (a frame of another site, the spare Chromium keeps ready). */
const renderer = (byType: ByType): number => byType.renderer?.privateMb ?? byType.renderer?.mb ?? 0;

export async function rendererMb(page: Page): Promise<number> {
  return renderer(await processMemory(page));
}

async function pageMetrics(page: Page, gc: boolean): Promise<Record<string, number>> {
  const cdp = await page.context().newCDPSession(page);
  try {
    if (gc)
      await cdp.send('HeapProfiler.collectGarbage');
    await cdp.send('Performance.enable');
    const {metrics} = await cdp.send('Performance.getMetrics') as {metrics: {name: string; value: number}[]};
    await cdp.send('Performance.disable');
    const of = (name: string) => metrics.find((m) => m.name === name)?.value ?? 0;
    return {heapUsedMb: Math.round(of('JSHeapUsedSize') / MB), heapTotalMb: Math.round(of('JSHeapTotalSize') / MB),
      nodes: of('Nodes'), documents: of('Documents'), listeners: of('JSEventListeners'), frames: of('Frames')};
  }
  finally {
    await cdp.detach().catch(() => undefined);
  }
}

/** One line per feature end, and the renderer's MB it read; a reading that fails is written as its
 * error, never fails the feature. */
export async function logMemory(page: Page, fields: Record<string, unknown>): Promise<number | undefined> {
  const file = process.env.BDD_MEMORY_LOG;
  if (!file || page.isClosed())
    return undefined;
  const start = Date.now();
  const line: Record<string, unknown> = {t: new Date().toISOString(), worker: process.env.TEST_WORKER_INDEX,
    project: process.env.BDD_PROJECT, ...fields};
  let mb: number | undefined;
  try {
    const processes = await processMemory(page);
    line.processes = processes;
    mb = renderer(processes);
    line.page = await pageMetrics(page, false);
    if (process.env.BDD_MEMORY_GC === '1') {
      line.afterGc = await pageMetrics(page, true);
      line.processesAfterGc = await processMemory(page);
    }
  }
  catch (error) {
    line.error = String(error);
  }
  line.measureMs = Date.now() - start;
  await appendFile(file, JSON.stringify(line) + '\n').catch(() => undefined);
  return mb;
}
