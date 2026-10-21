// Bounds how long a turn may produce nothing, and kills agent process trees. A wedged agent may
// never act on a polite cancel, so the bound ends in a real kill.

import * as fs from 'node:fs';

/** Returned by `raceIdle` when the wrapped promise produced nothing in time. */
export const IDLE = Symbol('idle');

/** Resolves with the promise's value, or IDLE if it does not settle within `ms`.
 * Always clears its timer, so a long turn does not accumulate one per streamed event. */
export function raceIdle<T>(p: Promise<T>, ms: number): Promise<T | typeof IDLE> {
  let timer: NodeJS.Timeout;
  const idle = new Promise<typeof IDLE>((resolve) => {
    timer = setTimeout(() => resolve(IDLE), ms);
  });
  return Promise.race([p, idle]).finally(() => clearTimeout(timer));
}

/** Every process with its parent and command name, from /proc (Linux); empty anywhere else. */
export function processTable(): {pid: number; ppid: number; command: string}[] {
  const out: {pid: number; ppid: number; command: string}[] = [];
  let entries: string[];
  try {
    entries = fs.readdirSync('/proc');
  } catch {
    return out;
  }
  for (const entry of entries) {
    const pid = Number(entry);
    if (!pid) continue;
    try {
      const stat = fs.readFileSync(`/proc/${pid}/stat`, 'utf8');
      // The parenthesised command may contain spaces and ')': after the last ')' come state, then ppid.
      out.push({pid, ppid: Number(stat.slice(stat.lastIndexOf(')') + 2).split(' ')[1]),
        command: stat.slice(stat.indexOf('(') + 1, stat.lastIndexOf(')'))});
    } catch { /* exited while we were reading it */ }
  }
  return out;
}

export function descendantPids(root: number, table: {pid: number; ppid: number}[]): number[] {
  const children = new Map<number, number[]>();
  for (const {pid, ppid} of table)
    children.set(ppid, [...(children.get(ppid) ?? []), pid]);
  const out = [...(children.get(root) ?? [])];
  for (let i = 0; i < out.length; i++)
    out.push(...(children.get(out[i]) ?? []));
  return out;
}

/** SIGKILLs a process and its whole descendant tree. Hermes starts tool commands in sessions of
 * their own, so a process-group kill misses them; the tree is read before the first signal, while
 * the children still point at their parents. Returns the pids that could not be signalled. */
export function killTree(root: number): number[] {
  const pids = [root, ...descendantPids(root, processTable())];
  const failed: number[] = [];
  for (const pid of pids) {
    try {
      process.kill(pid, 'SIGKILL');
    } catch (e: any) {
      if (e?.code !== 'ESRCH')
        failed.push(pid);
    }
  }
  try {
    process.kill(-root, 'SIGKILL');
  } catch { /* not a group leader, or the group is gone */ }
  return failed;
}

/** SIGKILLs everything under a process (not the process itself). */
export function killDescendants(root: number): void {
  for (const pid of descendantPids(root, processTable())) {
    try {
      process.kill(pid, 'SIGKILL');
    } catch { /* already gone */ }
  }
}

function cmdline(pid: number): string {
  try {
    return fs.readFileSync(`/proc/${pid}/cmdline`, 'utf8').replace(/\0/g, ' ');
  } catch {
    return '';
  }
}

/** True for an agent process. Keyed on the module in the command line, not on the interpreter, so a
 * `python3` script a terminal tool call runs is never mistaken for an agent. */
export function isAgentProcess(commandLine: string): boolean {
  return /\s-m grokky_agent(\s|$)/.test(commandLine);
}

/** Reaps agent processes the session pool no longer tracks (a dispose that raced a spawn, a
 * kill the process survived): an untracked agent would hold memory with no turn left to time it out. */
export function startOrphanReaper(ownedPids: () => ReadonlySet<number>, intervalMs = 30000): NodeJS.Timeout {
  const timer = setInterval(() => {
    try {
      const owned = ownedPids();
      // Every process, not only this server's children: with tini as PID 1 an agent whose parent
      // link was lost is reparented to tini.
      for (const c of processTable()) {
        if (owned.has(c.pid) || !isAgentProcess(cmdline(c.pid)))
          continue;
        const failed = killTree(c.pid);
        console.warn(`reaper: SIGKILL orphaned ${c.command} (pid ${c.pid}) and its descendants — not owned by ` +
          `the session pool${failed.length ? `; could not signal ${failed.join(', ')}` : ''}`);
      }
    } catch (e: any) {
      console.warn('reaper: sweep failed:', e?.message ?? e);
    }
  }, intervalMs);
  timer.unref?.();
  return timer;
}
