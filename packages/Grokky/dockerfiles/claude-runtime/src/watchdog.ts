// Turn watchdog — bounds how long a single turn may produce nothing.
//
// Why this exists: the agent CLI can wedge. Observed repeatedly on help/grounding turns — its
// `rg --files` child finishes writing more than a pipe buffer of paths and blocks, and the CLI
// busy-spins instead of draining it. The process then sits at 100% CPU forever. `container.json`
// gives this container **1 CPU**, so one wedged turn starves the Node server that would otherwise
// process the browser's abort: the session hangs, and every other session on the container hangs
// with it. A client-side timeout cannot rescue that, because the thing that must act is starved.
//
// So the bound lives here, next to the SDK loop, and ends with a real kill rather than a polite
// abort — an abort the wedged process never reads is not a bound.

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

/** Direct children of this process. `startTime` is the kernel's starttime field, which is
 * monotonic per boot and so orders children by age without needing a wall clock.
 * Linux-only (this is a container); returns nothing anywhere else rather than throwing. */
export function childProcesses(parentPid = process.pid): {pid: number; command: string; startTime: number}[] {
  const out: {pid: number; command: string; startTime: number}[] = [];
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
      // stat field 4 is ppid; the comm field is parenthesised and may itself contain spaces,
      // so split after the last ')' rather than on whitespace.
      const stat = fs.readFileSync(`/proc/${pid}/stat`, 'utf8');
      const tail = stat.slice(stat.lastIndexOf(')') + 2).split(' ');
      if (Number(tail[1]) !== parentPid) continue;
      const comm = stat.slice(stat.indexOf('(') + 1, stat.lastIndexOf(')'));
      // tail[0] is state, so the stat fields shift by 2 from the documented numbering;
      // starttime is field 22, i.e. tail[19].
      out.push({pid, command: comm, startTime: Number(tail[19]) || 0});
    } catch { /* process exited while we were reading it */ }
  }
  return out;
}

/** Last resort after an abort the wedged process never acted on. Kills our direct children,
 * which is the agent CLI and anything it spawned. Returns how many were signalled. */
export function killStrayChildren(reason: string): number {
  const children = childProcesses();
  for (const c of children) {
    try {
      process.kill(c.pid, 'SIGKILL');
      console.warn(`watchdog: SIGKILL ${c.command} (pid ${c.pid}) — ${reason}`);
    } catch { /* already gone */ }
  }
  return children.length;
}

/** Reaps agent processes that outlived the turn that spawned them.
 *
 * The turn watchdog only bounds a turn that is *in flight*. A wedged CLI can also be orphaned:
 * the SDK query ends (or errors) while the process it spawned keeps spinning, so there is no turn
 * left to time out and nothing ever reclaims it. Observed holding 100% of this container's single
 * CPU for eight minutes with no active session, starving every later turn.
 *
 * The interlock is a **count**, not an idle check. The SDK runs one agent process per in-flight
 * query, so any surplus beyond `activeCount()` is by definition orphaned — and unlike "sweep only
 * when idle", that stays true under continuous load, which is exactly when a benchmark or a busy
 * user never leaves an idle window for an idle-gated reaper to act in. The surplus is killed
 * oldest-first, since a live query's process is always the most recently started. */
export function startOrphanReaper(activeCount: () => number, intervalMs = 30000): NodeJS.Timeout {
  const timer = setInterval(() => {
    // Everything here runs on a bare timer callback, where an escaping throw is an uncaught
    // exception that takes the runtime down. A janitor must never be able to kill the process
    // it is cleaning up after.
    try {
      const agents = childProcesses()
        .filter((c) => /claude/.test(c.command))
        .sort((a, b) => a.startTime - b.startTime);
      const surplus = agents.length - Math.max(0, activeCount());
      if (surplus <= 0) return;
      for (const c of agents.slice(0, surplus)) {
        try {
          process.kill(c.pid, 'SIGKILL');
          console.warn(`reaper: SIGKILL orphaned ${c.command} (pid ${c.pid}) — ` +
            `${agents.length} agent process(es), ${activeCount()} query/queries in flight`);
        } catch { /* already gone */ }
      }
    } catch (e: any) {
      console.warn('reaper: sweep failed:', e?.message ?? e);
    }
  }, intervalMs);
  // Must not keep the process alive on its own.
  timer.unref?.();
  return timer;
}
