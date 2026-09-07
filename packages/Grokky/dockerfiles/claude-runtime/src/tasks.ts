/** One record per queued Grokky:aiChatTurnTask call admitting one chat turn; the queue carries
 *  only this admission, the turn streams over the browser WebSocket. Fail-open: a turn on a stand
 *  with no live workers is delayed by CLAIM_WAIT_MS, never blocked — but while workers are seen
 *  polling, an unclaimed task means the queue is busy, not broken, and the turn keeps waiting. */

export interface TaskEnd {status: 'ended' | 'released'}

interface Task {
  sessionId: string;
  claimed: boolean;
  end?: TaskEnd;
  claimWaiters: ((claimed: boolean) => void)[];
  endWaiters: ((end: TaskEnd) => void)[];
  deleteTimer?: NodeJS.Timeout;
}

// Backstop for a task whose turn never arrives (worker claimed, browser send failed).
const ORPHAN_TTL_MS = 30 * 60_000;
// How long an ended task lingers so a late worker poll still sees the terminal status.
const END_LINGER_MS = 5 * 60_000;
// One long-poll leg; the worker re-polls on {status: 'running'}. Kept well under proxy timeouts.
const CLAIM_POLL_MS = 25_000;
const CLAIM_WAIT_MS = 60_000;
// Any /task/claim poll this recent proves the worker pool is alive (a slot-holding worker
// re-polls every CLAIM_POLL_MS, so a busy queue refreshes this continuously).
const WORKER_LIVE_MS = 90_000;
// Bounds a wait whose workers died mid-queue; matches the worker's own HOLD_CAP_MS.
const CLAIM_WAIT_CAP_MS = 30 * 60_000;

let lastWorkerPollAt = 0;

const tasks = new Map<string, Task>();

function armDelete(taskId: string, t: Task, ms: number): void {
  clearTimeout(t.deleteTimer);
  t.deleteTimer = setTimeout(() => tasks.delete(taskId), ms);
  t.deleteTimer.unref?.();
}

function getOrCreate(taskId: string, sessionId: string): Task {
  let t = tasks.get(taskId);
  if (!t) {
    t = {sessionId, claimed: false, claimWaiters: [], endWaiters: []};
    tasks.set(taskId, t);
    armDelete(taskId, t, ORPHAN_TTL_MS);
  }
  return t;
}

/** One worker long-poll leg: claim (idempotent), then wait up to CLAIM_POLL_MS for the turn to end;
 * {status: 'running'} means re-poll. */
export function claimTask(taskId: string, sessionId: string): Promise<TaskEnd | {status: 'running'}> {
  lastWorkerPollAt = Date.now();
  const t = getOrCreate(taskId, sessionId);
  if (!t.claimed) {
    t.claimed = true;
    t.claimWaiters.splice(0).forEach((w) => w(true));
  }
  if (t.end)
    return Promise.resolve(t.end);
  return new Promise((resolve) => {
    const onEnd = (end: TaskEnd) => {
      clearTimeout(timer);
      resolve(end);
    };
    const timer = setTimeout(() => {
      const i = t.endWaiters.indexOf(onEnd);
      if (i >= 0)
        t.endWaiters.splice(i, 1);
      resolve({status: 'running'});
    }, CLAIM_POLL_MS);
    t.endWaiters.push(onEnd);
  });
}

/** Turn side: hold until the worker claims the task; false on fail-open timeout or abort. */
export function waitForClaim(taskId: string, sessionId: string, signal: AbortSignal): Promise<boolean> {
  const t = getOrCreate(taskId, sessionId);
  // A turn may outlive the orphan backstop; endTask (runTurn's finally) is the terminal now.
  clearTimeout(t.deleteTimer);
  if (t.claimed || signal.aborted)
    return Promise.resolve(t.claimed);
  const startedAt = Date.now();
  return new Promise((resolve) => {
    let timer: NodeJS.Timeout;
    const done = (claimed: boolean) => {
      clearTimeout(timer);
      resolve(claimed);
    };
    const arm = () => {
      timer = setTimeout(() => {
        if (Date.now() - lastWorkerPollAt < WORKER_LIVE_MS && Date.now() - startedAt < CLAIM_WAIT_CAP_MS)
          arm();
        else
          done(false);
      }, CLAIM_WAIT_MS);
    };
    arm();
    t.claimWaiters.push(done);
    signal.addEventListener('abort', () => done(false), {once: true});
  });
}

/** Terminal for every turn path — called from runTurn's finally. Idempotent. */
export function endTask(taskId: string): void {
  const t = tasks.get(taskId);
  if (!t || t.end)
    return;
  t.end = {status: 'ended'};
  armDelete(taskId, t, END_LINGER_MS);
  t.endWaiters.splice(0).forEach((w) => w(t.end!));
}

/** Platform-side cancel (task revoked): returns the sessionId whose turn should be aborted, if still live. */
export function releaseTask(taskId: string): string | null {
  const t = tasks.get(taskId);
  if (!t || t.end)
    return null;
  t.end = {status: 'released'};
  armDelete(taskId, t, END_LINGER_MS);
  t.endWaiters.splice(0).forEach((w) => w(t.end!));
  t.claimWaiters.splice(0).forEach((w) => w(false));
  return t.sessionId;
}
