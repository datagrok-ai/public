/** Per-call task context: cooperative cancellation flag + progress reporting
 *  (mirrors DatagrokTask.update_state in datagrok_task.py). */
import {AsyncLocalStorage} from 'node:async_hooks';

import {FanoutPublisher} from './fanout-publisher';
import {Const, PipeClient} from './pipe-client';
import {logWarn} from './logger';

export const taskContextStorage = new AsyncLocalStorage<TaskContext>();

(globalThis as any).DG_TASK_PROGRESS = (percent: number | null, description: string): void =>
  taskContextStorage.getStore()?.progress(percent, description);

export class CancelledError extends Error {
  constructor(callId: string) {
    super(`Call ${callId} was cancelled`);
    this.name = 'CancelledError';
  }
}

export class TaskContext {
  readonly callId: string;
  cancelled = false;
  /** Set when the Canceled CALL was already published (by pidbox) so the runner
   *  suppresses the duplicate result. */
  cancelPublished = false;
  pipe: PipeClient | null = null;

  private readonly publisher: FanoutPublisher;

  constructor(callId: string, publisher: FanoutPublisher) {
    this.callId = callId;
    this.publisher = publisher;
  }

  /** Reports progress: 'PROGRESS {json}' over the pipe when one is open, otherwise a
   *  fanout 'progress' message (datlas parses {percent, description} either way).
   *  Throws CancelledError when the call was revoked; other errors are swallowed. */
  progress(percent: number | null, description: string): void {
    if (this.cancelled)
      throw new CancelledError(this.callId);
    try {
      const message = {'percent': percent, 'description': description};
      if (this.pipe != null && this.pipe.isConnected)
        this.pipe.sendText(`${Const.PROGRESS} ${JSON.stringify(message)}`);
      else
        this.publisher.publish(message, this.callId, 'progress').catch(() => {});
    }
    catch (e: any) {
      logWarn(`Progress update failed: ${e?.message ?? e}`, this.callId);
    }
  }
}
