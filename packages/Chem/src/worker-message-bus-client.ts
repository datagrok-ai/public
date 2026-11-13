export class WorkerMessageBusClient {
  _worker: Worker;
  /** In-flight calls, so {@link terminate} can reject them instead of leaving callers hung forever. */
  private _pending = new Set<{reject: (reason?: any) => void, port: MessagePort}>();
  /** Resolves once this worker's module is initialized; gates calls so none reach a still-loading worker
   * (e.g. one that's reloading after a cancel-restart). Reset by each moduleInit. */
  protected _ready: Promise<unknown> = Promise.resolve();
  /** Set by {@link terminate}: a call made after it rejects with this at once (the worker would never answer). */
  private _terminated: string | null = null;
  private _idleWaiters: (() => void)[] = [];
  constructor(worker: Worker) {
    this._worker = worker;
  }

  /** The calls sent and not answered yet (the Crux service gives a long job to the worker with the fewest). */
  get pendingCalls(): number {
    return this._pending.size;
  }

  /** Resolves once no call is in flight (at once when none is), or when the worker is terminated. */
  whenIdle(): Promise<void> {
    return this._pending.size === 0 ? Promise.resolve() : new Promise((resolve) => this._idleWaiters.push(resolve));
  }

  private notifyIdle(): void {
    const waiters = this._idleWaiters;
    this._idleWaiters = [];
    for (const resolve of waiters)
      resolve();
  }
  async call(op: string, args: any[] = [], gated = true): Promise<any> {
    if (gated)
      await this._ready; // module init itself runs ungated, else it would wait on its own result
    return new Promise((res, rej) => {
      if (this._terminated !== null) { // a call that waited for the module, or came, after the worker was ended
        rej(new Error(this._terminated));
        return;
      }
      // https://advancedweb.hu/how-to-use-async-await-with-postmessage/
      // {op, args} -> {op, retval} | {error}
      const channel = new MessageChannel();
      const pending = {reject: rej, port: channel.port1};
      this._pending.add(pending);
      channel.port1.onmessage = ({data}) => {
        this._pending.delete(pending);
        channel.port1.close();
        if (this._pending.size === 0)
          this.notifyIdle();
        if (data.error)
          rej(data.error);
        else
          res(data.retval);
      };
      this._worker.postMessage({op: op, args: args}, [channel.port2]);
    });
  }

  /** Ends the worker; its in-flight calls, and any made after, reject with `reason` ('Worker terminated' by
   * default). */
  terminate(reason = 'Worker terminated') {
    this._terminated = reason;
    // Reject in-flight calls so awaiting callers fail fast instead of hanging once the worker is gone.
    for (const pending of this._pending) {
      pending.port.close();
      pending.reject(new Error(reason));
    }
    this._pending.clear();
    this.notifyIdle();
    this._worker.terminate();
  }
}
