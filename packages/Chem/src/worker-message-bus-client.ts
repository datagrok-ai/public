export class WorkerMessageBusClient {
  _worker: Worker;
  /** In-flight calls, so {@link terminate} can reject them instead of leaving callers hung forever. */
  private _pending = new Set<{reject: (reason?: any) => void, port: MessagePort}>();
  /** Resolves once this worker's module is initialized; gates calls so none reach a still-loading worker
   * (e.g. one that's reloading after a cancel-restart). Reset by each moduleInit. */
  protected _ready: Promise<unknown> = Promise.resolve();
  constructor(worker: Worker) {
    this._worker = worker;
  }
  async call(op: string, args: any[] = [], gated = true): Promise<any> {
    if (gated)
      await this._ready; // module init itself runs ungated, else it would wait on its own result
    return new Promise((res, rej) => {
      // https://advancedweb.hu/how-to-use-async-await-with-postmessage/
      // {op, args} -> {op, retval} | {error}
      const channel = new MessageChannel();
      const pending = {reject: rej, port: channel.port1};
      this._pending.add(pending);
      channel.port1.onmessage = ({data}) => {
        this._pending.delete(pending);
        channel.port1.close();
        if (data.error)
          rej(data.error);
        else
          res(data.retval);
      };
      this._worker.postMessage({op: op, args: args}, [channel.port2]);
    });
  }

  terminate() {
    // Reject in-flight calls so awaiting callers fail fast instead of hanging once the worker is gone.
    for (const pending of this._pending) {
      pending.port.close();
      pending.reject(new Error('Worker terminated'));
    }
    this._pending.clear();
    this._worker.terminate();
  }
}
