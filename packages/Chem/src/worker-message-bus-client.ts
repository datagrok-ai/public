/** Thrown into in-flight {@link WorkerMessageBusClient.call} promises when their worker is killed,
 * so awaiting callers can recognize cancellation and unwind quietly instead of erroring. */
export class WorkerCancelledError extends Error {
  constructor(reason = 'Worker terminated') {
    super(reason);
    this.name = 'WorkerCancelledError';
  }
}

export class WorkerMessageBusClient {
  _worker: Worker;
  /** In-flight calls, so {@link terminate} can reject them instead of leaving callers hung forever. */
  private _pending = new Set<{reject: (reason?: any) => void, port: MessagePort}>();
  constructor(worker: Worker) {
    this._worker = worker;
  }
  async call(op: string, args: any[] = []): Promise<any> {
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
      pending.reject(new WorkerCancelledError());
    }
    this._pending.clear();
    this._worker.terminate();
  }
}
