export class WorkerMessageBusClient {
  _worker: Worker;
  constructor(worker: Worker) {
    this._worker = worker;
  }
  async call(op: string, args: any[] = []) {
    return new Promise((res, rej) => {
      // https://advancedweb.hu/how-to-use-async-await-with-postmessage/
      // {op, args} -> {op, retval} | {error}
      const channel = new MessageChannel();
      channel.port1.onmessage = ({data}) => {
        channel.port1.close();
        if (data.error) {
          rej(data.error);
        } else {
          res(data.retval);
        }
      };
      this._worker.postMessage({op: op, args: args}, [channel.port2]);
    });
  }
}
