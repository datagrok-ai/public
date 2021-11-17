import Worker from "./rdkit.worker.ts"; // .ts!

export class RdKitServiceWorkerHandler {

  worker: Worker;
  constructor() {
    this.worker = new Worker();
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
      this.worker.postMessage({op: op, args: args}, [channel.port2]);

    });
  }

  moduleInit = async (pathToRdkit: string) =>
    this.call('module::init', [pathToRdkit]);
  substructInit = async (dict: any) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query: string, querySmarts: string) =>
    this.call('substructLibrary::search', [query, querySmarts]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');

}