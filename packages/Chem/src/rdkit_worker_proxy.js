import WorkerClass from './rdkit_worker.js';

export class RdKitWorkerProxy {

  constructor(basePath, path = 'rdkit_worker.js') {
    // this.worker = new Worker(basePath + 'src/' + path);
    // https://webpack.js.org/guides/web-workers/
    // this.worker = new Worker('./' + path, { type: 'module' });
    // https://dannadori.medium.com/how-to-bundle-webworker-in-npm-package-620dcec922e1
    this.worker = WorkerClass();
  }

  async call(op, args = []) {
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

  moduleInit = async (pathToRdkit) =>
    this.call('module::init', [pathToRdkit]);
  substructInit = async (dict) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query, querySmarts) =>
    this.call('substructLibrary::search', [query, querySmarts]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');

}