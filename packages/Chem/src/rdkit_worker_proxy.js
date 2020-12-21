class RdKitWorkerProxy {

  constructor(basePath, path = 'rdkit_worker.js') {
    this.worker = new Worker(basePath + 'src/' + path);
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

  moduleInit = async () =>
    this.call('module::init');
  substructInit = async (dict) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query) =>
    this.call('substructLibrary::search', [query]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');

}