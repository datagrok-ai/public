class RdKitWorkerProxy {

  constructor(path = 'http://localhost:8080/api/packages/published/files/Chem/0.3.0-8b678c/src/rdkit_worker.js') {
    
    // https://advancedweb.hu/how-to-use-async-await-with-postmessage/
    //this.channel = new MessageChannel();
    this.worker = new Worker(path);
    
  }

  async call(op, args = []) { return new Promise((res, rej) => {
    
    // {op, args} -> {op, retval} | {error}
    const channel = new MessageChannel(); 
    channel.port1.onmessage = ({ data }) => {
      channel.port1.close();
      if (data.error) {
        rej(data.error);
      } else {
        res(data.retval);
      }
    };
    this.worker.postMessage({ op: op, args: args }, [channel.port2]);
    
  }); }
  
  moduleInit = async () =>
    this.call('module::init');
  substructInit = async (dict) =>
    this.call('substructLibrary::init', [dict]);
  substructSearch = async (query) =>
    this.call('substructLibrary::search', [query]);
  substructDeinit = async () =>
    this.call('substructLibrary::deinit');
  
}