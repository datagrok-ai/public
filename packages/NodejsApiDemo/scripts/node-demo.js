//name: NodeDemo
//description: nodejs script using the js-api (grok/DG) + a package python script + package files
//language: nodejs
//output: string user
//output: int cells
//output: double mean
//output: string fileText

// js-api bootstrap: resolves globally on images that ship datagrok-api;
// falls back to the staged copy on this stand
let dgapi;
try { dgapi = require('datagrok-api/datagrok'); }
catch (_) { dgapi = require('/tmp/datagrok-api/datagrok.js'); }
if (!globalThis.__dgReady)
  globalThis.__dgReady = dgapi.startDatagrok({apiUrl: @DATAGROK_API_URL, apiToken: @USER_API_KEY, detached: true});
await globalThis.__dgReady;
grok.dapi.token = @USER_API_KEY;
grok.dapi.root = @DATAGROK_API_URL;

user = (await grok.dapi.users.current()).friendlyName;
const t = DG.DataFrame.fromCsv('x\n1\n2\n3\n4');
cells = t.rowCount * t.columns.length;
mean = await grok.functions.call('NodejsApiDemo:PyMean', {df: t});
fileText = (await grok.dapi.files.readAsText('System:AppData/NodejsApiDemo/hello.txt')).trim();
