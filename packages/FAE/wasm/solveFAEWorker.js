// The following code is generated automatically.

import {exportFAE} from '../wasm/FAEForWebWorker';
import {cppWrapper} from '../wasm/callWasmForWebWorker.js';

onmessage = async function (evt) {
  exportFAE().then(module =>
    {
      let args = evt.data;
      let result = cppWrapper(module, args, 'solveFAE', 'number');
      postMessage({'callResult': result, 'args': args});
    } )
}
