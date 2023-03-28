// The following code is generated automatically.

import {exportSVMlib} from '../wasm/SVMlibForWebWorker';
import {cppWrapper} from '../wasm/callWasmForWebWorker.js';

onmessage = async function (evt) {
  exportSVMlib().then(module => 
    {
      let args = evt.data;
      let result = cppWrapper(module, args, 'demoLinearKernel', 'number');
      postMessage({'callResult': result, 'args': args});
    } )
}