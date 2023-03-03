// The following code is generated automatically.

import {exportCppLib} from '../wasm/CppLibForWebWorker';
import {cppWrapper} from '../wasm/callWasmForWebWorker.js';

onmessage = async function (evt) {
  exportCppLib().then(module => 
    {
      let args = evt.data;
      let result = cppWrapper(module, args, 'addIntCols', 'number');
      postMessage({'callResult': result, 'args': args});
    } )
}