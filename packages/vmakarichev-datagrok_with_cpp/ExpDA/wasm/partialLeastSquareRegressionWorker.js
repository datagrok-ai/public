// The following code is generated automatically.

import {exportEDALib} from '../wasm/EDALibForWebWorker';
import {cppWrapper} from '../wasm/callWasmForWebWorker.js';

onmessage = async function (evt) {
  exportEDALib().then(module => 
    {
      let args = evt.data;
      let result = cppWrapper(module, args, 'partialLeastSquareRegression', 'number');
      postMessage({'callResult': result, 'args': args});
    } )
}