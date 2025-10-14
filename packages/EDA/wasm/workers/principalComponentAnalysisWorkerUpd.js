// The following code is generated automatically.

import {exportEDA} from '../EDAForWebWorker';
import {cppWrapper} from '../callWasmForWebWorker.js';

onmessage = async function(evt) {
  exportEDA().then((module) => {
    const args = evt.data;
    try {
      const result = cppWrapper(module, args, 'principalComponentAnalysis', 'number');
      postMessage({'callResult': result, 'args': args});
    } catch (e) {
      postMessage({'callResult': -1});
    }
  } );
};
