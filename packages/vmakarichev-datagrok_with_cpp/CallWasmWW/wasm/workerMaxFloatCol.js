import {exportLib} from '../wasm/LibForWebWorker';

// Import utilities
import {cppWrapper} from './callWasmUtils';

onmessage = async function (evt) {
  exportLib().then(module => 
    {     
        let args = evt.data
        let result = cppWrapper(module, args, 'maxFloatCol', 'number');
        postMessage({'callResult': result, 'args': args});
    } )
}