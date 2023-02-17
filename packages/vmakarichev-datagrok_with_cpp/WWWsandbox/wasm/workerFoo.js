import {Module} from '../wasm/foo';

onmessage = async function (evt) {

    Module().then(foo => 
        postMessage(foo._foo(evt.data)));
}