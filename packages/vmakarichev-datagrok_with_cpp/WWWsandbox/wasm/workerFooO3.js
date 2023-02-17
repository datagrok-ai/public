import {Module} from '../wasm/fooO3';

onmessage = async function (evt) {

    Module().then(foo => 
        postMessage(foo.ccall('foo', 'number', ['number'], [evt.data])));
}