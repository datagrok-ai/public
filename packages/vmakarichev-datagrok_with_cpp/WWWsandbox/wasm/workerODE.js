import {Module} from '../wasm/testODE';

onmessage = async function (evt) {

    console.log( await Module());

    postMessage('hi');
}