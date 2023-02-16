import {Module} from '../fibs/fib';

onmessage = function(e) {

    (async () => {
        console.log('WORKER!');

        Module().then(inst => inst._fib(10));        

        postMessage('Hi');
    })();         
}