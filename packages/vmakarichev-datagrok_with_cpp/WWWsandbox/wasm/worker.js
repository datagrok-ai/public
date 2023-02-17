onmessage = function(e) {

    (async () => {
        const fetchPromise = fetch(new URL('../wasm/fib.wasm', import.meta.url));
        const module = await WebAssembly.compileStreaming(fetchPromise);
        const instance = await WebAssembly.instantiate(module);
        
        const func = instance.exports.fib; // this is the function required
        
        postMessage(func(e.data));
    })();         
}
