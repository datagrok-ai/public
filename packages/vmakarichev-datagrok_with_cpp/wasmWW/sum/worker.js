onmessage = function(e) {

    (async () => {
        const fetchPromise = fetch(new URL('../sum/sum.wasm', import.meta.url));
        const module = await WebAssembly.compileStreaming(fetchPromise);
        const instance = await WebAssembly.instantiate(module);
        
        const sum = instance.exports.sum;
        const memory = instance.exports.memory;

        let arr = new Int32Array(memory.buffer, 0, e.data.length);
        arr.set(e.data);
        
        postMessage(sum(arr.byteOffset, arr.length));
    })();         
}