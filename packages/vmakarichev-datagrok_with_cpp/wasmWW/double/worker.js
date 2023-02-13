onmessage = function(e) {

    (async () => {
        const fetchPromise = fetch(new URL('../double/double.wasm', import.meta.url));
        const module = await WebAssembly.compileStreaming(fetchPromise);
        const instance = await WebAssembly.instantiate(module);
        
        const doubleArray = instance.exports.doubleArray;
        const memory = instance.exports.memory;

        let arr = new Int32Array(memory.buffer, 0, e.data.length);
        arr.set(e.data);

        doubleArray(arr.byteOffset, arr.length)
        
        postMessage(arr);
    })();         
}