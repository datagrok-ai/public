(async () => {
    const fetchPromise = fetch('fib.wasm');
    const module = await WebAssembly.compileStreaming(fetchPromise);
    const instance = await WebAssembly.instantiate(module);

    const func = instance.exports.fib; // this is the function required
 
    for(let i = 1; i <= 10; i++)
      console.log(i + " <-> " + fib(i));
   
})();