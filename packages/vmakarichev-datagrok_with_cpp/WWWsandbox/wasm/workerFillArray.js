// Here, we use an approach from: https://dzone.com/articles/webassembly-web-workers 

onmessage = function (evt) {
    fetch(new URL('../wasm/fillArray.wasm', import.meta.url)).then(response => 
      response.arrayBuffer()
    ).then(bytes =>
      WebAssembly.compile(bytes)
    ).then(WasmModule => {
        var importObject = {
          'env': {
            'memoryBase': 0,
            'tableBase': 0,
            'memory': new WebAssembly.Memory({ initial: 256 }),
            'table': new WebAssembly.Table({ initial: 0, element: 'anyfunc' })
          }
        };
  
        var objInstance = null;
  
        WebAssembly.instantiate(WasmModule, importObject
          ).then(instance => {
            objInstance = instance;
            console.log('Array in worker:');
            console.log(evt.data.array);
            console.log(objInstance);
            console.log(objInstance.exports.a);

            const fillArray = instance.exports.c;
            const memory = instance.exports.a;

            let arr = new Float32Array(memory.buffer, 0, evt.data.array.length);

            fillArray(arr.byteOffset, arr.length);
            
            console.log(arr);

            postMessage(/*objInstance.exports.fib(evt.data)*/'Hi');
          }
        );
      }
    );
  }