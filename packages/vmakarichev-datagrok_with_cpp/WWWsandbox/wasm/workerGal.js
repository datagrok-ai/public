// Here, we use an approach from: https://dzone.com/articles/webassembly-web-workers 

onmessage = function (evt) {
  fetch(new URL('../wasm/fib.wasm', import.meta.url)).then(response => 
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
          //console.log(objInstance);
          postMessage(objInstance.exports.fib(evt.data));
        }
      );
    }
  );
}