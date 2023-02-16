var g_importObject = {
    'env': {
      'memoryBase': 0,       
      'tableBase': 0,       
      'memory': new WebAssembly.Memory({ initial: 256 }),       
      'table': new WebAssembly.Table({ initial: 0, element: 'anyfunc' })
    }
  };    
/*  
  // The WebAssembly module instance that we'll be working with   
  var g_objInstance = null;     
  

  // Listen for messages from the main thread. Because all messages
  // to this thread come through this method, we need a way to know
  // what is being asked of us which is why we included the
  // MessagePurpose property.   
  onmessage = function (evt) {

    // If we've been asked to call the module's Add method then...     
    var objData = evt.data;     
    var sMessagePurpose = objData.MessagePurpose;  

    if (sMessagePurpose === "Compute") {

      // Call the add method in the WebAssembly module and pass the       
      // result back to the main thread

      var iResult = g_objInstance.exports._fib(objData.Val); 

      console.log(`Fib is ${iResult.toString()}.`);

      postMessage(`Fib is ${iResult.toString()}.`);
    } // If we've been passed a compiled WebAssembly module then... 

    else if (sMessagePurpose === "CompiledModule") {

      // NOTE: Unlike when we pass in the bytes to instantiate, we
      // don't have a separate 'instance' and 'modules' object
      // returned in this case since we started out with the module
      // object. We're only passed back the instance in this case.       
      WebAssembly.instantiate(objData.WasmModule, g_importObject).then(instance =>
        // Hold onto the module's instance so that we can reuse it
        g_objInstance = instance
      );

    }
  }
  */

onmessage = function(e) {

    (async () => {
      
        
        const g_objInstance = await WebAssembly.instantiate(module);//, g_importObject);
        
        const func = g_objInstance.exports.fib; // this is the function required
        
        postMessage(func(e.data));
    })();        
}