
var exportEDA = (() => {
  var _scriptDir = typeof document !== 'undefined' && document.currentScript ? document.currentScript.src : undefined;
  
  return (
function(exportEDA) {
  exportEDA = exportEDA || {};

var Module=typeof exportEDA!="undefined"?exportEDA:{};var readyPromiseResolve,readyPromiseReject;Module["ready"]=new Promise(function(resolve,reject){readyPromiseResolve=resolve;readyPromiseReject=reject});var moduleOverrides=Object.assign({},Module);var arguments_=[];var thisProgram="./this.program";var quit_=(status,toThrow)=>{throw toThrow};var ENVIRONMENT_IS_WEB=typeof window=="object";var ENVIRONMENT_IS_WORKER=typeof importScripts=="function";var ENVIRONMENT_IS_NODE=typeof process=="object"&&typeof process.versions=="object"&&typeof process.versions.node=="string";var scriptDirectory="";function locateFile(path){if(Module["locateFile"]){return Module["locateFile"](path,scriptDirectory)}return scriptDirectory+path}var read_,readAsync,readBinary,setWindowTitle;if(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER){if(ENVIRONMENT_IS_WORKER){scriptDirectory=self.location.href}else if(typeof document!="undefined"&&document.currentScript){scriptDirectory=document.currentScript.src}if(_scriptDir){scriptDirectory=_scriptDir}if(scriptDirectory.indexOf("blob:")!==0){scriptDirectory=scriptDirectory.substr(0,scriptDirectory.replace(/[?#].*/,"").lastIndexOf("/")+1)}else{scriptDirectory=""}{read_=url=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.send(null);return xhr.responseText};if(ENVIRONMENT_IS_WORKER){readBinary=url=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.responseType="arraybuffer";xhr.send(null);return new Uint8Array(xhr.response)}}readAsync=(url,onload,onerror)=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,true);xhr.responseType="arraybuffer";xhr.onload=()=>{if(xhr.status==200||xhr.status==0&&xhr.response){onload(xhr.response);return}onerror()};xhr.onerror=onerror;xhr.send(null)}}setWindowTitle=title=>document.title=title}else{}var out=Module["print"]||console.log.bind(console);var err=Module["printErr"]||console.warn.bind(console);Object.assign(Module,moduleOverrides);moduleOverrides=null;if(Module["arguments"])arguments_=Module["arguments"];if(Module["thisProgram"])thisProgram=Module["thisProgram"];if(Module["quit"])quit_=Module["quit"];var wasmBinary;if(Module["wasmBinary"])wasmBinary=Module["wasmBinary"];var noExitRuntime=Module["noExitRuntime"]||true;if(typeof WebAssembly!="object"){abort("no native wasm support detected")}var wasmMemory;var ABORT=false;var EXITSTATUS;var UTF8Decoder=typeof TextDecoder!="undefined"?new TextDecoder("utf8"):undefined;function UTF8ArrayToString(heapOrArray,idx,maxBytesToRead){var endIdx=idx+maxBytesToRead;var endPtr=idx;while(heapOrArray[endPtr]&&!(endPtr>=endIdx))++endPtr;if(endPtr-idx>16&&heapOrArray.buffer&&UTF8Decoder){return UTF8Decoder.decode(heapOrArray.subarray(idx,endPtr))}var str="";while(idx<endPtr){var u0=heapOrArray[idx++];if(!(u0&128)){str+=String.fromCharCode(u0);continue}var u1=heapOrArray[idx++]&63;if((u0&224)==192){str+=String.fromCharCode((u0&31)<<6|u1);continue}var u2=heapOrArray[idx++]&63;if((u0&240)==224){u0=(u0&15)<<12|u1<<6|u2}else{u0=(u0&7)<<18|u1<<12|u2<<6|heapOrArray[idx++]&63}if(u0<65536){str+=String.fromCharCode(u0)}else{var ch=u0-65536;str+=String.fromCharCode(55296|ch>>10,56320|ch&1023)}}return str}function UTF8ToString(ptr,maxBytesToRead){return ptr?UTF8ArrayToString(HEAPU8,ptr,maxBytesToRead):""}function stringToUTF8Array(str,heap,outIdx,maxBytesToWrite){if(!(maxBytesToWrite>0))return 0;var startIdx=outIdx;var endIdx=outIdx+maxBytesToWrite-1;for(var i=0;i<str.length;++i){var u=str.charCodeAt(i);if(u>=55296&&u<=57343){var u1=str.charCodeAt(++i);u=65536+((u&1023)<<10)|u1&1023}if(u<=127){if(outIdx>=endIdx)break;heap[outIdx++]=u}else if(u<=2047){if(outIdx+1>=endIdx)break;heap[outIdx++]=192|u>>6;heap[outIdx++]=128|u&63}else if(u<=65535){if(outIdx+2>=endIdx)break;heap[outIdx++]=224|u>>12;heap[outIdx++]=128|u>>6&63;heap[outIdx++]=128|u&63}else{if(outIdx+3>=endIdx)break;heap[outIdx++]=240|u>>18;heap[outIdx++]=128|u>>12&63;heap[outIdx++]=128|u>>6&63;heap[outIdx++]=128|u&63}}heap[outIdx]=0;return outIdx-startIdx}function stringToUTF8(str,outPtr,maxBytesToWrite){return stringToUTF8Array(str,HEAPU8,outPtr,maxBytesToWrite)}var buffer,HEAP8,HEAPU8,HEAP16,HEAPU16,HEAP32,HEAPU32,HEAPF32,HEAPF64;function updateGlobalBufferAndViews(buf){buffer=buf;Module["HEAP8"]=HEAP8=new Int8Array(buf);Module["HEAP16"]=HEAP16=new Int16Array(buf);Module["HEAP32"]=HEAP32=new Int32Array(buf);Module["HEAPU8"]=HEAPU8=new Uint8Array(buf);Module["HEAPU16"]=HEAPU16=new Uint16Array(buf);Module["HEAPU32"]=HEAPU32=new Uint32Array(buf);Module["HEAPF32"]=HEAPF32=new Float32Array(buf);Module["HEAPF64"]=HEAPF64=new Float64Array(buf)}var INITIAL_MEMORY=Module["INITIAL_MEMORY"]||268435456;var wasmTable;var __ATPRERUN__=[];var __ATINIT__=[];var __ATPOSTRUN__=[];var runtimeInitialized=false;function preRun(){if(Module["preRun"]){if(typeof Module["preRun"]=="function")Module["preRun"]=[Module["preRun"]];while(Module["preRun"].length){addOnPreRun(Module["preRun"].shift())}}callRuntimeCallbacks(__ATPRERUN__)}function initRuntime(){runtimeInitialized=true;callRuntimeCallbacks(__ATINIT__)}function postRun(){if(Module["postRun"]){if(typeof Module["postRun"]=="function")Module["postRun"]=[Module["postRun"]];while(Module["postRun"].length){addOnPostRun(Module["postRun"].shift())}}callRuntimeCallbacks(__ATPOSTRUN__)}function addOnPreRun(cb){__ATPRERUN__.unshift(cb)}function addOnInit(cb){__ATINIT__.unshift(cb)}function addOnPostRun(cb){__ATPOSTRUN__.unshift(cb)}var runDependencies=0;var runDependencyWatcher=null;var dependenciesFulfilled=null;function addRunDependency(id){runDependencies++;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}}function removeRunDependency(id){runDependencies--;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}if(runDependencies==0){if(runDependencyWatcher!==null){clearInterval(runDependencyWatcher);runDependencyWatcher=null}if(dependenciesFulfilled){var callback=dependenciesFulfilled;dependenciesFulfilled=null;callback()}}}function abort(what){if(Module["onAbort"]){Module["onAbort"](what)}what="Aborted("+what+")";err(what);ABORT=true;EXITSTATUS=1;what+=". Build with -sASSERTIONS for more info.";var e=new WebAssembly.RuntimeError(what);readyPromiseReject(e);throw e}var dataURIPrefix="data:application/octet-stream;base64,";function isDataURI(filename){return filename.startsWith(dataURIPrefix)}var wasmBinaryFile;wasmBinaryFile="EDA.wasm";if(!isDataURI(wasmBinaryFile)){wasmBinaryFile=locateFile(wasmBinaryFile)}function getBinary(file){try{if(file==wasmBinaryFile&&wasmBinary){return new Uint8Array(wasmBinary)}if(readBinary){return readBinary(file)}throw"both async and sync fetching of the wasm failed"}catch(err){abort(err)}}function getBinaryPromise(){if(!wasmBinary&&(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER)){if(typeof fetch=="function"){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){if(!response["ok"]){throw"failed to load wasm binary file at '"+wasmBinaryFile+"'"}return response["arrayBuffer"]()}).catch(function(){return getBinary(wasmBinaryFile)})}}return Promise.resolve().then(function(){return getBinary(wasmBinaryFile)})}function createWasm(){var info={"a":asmLibraryArg};function receiveInstance(instance,module){var exports=instance.exports;Module["asm"]=exports;wasmMemory=Module["asm"]["f"];updateGlobalBufferAndViews(wasmMemory.buffer);wasmTable=Module["asm"]["k"];addOnInit(Module["asm"]["g"]);removeRunDependency("wasm-instantiate")}addRunDependency("wasm-instantiate");function receiveInstantiationResult(result){receiveInstance(result["instance"])}function instantiateArrayBuffer(receiver){return getBinaryPromise().then(function(binary){return WebAssembly.instantiate(binary,info)}).then(function(instance){return instance}).then(receiver,function(reason){err("failed to asynchronously prepare wasm: "+reason);abort(reason)})}function instantiateAsync(){if(!wasmBinary&&typeof WebAssembly.instantiateStreaming=="function"&&!isDataURI(wasmBinaryFile)&&typeof fetch=="function"){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){var result=WebAssembly.instantiateStreaming(response,info);return result.then(receiveInstantiationResult,function(reason){err("wasm streaming compile failed: "+reason);err("falling back to ArrayBuffer instantiation");return instantiateArrayBuffer(receiveInstantiationResult)})})}else{return instantiateArrayBuffer(receiveInstantiationResult)}}if(Module["instantiateWasm"]){try{var exports=Module["instantiateWasm"](info,receiveInstance);return exports}catch(e){err("Module.instantiateWasm callback failed with error: "+e);readyPromiseReject(e)}}instantiateAsync().catch(readyPromiseReject);return{}}function callRuntimeCallbacks(callbacks){while(callbacks.length>0){callbacks.shift()(Module)}}function ___assert_fail(condition,filename,line,func){abort("Assertion failed: "+UTF8ToString(condition)+", at: "+[filename?UTF8ToString(filename):"unknown filename",line,func?UTF8ToString(func):"unknown function"])}function ExceptionInfo(excPtr){this.excPtr=excPtr;this.ptr=excPtr-24;this.set_type=function(type){HEAPU32[this.ptr+4>>2]=type};this.get_type=function(){return HEAPU32[this.ptr+4>>2]};this.set_destructor=function(destructor){HEAPU32[this.ptr+8>>2]=destructor};this.get_destructor=function(){return HEAPU32[this.ptr+8>>2]};this.set_refcount=function(refcount){HEAP32[this.ptr>>2]=refcount};this.set_caught=function(caught){caught=caught?1:0;HEAP8[this.ptr+12>>0]=caught};this.get_caught=function(){return HEAP8[this.ptr+12>>0]!=0};this.set_rethrown=function(rethrown){rethrown=rethrown?1:0;HEAP8[this.ptr+13>>0]=rethrown};this.get_rethrown=function(){return HEAP8[this.ptr+13>>0]!=0};this.init=function(type,destructor){this.set_adjusted_ptr(0);this.set_type(type);this.set_destructor(destructor);this.set_refcount(0);this.set_caught(false);this.set_rethrown(false)};this.add_ref=function(){var value=HEAP32[this.ptr>>2];HEAP32[this.ptr>>2]=value+1};this.release_ref=function(){var prev=HEAP32[this.ptr>>2];HEAP32[this.ptr>>2]=prev-1;return prev===1};this.set_adjusted_ptr=function(adjustedPtr){HEAPU32[this.ptr+16>>2]=adjustedPtr};this.get_adjusted_ptr=function(){return HEAPU32[this.ptr+16>>2]};this.get_exception_ptr=function(){var isPointer=___cxa_is_pointer_type(this.get_type());if(isPointer){return HEAPU32[this.excPtr>>2]}var adjusted=this.get_adjusted_ptr();if(adjusted!==0)return adjusted;return this.excPtr}}var exceptionLast=0;var uncaughtExceptionCount=0;function ___cxa_throw(ptr,type,destructor){var info=new ExceptionInfo(ptr);info.init(type,destructor);exceptionLast=ptr;uncaughtExceptionCount++;throw ptr}function _abort(){abort("")}function _emscripten_memcpy_big(dest,src,num){HEAPU8.copyWithin(dest,src,src+num)}function getHeapMax(){return 2147483648}function emscripten_realloc_buffer(size){try{wasmMemory.grow(size-buffer.byteLength+65535>>>16);updateGlobalBufferAndViews(wasmMemory.buffer);return 1}catch(e){}}function _emscripten_resize_heap(requestedSize){var oldSize=HEAPU8.length;requestedSize=requestedSize>>>0;var maxHeapSize=getHeapMax();if(requestedSize>maxHeapSize){return false}let alignUp=(x,multiple)=>x+(multiple-x%multiple)%multiple;for(var cutDown=1;cutDown<=4;cutDown*=2){var overGrownHeapSize=oldSize*(1+.2/cutDown);overGrownHeapSize=Math.min(overGrownHeapSize,requestedSize+100663296);var newSize=Math.min(maxHeapSize,alignUp(Math.max(requestedSize,overGrownHeapSize),65536));var replacement=emscripten_realloc_buffer(newSize);if(replacement){return true}}return false}function getCFunc(ident){var func=Module["_"+ident];return func}function writeArrayToMemory(array,buffer){HEAP8.set(array,buffer)}function ccall(ident,returnType,argTypes,args,opts){var toC={"string":str=>{var ret=0;if(str!==null&&str!==undefined&&str!==0){var len=(str.length<<2)+1;ret=stackAlloc(len);stringToUTF8(str,ret,len)}return ret},"array":arr=>{var ret=stackAlloc(arr.length);writeArrayToMemory(arr,ret);return ret}};function convertReturnValue(ret){if(returnType==="string"){return UTF8ToString(ret)}if(returnType==="boolean")return Boolean(ret);return ret}var func=getCFunc(ident);var cArgs=[];var stack=0;if(args){for(var i=0;i<args.length;i++){var converter=toC[argTypes[i]];if(converter){if(stack===0)stack=stackSave();cArgs[i]=converter(args[i])}else{cArgs[i]=args[i]}}}var ret=func.apply(null,cArgs);function onDone(ret){if(stack!==0)stackRestore(stack);return convertReturnValue(ret)}ret=onDone(ret);return ret}function cwrap(ident,returnType,argTypes,opts){argTypes=argTypes||[];var numericArgs=argTypes.every(type=>type==="number"||type==="boolean");var numericRet=returnType!=="string";if(numericRet&&numericArgs&&!opts){return getCFunc(ident)}return function(){return ccall(ident,returnType,argTypes,arguments,opts)}}var asmLibraryArg={"a":___assert_fail,"b":___cxa_throw,"c":_abort,"e":_emscripten_memcpy_big,"d":_emscripten_resize_heap};var asm=createWasm();var ___wasm_call_ctors=Module["___wasm_call_ctors"]=function(){return(___wasm_call_ctors=Module["___wasm_call_ctors"]=Module["asm"]["g"]).apply(null,arguments)};var _principalComponentAnalysis=Module["_principalComponentAnalysis"]=function(){return(_principalComponentAnalysis=Module["_principalComponentAnalysis"]=Module["asm"]["h"]).apply(null,arguments)};var _error=Module["_error"]=function(){return(_error=Module["_error"]=Module["asm"]["i"]).apply(null,arguments)};var _principalComponentAnalysisNipals=Module["_principalComponentAnalysisNipals"]=function(){return(_principalComponentAnalysisNipals=Module["_principalComponentAnalysisNipals"]=Module["asm"]["j"]).apply(null,arguments)};var _free=Module["_free"]=function(){return(_free=Module["_free"]=Module["asm"]["l"]).apply(null,arguments)};var _malloc=Module["_malloc"]=function(){return(_malloc=Module["_malloc"]=Module["asm"]["m"]).apply(null,arguments)};var _partialLeastSquareRegression=Module["_partialLeastSquareRegression"]=function(){return(_partialLeastSquareRegression=Module["_partialLeastSquareRegression"]=Module["asm"]["n"]).apply(null,arguments)};var _generateDataset=Module["_generateDataset"]=function(){return(_generateDataset=Module["_generateDataset"]=Module["asm"]["o"]).apply(null,arguments)};var _normalizeDataset=Module["_normalizeDataset"]=function(){return(_normalizeDataset=Module["_normalizeDataset"]=Module["asm"]["p"]).apply(null,arguments)};var _trainLSSVM=Module["_trainLSSVM"]=function(){return(_trainLSSVM=Module["_trainLSSVM"]=Module["asm"]["q"]).apply(null,arguments)};var _predictByLSSVM=Module["_predictByLSSVM"]=function(){return(_predictByLSSVM=Module["_predictByLSSVM"]=Module["asm"]["r"]).apply(null,arguments)};var _trainAndAnalyzeLSSVM=Module["_trainAndAnalyzeLSSVM"]=function(){return(_trainAndAnalyzeLSSVM=Module["_trainAndAnalyzeLSSVM"]=Module["asm"]["s"]).apply(null,arguments)};var _fitLinearRegressionParamsWithDataNormalizing=Module["_fitLinearRegressionParamsWithDataNormalizing"]=function(){return(_fitLinearRegressionParamsWithDataNormalizing=Module["_fitLinearRegressionParamsWithDataNormalizing"]=Module["asm"]["t"]).apply(null,arguments)};var _fitLinearRegressionParams=Module["_fitLinearRegressionParams"]=function(){return(_fitLinearRegressionParams=Module["_fitLinearRegressionParams"]=Module["asm"]["u"]).apply(null,arguments)};var _fitSoftmax=Module["_fitSoftmax"]=function(){return(_fitSoftmax=Module["_fitSoftmax"]=Module["asm"]["v"]).apply(null,arguments)};var stackSave=Module["stackSave"]=function(){return(stackSave=Module["stackSave"]=Module["asm"]["w"]).apply(null,arguments)};var stackRestore=Module["stackRestore"]=function(){return(stackRestore=Module["stackRestore"]=Module["asm"]["x"]).apply(null,arguments)};var stackAlloc=Module["stackAlloc"]=function(){return(stackAlloc=Module["stackAlloc"]=Module["asm"]["y"]).apply(null,arguments)};var ___cxa_is_pointer_type=Module["___cxa_is_pointer_type"]=function(){return(___cxa_is_pointer_type=Module["___cxa_is_pointer_type"]=Module["asm"]["z"]).apply(null,arguments)};Module["ccall"]=ccall;Module["cwrap"]=cwrap;var calledRun;dependenciesFulfilled=function runCaller(){if(!calledRun)run();if(!calledRun)dependenciesFulfilled=runCaller};function run(args){args=args||arguments_;if(runDependencies>0){return}preRun();if(runDependencies>0){return}function doRun(){if(calledRun)return;calledRun=true;Module["calledRun"]=true;if(ABORT)return;initRuntime();readyPromiseResolve(Module);if(Module["onRuntimeInitialized"])Module["onRuntimeInitialized"]();postRun()}if(Module["setStatus"]){Module["setStatus"]("Running...");setTimeout(function(){setTimeout(function(){Module["setStatus"]("")},1);doRun()},1)}else{doRun()}}if(Module["preInit"]){if(typeof Module["preInit"]=="function")Module["preInit"]=[Module["preInit"]];while(Module["preInit"].length>0){Module["preInit"].pop()()}}run();


  return exportEDA.ready
}
);
})();
if (typeof exports === 'object' && typeof module === 'object')
  module.exports = exportEDA;
else if (typeof define === 'function' && define['amd'])
  define([], function() { return exportEDA; });
else if (typeof exports === 'object')
  exports["exportEDA"] = exportEDA;


var principalComponentAnalysis = {
  arguments: {
    columns: {
      type: 'floatColumns'
    },
    componentsCount: {
      type: 'num'
    },
    centerNum: {
      type: 'num'
    },
    scaleNum: {
      type: 'num'
    },
    components: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'columns',
        value: 'numOfRows'
      },
      numOfColumns: {
        ref: 'componentsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'tableFromColumns',
    source: 'components'
  }
}; // principalComponentAnalysis

var error = {
  arguments: {
    col1: {
      type: 'floatColumn'
    },
    col2: {
      type: 'floatColumn'
    }
  },
  output: {
    type: 'double',
    source: '_callResult'
  }
}; // error

var principalComponentAnalysisNipals = {
  arguments: {
    columns: {
      type: 'floatColumns'
    },
    componentsCount: {
      type: 'num'
    },
    components: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'columns',
        value: 'numOfRows'
      },
      numOfColumns: {
        ref: 'componentsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'tableFromColumns',
    source: 'components'
  }
}; // principalComponentAnalysisNipals

var partialLeastSquareRegression = {
  arguments: {
    features: {
      type: 'floatColumns'
    },
    predict: {
      type: 'floatColumn'
    },
    componentsCount: {
      type: 'num'
    },
    prediction: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'predict',
        value: 'numOfRows'
      }
    },
    regressionCoefficients: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'features',
        value: 'numOfColumns'
      }
    },
    tScores: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'predict',
        value: 'numOfRows'
      },
      numOfColumns: {
        ref: 'componentsCount',
        value: 'data'
      }
    },
    uScores: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'predict',
        value: 'numOfRows'
      },
      numOfColumns: {
        ref: 'componentsCount',
        value: 'data'
      }
    },
    xLoadings: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'features',
        value: 'numOfColumns'
      },
      numOfColumns: {
        ref: 'componentsCount',
        value: 'data'
      }
    },
    yLoadings: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'componentsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'objects',
    source: ['prediction', 'regressionCoefficients', 'tScores', 'uScores', 'xLoadings', 'yLoadings']
  }
}; // partialLeastSquareRegression

var generateDataset = {
  arguments: {
    kernel: {
      type: 'num'
    },
    kernelParams: {
      type: 'floatColumn'
    },
    samplesCount: {
      type: 'num'
    },
    featuresCount: {
      type: 'num'
    },
    min: {
      type: 'num'
    },
    max: {
      type: 'num'
    },
    violatorsPercentage: {
      type: 'num'
    },
    dataset: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'samplesCount',
        value: 'data'
      },
      numOfColumns: {
        ref: 'featuresCount',
        value: 'data'
      }
    },
    labels: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'samplesCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'objects',
    source: ['dataset', 'labels']
  }
}; // generateDataset

var normalizeDataset = {
  arguments: {
    data: {
      type: 'floatColumns'
    },
    normalizedData: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'data',
        value: 'numOfColumns'
      },
      numOfColumns: {
        ref: 'data',
        value: 'numOfRows'
      }
    },
    means: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'data',
        value: 'numOfColumns'
      }
    },
    stdDevs: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'data',
        value: 'numOfColumns'
      }
    }
  },
  output: {
    type: 'objects',
    source: ['normalizedData', 'means', 'stdDevs']
  }
}; // normalizeDataset

var trainLSSVM = {
  arguments: {
    gamma: {
      type: 'num'
    },
    kernel: {
      type: 'num'
    },
    kernelParams: {
      type: 'floatColumn'
    },
    modelParamsCount: {
      type: 'num'
    },
    precomputedWeightsCount: {
      type: 'num'
    },
    dataset: {
      type: 'floatColumns'
    },
    labels: {
      type: 'floatColumn'
    },
    normalizedData: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      },
      numOfColumns: {
        ref: 'dataset',
        value: 'numOfRows'
      }
    },
    means: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      }
    },
    stdDevs: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      }
    },
    modelParams: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'modelParamsCount',
        value: 'data'
      }
    },
    precomputedWeights: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'precomputedWeightsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'objects',
    source: ['normalizedData', 'means', 'stdDevs', 'modelParams', 'precomputedWeights']
  }
}; // trainLSSVM

var predictByLSSVM = {
  arguments: {
    kernel: {
      type: 'num'
    },
    kernelParams: {
      type: 'floatColumn'
    },
    normalizedData: {
      type: 'floatColumns'
    },
    labels: {
      type: 'floatColumn'
    },
    means: {
      type: 'floatColumn'
    },
    stdDevs: {
      type: 'floatColumn'
    },
    modelParams: {
      type: 'floatColumn'
    },
    precomputedWeights: {
      type: 'floatColumn'
    },
    targetData: {
      type: 'floatColumns'
    },
    prediction: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'targetData',
        value: 'numOfRows'
      }
    }
  },
  output: {
    type: 'column',
    source: 'prediction'
  }
}; // predictByLSSVM

var trainAndAnalyzeLSSVM = {
  arguments: {
    gamma: {
      type: 'num'
    },
    kernel: {
      type: 'num'
    },
    kernelParams: {
      type: 'floatColumn'
    },
    modelParamsCount: {
      type: 'num'
    },
    precomputedWeightsCount: {
      type: 'num'
    },
    confusionMatrixElementsCount: {
      type: 'num'
    },
    dataset: {
      type: 'floatColumns'
    },
    labels: {
      type: 'floatColumn'
    },
    normalizedData: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      },
      numOfColumns: {
        ref: 'dataset',
        value: 'numOfRows'
      }
    },
    means: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      }
    },
    stdDevs: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfColumns'
      }
    },
    modelParams: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'modelParamsCount',
        value: 'data'
      }
    },
    precomputedWeights: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'precomputedWeightsCount',
        value: 'data'
      }
    },
    predictedLabels: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfRows'
      }
    },
    correctness: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'dataset',
        value: 'numOfRows'
      }
    },
    consfusionMatrix: {
      type: 'newIntColumn',
      numOfRows: {
        ref: 'confusionMatrixElementsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'objects',
    source: ['normalizedData', 'means', 'stdDevs', 'modelParams', 'precomputedWeights', 'predictedLabels', 'correctness', 'consfusionMatrix']
  }
}; // trainAndAnalyzeLSSVM

var fitLinearRegressionParamsWithDataNormalizing = {
  arguments: {
    features: {
      type: 'floatColumns'
    },
    featureAvgs: {
      type: 'floatColumn'
    },
    featureStdDevs: {
      type: 'floatColumn'
    },
    targets: {
      type: 'floatColumn'
    },
    targetsAvg: {
      type: 'num'
    },
    targetsStdDev: {
      type: 'num'
    },
    paramsCount: {
      type: 'num'
    },
    params: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'paramsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'column',
    source: 'params'
  }
}; // fitLinearRegressionParamsWithDataNormalizing

var fitLinearRegressionParams = {
  arguments: {
    features: {
      type: 'floatColumns'
    },
    targets: {
      type: 'floatColumn'
    },
    paramsCount: {
      type: 'num'
    },
    params: {
      type: 'newFloatColumn',
      numOfRows: {
        ref: 'paramsCount',
        value: 'data'
      }
    }
  },
  output: {
    type: 'column',
    source: 'params'
  }
}; // fitLinearRegressionParams

var fitSoftmax = {
  arguments: {
    features: {
      type: 'floatColumns'
    },
    featureAvgs: {
      type: 'floatColumn'
    },
    featureStdDevs: {
      type: 'floatColumn'
    },
    targets: {
      type: 'intColumn'
    },
    classesCount: {
      type: 'num'
    },
    iterCount: {
      type: 'num'
    },
    learningRate: {
      type: 'num'
    },
    penalty: {
      type: 'num'
    },
    tolerance: {
      type: 'num'
    },
    paramsRows: {
      type: 'num'
    },
    paramsCols: {
      type: 'num'
    },
    params: {
      type: 'newFloatColumns',
      numOfRows: {
        ref: 'paramsRows',
        value: 'data'
      },
      numOfColumns: {
        ref: 'paramsCols',
        value: 'data'
      }
    }
  },
  output: {
    type: 'tableFromColumns',
    source: 'params'
  }
}; // fitSoftmax

var EDA = undefined;

async function initEDA() {
  if (EDA === undefined) {
    console.log("Wasm not Loaded, Loading");
    EDA = await exportEDA();
    EDA.principalComponentAnalysis = principalComponentAnalysis;
    EDA.error = error;
    EDA.principalComponentAnalysisNipals = principalComponentAnalysisNipals;
    EDA.partialLeastSquareRegression = partialLeastSquareRegression;
    EDA.generateDataset = generateDataset;
    EDA.normalizeDataset = normalizeDataset;
    EDA.trainLSSVM = trainLSSVM;
    EDA.predictByLSSVM = predictByLSSVM;
    EDA.trainAndAnalyzeLSSVM = trainAndAnalyzeLSSVM;
    EDA.fitLinearRegressionParamsWithDataNormalizing = fitLinearRegressionParamsWithDataNormalizing;
    EDA.fitLinearRegressionParams = fitLinearRegressionParams;
    EDA.fitSoftmax = fitSoftmax;
  } else {
    console.log("Wasm Loaded, Passing");
  }
}