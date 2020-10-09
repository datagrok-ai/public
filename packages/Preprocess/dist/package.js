var noisefilters =
/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, { enumerable: true, get: getter });
/******/ 		}
/******/ 	};
/******/
/******/ 	// define __esModule on exports
/******/ 	__webpack_require__.r = function(exports) {
/******/ 		if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 			Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 		}
/******/ 		Object.defineProperty(exports, '__esModule', { value: true });
/******/ 	};
/******/
/******/ 	// create a fake namespace object
/******/ 	// mode & 1: value is a module id, require it
/******/ 	// mode & 2: merge all properties of value into the ns
/******/ 	// mode & 4: return value when already ns object
/******/ 	// mode & 8|1: behave like require
/******/ 	__webpack_require__.t = function(value, mode) {
/******/ 		if(mode & 1) value = __webpack_require__(value);
/******/ 		if(mode & 8) return value;
/******/ 		if((mode & 4) && typeof value === 'object' && value && value.__esModule) return value;
/******/ 		var ns = Object.create(null);
/******/ 		__webpack_require__.r(ns);
/******/ 		Object.defineProperty(ns, 'default', { enumerable: true, value: value });
/******/ 		if(mode & 2 && typeof value != 'string') for(var key in value) __webpack_require__.d(ns, key, function(key) { return value[key]; }.bind(null, key));
/******/ 		return ns;
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = "./src/package.js");
/******/ })
/************************************************************************/
/******/ ({

/***/ "./src/package.js":
/*!************************!*\
  !*** ./src/package.js ***!
  \************************/
/*! exports provided: _package, SMA_filter, Exp_filter, Kalman_filter, MinMax_transform, Zscore_transform, box_cox_transform, fourier_filter */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "_package", function() { return _package; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "SMA_filter", function() { return SMA_filter; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "Exp_filter", function() { return Exp_filter; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "Kalman_filter", function() { return Kalman_filter; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "MinMax_transform", function() { return MinMax_transform; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "Zscore_transform", function() { return Zscore_transform; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "box_cox_transform", function() { return box_cox_transform; });
/* harmony export (binding) */ __webpack_require__.d(__webpack_exports__, "fourier_filter", function() { return fourier_filter; });
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/grok */ "datagrok-api/grok");
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! datagrok-api/ui */ "datagrok-api/ui");
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__);
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */




let _package = new datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Package"]();
fetch('algs.wasm').then(response => response.arrayBuffer())
.then(bytes => instantiate(bytes, importObject))
.then(instance => instance.exports.e());


//name: Moving Avarage Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: int window_size [SMA Window Size]
function SMA_filter(data,column_to_filter,window_size) {
    let js_wrapped_sma = Module.cwrap("sma", "null", ["number","number","number"]);
    let column_name = column_to_filter.name +' SMA Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length,window_size);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Exponential Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double filter_ratio [Exponential Filter Parameter]
function Exp_filter(data,column_to_filter,filter_ratio) {
    let js_wrapped_exp = Module.cwrap("exps", "null", ["number","number","number"]);
    let column_name = column_to_filter.name +' Exponentially Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    console.log(dataPtr);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_exp(dataPtr,filter_array.length,filter_ratio);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}


//name: Kalman Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double Q [Covariance of the process noise]
//input: double R [Covariance of the observation noise]
//input: double P [a posteriori estimate covariance]
function Kalman_filter(data,column_to_filter,Q,R,P) {
    let js_wrapped_kalm = Module.cwrap("kalman", "null", ["number","number","number","number","number"]);
    let column_name = column_to_filter.name +' Kalman Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_kalm(dataPtr,filter_array.length,Q,R,P)
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) { column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Min Max Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
function MinMax_transform(data,column_to_filter) {
    let js_wrapped_sma = Module.cwrap("minmax", "null", ["number","number"]);
    let column_name = column_to_filter.name +' Min Max Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Z-score Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
function Zscore_transform(data,column_to_filter) {
    let js_wrapped_sma = Module.cwrap("zscore", "null", ["number","number"]);
    let column_name = column_to_filter.name +' Z-score Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Box Cox Transform
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lambda
//input: double ofset
function box_cox_transform(data,column_to_filter,lambda, ofset) {
    let js_wrapped_sma = Module.cwrap("boxcox", "null", ["number","number","number","number"]);
    let column_name = column_to_filter.name +' Box Cox Transformed';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length,lambda,ofset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Fourier Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lowcut
//input: double hicut
function fourier_filter(data,column_to_filter,lowcut, hicut) {
    let js_wrapped_ff = Module.cwrap("ffilter", "null", ["number","number","number","number"]);
    let column_name = column_to_filter.name +' Fourier Filtered (L: '+lowcut + '; H: ' + hicut + ')';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_ff(dataPtr,filter_array.length,lowcut,hicut);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__["Column"].fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

/***/ }),

/***/ "datagrok-api/dg":
/*!*********************!*\
  !*** external "DG" ***!
  \*********************/
/*! no static exports found */
/***/ (function(module, exports) {

module.exports = DG;

/***/ }),

/***/ "datagrok-api/grok":
/*!***********************!*\
  !*** external "grok" ***!
  \***********************/
/*! no static exports found */
/***/ (function(module, exports) {

module.exports = grok;

/***/ }),

/***/ "datagrok-api/ui":
/*!*********************!*\
  !*** external "ui" ***!
  \*********************/
/*! no static exports found */
/***/ (function(module, exports) {

module.exports = ui;

/***/ })

/******/ });
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIndlYnBhY2s6Ly9ub2lzZWZpbHRlcnMvd2VicGFjay9ib290c3RyYXAiLCJ3ZWJwYWNrOi8vbm9pc2VmaWx0ZXJzLy4vc3JjL3BhY2thZ2UuanMiLCJ3ZWJwYWNrOi8vbm9pc2VmaWx0ZXJzL2V4dGVybmFsIFwiREdcIiIsIndlYnBhY2s6Ly9ub2lzZWZpbHRlcnMvZXh0ZXJuYWwgXCJncm9rXCIiLCJ3ZWJwYWNrOi8vbm9pc2VmaWx0ZXJzL2V4dGVybmFsIFwidWlcIiJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOztRQUFBO1FBQ0E7O1FBRUE7UUFDQTs7UUFFQTtRQUNBO1FBQ0E7UUFDQTtRQUNBO1FBQ0E7UUFDQTtRQUNBO1FBQ0E7UUFDQTs7UUFFQTtRQUNBOztRQUVBO1FBQ0E7O1FBRUE7UUFDQTtRQUNBOzs7UUFHQTtRQUNBOztRQUVBO1FBQ0E7O1FBRUE7UUFDQTtRQUNBO1FBQ0EsMENBQTBDLGdDQUFnQztRQUMxRTtRQUNBOztRQUVBO1FBQ0E7UUFDQTtRQUNBLHdEQUF3RCxrQkFBa0I7UUFDMUU7UUFDQSxpREFBaUQsY0FBYztRQUMvRDs7UUFFQTtRQUNBO1FBQ0E7UUFDQTtRQUNBO1FBQ0E7UUFDQTtRQUNBO1FBQ0E7UUFDQTtRQUNBO1FBQ0EseUNBQXlDLGlDQUFpQztRQUMxRSxnSEFBZ0gsbUJBQW1CLEVBQUU7UUFDckk7UUFDQTs7UUFFQTtRQUNBO1FBQ0E7UUFDQSwyQkFBMkIsMEJBQTBCLEVBQUU7UUFDdkQsaUNBQWlDLGVBQWU7UUFDaEQ7UUFDQTtRQUNBOztRQUVBO1FBQ0Esc0RBQXNELCtEQUErRDs7UUFFckg7UUFDQTs7O1FBR0E7UUFDQTs7Ozs7Ozs7Ozs7OztBQ2xGQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUMwQztBQUNKO0FBQ0E7O0FBRS9CLG1CQUFtQix1REFBVTtBQUNwQztBQUNBO0FBQ0E7OztBQUdBO0FBQ0E7QUFDQTtBQUNBO0FBQ087QUFDUDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZUFBZSx5QkFBeUIsT0FBTztBQUMvQyxxQkFBcUIsc0RBQVM7QUFDOUI7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZUFBZSx5QkFBeUIsT0FBTztBQUMvQyxxQkFBcUIsc0RBQVM7QUFDOUI7QUFDQTs7O0FBR0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ087QUFDUDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZUFBZSx5QkFBeUIsT0FBTztBQUMvQyxxQkFBcUIsc0RBQVM7QUFDOUI7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxlQUFlLHlCQUF5QixPQUFPO0FBQy9DLHFCQUFxQixzREFBUztBQUM5QjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGVBQWUseUJBQXlCLE9BQU87QUFDL0MscUJBQXFCLHNEQUFTO0FBQzlCO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGVBQWUseUJBQXlCLE9BQU87QUFDL0MscUJBQXFCLHNEQUFTO0FBQzlCO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNPO0FBQ1A7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZUFBZSx5QkFBeUIsT0FBTztBQUMvQyxxQkFBcUIsc0RBQVM7QUFDOUI7QUFDQSxDOzs7Ozs7Ozs7OztBQ2hLQSxvQjs7Ozs7Ozs7Ozs7QUNBQSxzQjs7Ozs7Ozs7Ozs7QUNBQSxvQiIsImZpbGUiOiJwYWNrYWdlLmpzIiwic291cmNlc0NvbnRlbnQiOlsiIFx0Ly8gVGhlIG1vZHVsZSBjYWNoZVxuIFx0dmFyIGluc3RhbGxlZE1vZHVsZXMgPSB7fTtcblxuIFx0Ly8gVGhlIHJlcXVpcmUgZnVuY3Rpb25cbiBcdGZ1bmN0aW9uIF9fd2VicGFja19yZXF1aXJlX18obW9kdWxlSWQpIHtcblxuIFx0XHQvLyBDaGVjayBpZiBtb2R1bGUgaXMgaW4gY2FjaGVcbiBcdFx0aWYoaW5zdGFsbGVkTW9kdWxlc1ttb2R1bGVJZF0pIHtcbiBcdFx0XHRyZXR1cm4gaW5zdGFsbGVkTW9kdWxlc1ttb2R1bGVJZF0uZXhwb3J0cztcbiBcdFx0fVxuIFx0XHQvLyBDcmVhdGUgYSBuZXcgbW9kdWxlIChhbmQgcHV0IGl0IGludG8gdGhlIGNhY2hlKVxuIFx0XHR2YXIgbW9kdWxlID0gaW5zdGFsbGVkTW9kdWxlc1ttb2R1bGVJZF0gPSB7XG4gXHRcdFx0aTogbW9kdWxlSWQsXG4gXHRcdFx0bDogZmFsc2UsXG4gXHRcdFx0ZXhwb3J0czoge31cbiBcdFx0fTtcblxuIFx0XHQvLyBFeGVjdXRlIHRoZSBtb2R1bGUgZnVuY3Rpb25cbiBcdFx0bW9kdWxlc1ttb2R1bGVJZF0uY2FsbChtb2R1bGUuZXhwb3J0cywgbW9kdWxlLCBtb2R1bGUuZXhwb3J0cywgX193ZWJwYWNrX3JlcXVpcmVfXyk7XG5cbiBcdFx0Ly8gRmxhZyB0aGUgbW9kdWxlIGFzIGxvYWRlZFxuIFx0XHRtb2R1bGUubCA9IHRydWU7XG5cbiBcdFx0Ly8gUmV0dXJuIHRoZSBleHBvcnRzIG9mIHRoZSBtb2R1bGVcbiBcdFx0cmV0dXJuIG1vZHVsZS5leHBvcnRzO1xuIFx0fVxuXG5cbiBcdC8vIGV4cG9zZSB0aGUgbW9kdWxlcyBvYmplY3QgKF9fd2VicGFja19tb2R1bGVzX18pXG4gXHRfX3dlYnBhY2tfcmVxdWlyZV9fLm0gPSBtb2R1bGVzO1xuXG4gXHQvLyBleHBvc2UgdGhlIG1vZHVsZSBjYWNoZVxuIFx0X193ZWJwYWNrX3JlcXVpcmVfXy5jID0gaW5zdGFsbGVkTW9kdWxlcztcblxuIFx0Ly8gZGVmaW5lIGdldHRlciBmdW5jdGlvbiBmb3IgaGFybW9ueSBleHBvcnRzXG4gXHRfX3dlYnBhY2tfcmVxdWlyZV9fLmQgPSBmdW5jdGlvbihleHBvcnRzLCBuYW1lLCBnZXR0ZXIpIHtcbiBcdFx0aWYoIV9fd2VicGFja19yZXF1aXJlX18ubyhleHBvcnRzLCBuYW1lKSkge1xuIFx0XHRcdE9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBuYW1lLCB7IGVudW1lcmFibGU6IHRydWUsIGdldDogZ2V0dGVyIH0pO1xuIFx0XHR9XG4gXHR9O1xuXG4gXHQvLyBkZWZpbmUgX19lc01vZHVsZSBvbiBleHBvcnRzXG4gXHRfX3dlYnBhY2tfcmVxdWlyZV9fLnIgPSBmdW5jdGlvbihleHBvcnRzKSB7XG4gXHRcdGlmKHR5cGVvZiBTeW1ib2wgIT09ICd1bmRlZmluZWQnICYmIFN5bWJvbC50b1N0cmluZ1RhZykge1xuIFx0XHRcdE9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBTeW1ib2wudG9TdHJpbmdUYWcsIHsgdmFsdWU6ICdNb2R1bGUnIH0pO1xuIFx0XHR9XG4gXHRcdE9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCAnX19lc01vZHVsZScsIHsgdmFsdWU6IHRydWUgfSk7XG4gXHR9O1xuXG4gXHQvLyBjcmVhdGUgYSBmYWtlIG5hbWVzcGFjZSBvYmplY3RcbiBcdC8vIG1vZGUgJiAxOiB2YWx1ZSBpcyBhIG1vZHVsZSBpZCwgcmVxdWlyZSBpdFxuIFx0Ly8gbW9kZSAmIDI6IG1lcmdlIGFsbCBwcm9wZXJ0aWVzIG9mIHZhbHVlIGludG8gdGhlIG5zXG4gXHQvLyBtb2RlICYgNDogcmV0dXJuIHZhbHVlIHdoZW4gYWxyZWFkeSBucyBvYmplY3RcbiBcdC8vIG1vZGUgJiA4fDE6IGJlaGF2ZSBsaWtlIHJlcXVpcmVcbiBcdF9fd2VicGFja19yZXF1aXJlX18udCA9IGZ1bmN0aW9uKHZhbHVlLCBtb2RlKSB7XG4gXHRcdGlmKG1vZGUgJiAxKSB2YWx1ZSA9IF9fd2VicGFja19yZXF1aXJlX18odmFsdWUpO1xuIFx0XHRpZihtb2RlICYgOCkgcmV0dXJuIHZhbHVlO1xuIFx0XHRpZigobW9kZSAmIDQpICYmIHR5cGVvZiB2YWx1ZSA9PT0gJ29iamVjdCcgJiYgdmFsdWUgJiYgdmFsdWUuX19lc01vZHVsZSkgcmV0dXJuIHZhbHVlO1xuIFx0XHR2YXIgbnMgPSBPYmplY3QuY3JlYXRlKG51bGwpO1xuIFx0XHRfX3dlYnBhY2tfcmVxdWlyZV9fLnIobnMpO1xuIFx0XHRPYmplY3QuZGVmaW5lUHJvcGVydHkobnMsICdkZWZhdWx0JywgeyBlbnVtZXJhYmxlOiB0cnVlLCB2YWx1ZTogdmFsdWUgfSk7XG4gXHRcdGlmKG1vZGUgJiAyICYmIHR5cGVvZiB2YWx1ZSAhPSAnc3RyaW5nJykgZm9yKHZhciBrZXkgaW4gdmFsdWUpIF9fd2VicGFja19yZXF1aXJlX18uZChucywga2V5LCBmdW5jdGlvbihrZXkpIHsgcmV0dXJuIHZhbHVlW2tleV07IH0uYmluZChudWxsLCBrZXkpKTtcbiBcdFx0cmV0dXJuIG5zO1xuIFx0fTtcblxuIFx0Ly8gZ2V0RGVmYXVsdEV4cG9ydCBmdW5jdGlvbiBmb3IgY29tcGF0aWJpbGl0eSB3aXRoIG5vbi1oYXJtb255IG1vZHVsZXNcbiBcdF9fd2VicGFja19yZXF1aXJlX18ubiA9IGZ1bmN0aW9uKG1vZHVsZSkge1xuIFx0XHR2YXIgZ2V0dGVyID0gbW9kdWxlICYmIG1vZHVsZS5fX2VzTW9kdWxlID9cbiBcdFx0XHRmdW5jdGlvbiBnZXREZWZhdWx0KCkgeyByZXR1cm4gbW9kdWxlWydkZWZhdWx0J107IH0gOlxuIFx0XHRcdGZ1bmN0aW9uIGdldE1vZHVsZUV4cG9ydHMoKSB7IHJldHVybiBtb2R1bGU7IH07XG4gXHRcdF9fd2VicGFja19yZXF1aXJlX18uZChnZXR0ZXIsICdhJywgZ2V0dGVyKTtcbiBcdFx0cmV0dXJuIGdldHRlcjtcbiBcdH07XG5cbiBcdC8vIE9iamVjdC5wcm90b3R5cGUuaGFzT3duUHJvcGVydHkuY2FsbFxuIFx0X193ZWJwYWNrX3JlcXVpcmVfXy5vID0gZnVuY3Rpb24ob2JqZWN0LCBwcm9wZXJ0eSkgeyByZXR1cm4gT2JqZWN0LnByb3RvdHlwZS5oYXNPd25Qcm9wZXJ0eS5jYWxsKG9iamVjdCwgcHJvcGVydHkpOyB9O1xuXG4gXHQvLyBfX3dlYnBhY2tfcHVibGljX3BhdGhfX1xuIFx0X193ZWJwYWNrX3JlcXVpcmVfXy5wID0gXCJcIjtcblxuXG4gXHQvLyBMb2FkIGVudHJ5IG1vZHVsZSBhbmQgcmV0dXJuIGV4cG9ydHNcbiBcdHJldHVybiBfX3dlYnBhY2tfcmVxdWlyZV9fKF9fd2VicGFja19yZXF1aXJlX18ucyA9IFwiLi9zcmMvcGFja2FnZS5qc1wiKTtcbiIsIi8qIERvIG5vdCBjaGFuZ2UgdGhlc2UgaW1wb3J0IGxpbmVzLiBEYXRhZ3JvayB3aWxsIGltcG9ydCBBUEkgbGlicmFyeSBpbiBleGFjdGx5IHRoZSBzYW1lIG1hbm5lciAqL1xuaW1wb3J0ICogYXMgZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XG5pbXBvcnQgKiBhcyB1aSBmcm9tICdkYXRhZ3Jvay1hcGkvdWknO1xuaW1wb3J0ICogYXMgREcgZnJvbSBcImRhdGFncm9rLWFwaS9kZ1wiO1xuXG5leHBvcnQgbGV0IF9wYWNrYWdlID0gbmV3IERHLlBhY2thZ2UoKTtcbmZldGNoKCdhbGdzLndhc20nKS50aGVuKHJlc3BvbnNlID0+IHJlc3BvbnNlLmFycmF5QnVmZmVyKCkpXG4udGhlbihieXRlcyA9PiBpbnN0YW50aWF0ZShieXRlcywgaW1wb3J0T2JqZWN0KSlcbi50aGVuKGluc3RhbmNlID0+IGluc3RhbmNlLmV4cG9ydHMuZSgpKTtcblxuXG4vL25hbWU6IE1vdmluZyBBdmFyYWdlIEZpbHRlclxuLy9pbnB1dDogZGF0YWZyYW1lIGRhdGEgW0lucHV0IGRhdGEgdGFibGVdXG4vL2lucHV0OiBjb2x1bW4gY29sdW1uX3RvX2ZpbHRlciBbU2lnbmFsIHRvIEZpbHRlcl1cbi8vaW5wdXQ6IGludCB3aW5kb3dfc2l6ZSBbU01BIFdpbmRvdyBTaXplXVxuZXhwb3J0IGZ1bmN0aW9uIFNNQV9maWx0ZXIoZGF0YSxjb2x1bW5fdG9fZmlsdGVyLHdpbmRvd19zaXplKSB7XG4gICAgbGV0IGpzX3dyYXBwZWRfc21hID0gTW9kdWxlLmN3cmFwKFwic21hXCIsIFwibnVsbFwiLCBbXCJudW1iZXJcIixcIm51bWJlclwiLFwibnVtYmVyXCJdKTtcbiAgICBsZXQgY29sdW1uX25hbWUgPSBjb2x1bW5fdG9fZmlsdGVyLm5hbWUgKycgU01BIEZpbHRlcmVkJztcbiAgICBsZXQgZmlsdGVyX2FycmF5ID0gY29sdW1uX3RvX2ZpbHRlci5nZXRSYXdEYXRhKCk7XG4gICAgbGV0IG5EYXRhQnl0ZXMgPSBmaWx0ZXJfYXJyYXkubGVuZ3RoICogZmlsdGVyX2FycmF5LkJZVEVTX1BFUl9FTEVNRU5UO1xuICAgIGxldCBkYXRhUHRyID0gTW9kdWxlLl9tYWxsb2MobkRhdGFCeXRlcyk7XG4gICAgbGV0IGRhdGFIZWFwID0gbmV3IFVpbnQ4QXJyYXkoTW9kdWxlLkhFQVBVOC5idWZmZXIsIGRhdGFQdHIsIG5EYXRhQnl0ZXMpO1xuICAgIGRhdGFIZWFwLnNldChuZXcgVWludDhBcnJheShmaWx0ZXJfYXJyYXkuYnVmZmVyKSk7XG4gICAganNfd3JhcHBlZF9zbWEoZGF0YVB0cixmaWx0ZXJfYXJyYXkubGVuZ3RoLHdpbmRvd19zaXplKTtcbiAgICBsZXQgcmVzdWx0ID0gbmV3IEZsb2F0MzJBcnJheShkYXRhSGVhcC5idWZmZXIsIGRhdGFIZWFwLmJ5dGVPZmZzZXQsIGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBjb2x1bW4gPSBuZXcgRmxvYXQzMkFycmF5KGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBpPTA7XG4gICAgZm9yIChpID0gMDsgaSA8IGZpbHRlcl9hcnJheS5sZW5ndGg7IGkrKykge2NvbHVtbltpXT1yZXN1bHRbaV07fVxuICAgIGRhdGEuY29sdW1ucy5hZGQoREcuQ29sdW1uLmZyb21GbG9hdDMyQXJyYXkoY29sdW1uX25hbWUsIGNvbHVtbikpO1xuICAgIE1vZHVsZS5fZnJlZShkYXRhSGVhcC5ieXRlT2Zmc2V0KTtcbn1cblxuLy9uYW1lOiBFeHBvbmVudGlhbCBGaWx0ZXJcbi8vaW5wdXQ6IGRhdGFmcmFtZSBkYXRhIFtJbnB1dCBkYXRhIHRhYmxlXVxuLy9pbnB1dDogY29sdW1uIGNvbHVtbl90b19maWx0ZXIgW1NpZ25hbCB0byBGaWx0ZXJdXG4vL2lucHV0OiBkb3VibGUgZmlsdGVyX3JhdGlvIFtFeHBvbmVudGlhbCBGaWx0ZXIgUGFyYW1ldGVyXVxuZXhwb3J0IGZ1bmN0aW9uIEV4cF9maWx0ZXIoZGF0YSxjb2x1bW5fdG9fZmlsdGVyLGZpbHRlcl9yYXRpbykge1xuICAgIGxldCBqc193cmFwcGVkX2V4cCA9IE1vZHVsZS5jd3JhcChcImV4cHNcIiwgXCJudWxsXCIsIFtcIm51bWJlclwiLFwibnVtYmVyXCIsXCJudW1iZXJcIl0pO1xuICAgIGxldCBjb2x1bW5fbmFtZSA9IGNvbHVtbl90b19maWx0ZXIubmFtZSArJyBFeHBvbmVudGlhbGx5IEZpbHRlcmVkJztcbiAgICBsZXQgZmlsdGVyX2FycmF5ID0gY29sdW1uX3RvX2ZpbHRlci5nZXRSYXdEYXRhKCk7XG4gICAgbGV0IG5EYXRhQnl0ZXMgPSBmaWx0ZXJfYXJyYXkubGVuZ3RoICogZmlsdGVyX2FycmF5LkJZVEVTX1BFUl9FTEVNRU5UO1xuICAgIGxldCBkYXRhUHRyID0gTW9kdWxlLl9tYWxsb2MobkRhdGFCeXRlcyk7XG4gICAgY29uc29sZS5sb2coZGF0YVB0cik7XG4gICAgbGV0IGRhdGFIZWFwID0gbmV3IFVpbnQ4QXJyYXkoTW9kdWxlLkhFQVBVOC5idWZmZXIsIGRhdGFQdHIsIG5EYXRhQnl0ZXMpO1xuICAgIGRhdGFIZWFwLnNldChuZXcgVWludDhBcnJheShmaWx0ZXJfYXJyYXkuYnVmZmVyKSk7XG4gICAganNfd3JhcHBlZF9leHAoZGF0YVB0cixmaWx0ZXJfYXJyYXkubGVuZ3RoLGZpbHRlcl9yYXRpbyk7XG4gICAgbGV0IHJlc3VsdCA9IG5ldyBGbG9hdDMyQXJyYXkoZGF0YUhlYXAuYnVmZmVyLCBkYXRhSGVhcC5ieXRlT2Zmc2V0LCBmaWx0ZXJfYXJyYXkubGVuZ3RoKTtcbiAgICBsZXQgY29sdW1uID0gbmV3IEZsb2F0MzJBcnJheShmaWx0ZXJfYXJyYXkubGVuZ3RoKTtcbiAgICBsZXQgaT0wO1xuICAgIGZvciAoaSA9IDA7IGkgPCBmaWx0ZXJfYXJyYXkubGVuZ3RoOyBpKyspIHtjb2x1bW5baV09cmVzdWx0W2ldO31cbiAgICBkYXRhLmNvbHVtbnMuYWRkKERHLkNvbHVtbi5mcm9tRmxvYXQzMkFycmF5KGNvbHVtbl9uYW1lLCBjb2x1bW4pKTtcbiAgICBNb2R1bGUuX2ZyZWUoZGF0YUhlYXAuYnl0ZU9mZnNldCk7XG59XG5cblxuLy9uYW1lOiBLYWxtYW4gRmlsdGVyXG4vL2lucHV0OiBkYXRhZnJhbWUgZGF0YSBbSW5wdXQgZGF0YSB0YWJsZV1cbi8vaW5wdXQ6IGNvbHVtbiBjb2x1bW5fdG9fZmlsdGVyIFtTaWduYWwgdG8gRmlsdGVyXVxuLy9pbnB1dDogZG91YmxlIFEgW0NvdmFyaWFuY2Ugb2YgdGhlIHByb2Nlc3Mgbm9pc2VdXG4vL2lucHV0OiBkb3VibGUgUiBbQ292YXJpYW5jZSBvZiB0aGUgb2JzZXJ2YXRpb24gbm9pc2VdXG4vL2lucHV0OiBkb3VibGUgUCBbYSBwb3N0ZXJpb3JpIGVzdGltYXRlIGNvdmFyaWFuY2VdXG5leHBvcnQgZnVuY3Rpb24gS2FsbWFuX2ZpbHRlcihkYXRhLGNvbHVtbl90b19maWx0ZXIsUSxSLFApIHtcbiAgICBsZXQganNfd3JhcHBlZF9rYWxtID0gTW9kdWxlLmN3cmFwKFwia2FsbWFuXCIsIFwibnVsbFwiLCBbXCJudW1iZXJcIixcIm51bWJlclwiLFwibnVtYmVyXCIsXCJudW1iZXJcIixcIm51bWJlclwiXSk7XG4gICAgbGV0IGNvbHVtbl9uYW1lID0gY29sdW1uX3RvX2ZpbHRlci5uYW1lICsnIEthbG1hbiBGaWx0ZXJlZCc7XG4gICAgbGV0IGZpbHRlcl9hcnJheSA9IGNvbHVtbl90b19maWx0ZXIuZ2V0UmF3RGF0YSgpO1xuICAgIGxldCBuRGF0YUJ5dGVzID0gZmlsdGVyX2FycmF5Lmxlbmd0aCAqIGZpbHRlcl9hcnJheS5CWVRFU19QRVJfRUxFTUVOVDtcbiAgICBsZXQgZGF0YVB0ciA9IE1vZHVsZS5fbWFsbG9jKG5EYXRhQnl0ZXMpO1xuICAgIGxldCBkYXRhSGVhcCA9IG5ldyBVaW50OEFycmF5KE1vZHVsZS5IRUFQVTguYnVmZmVyLCBkYXRhUHRyLCBuRGF0YUJ5dGVzKTtcbiAgICBkYXRhSGVhcC5zZXQobmV3IFVpbnQ4QXJyYXkoZmlsdGVyX2FycmF5LmJ1ZmZlcikpO1xuICAgIGpzX3dyYXBwZWRfa2FsbShkYXRhUHRyLGZpbHRlcl9hcnJheS5sZW5ndGgsUSxSLFApXG4gICAgbGV0IHJlc3VsdCA9IG5ldyBGbG9hdDMyQXJyYXkoZGF0YUhlYXAuYnVmZmVyLCBkYXRhSGVhcC5ieXRlT2Zmc2V0LCBmaWx0ZXJfYXJyYXkubGVuZ3RoKTtcbiAgICBsZXQgY29sdW1uID0gbmV3IEZsb2F0MzJBcnJheShmaWx0ZXJfYXJyYXkubGVuZ3RoKTtcbiAgICBsZXQgaT0wO1xuICAgIGZvciAoaSA9IDA7IGkgPCBmaWx0ZXJfYXJyYXkubGVuZ3RoOyBpKyspIHsgY29sdW1uW2ldPXJlc3VsdFtpXTt9XG4gICAgZGF0YS5jb2x1bW5zLmFkZChERy5Db2x1bW4uZnJvbUZsb2F0MzJBcnJheShjb2x1bW5fbmFtZSwgY29sdW1uKSk7XG4gICAgTW9kdWxlLl9mcmVlKGRhdGFIZWFwLmJ5dGVPZmZzZXQpO1xufVxuXG4vL25hbWU6IE1pbiBNYXggTm9ybWFsaXphdGlvblxuLy9pbnB1dDogZGF0YWZyYW1lIGRhdGEgW0lucHV0IGRhdGEgdGFibGVdXG4vL2lucHV0OiBjb2x1bW4gY29sdW1uX3RvX2ZpbHRlciBbU2lnbmFsIHRvIEZpbHRlcl1cbmV4cG9ydCBmdW5jdGlvbiBNaW5NYXhfdHJhbnNmb3JtKGRhdGEsY29sdW1uX3RvX2ZpbHRlcikge1xuICAgIGxldCBqc193cmFwcGVkX3NtYSA9IE1vZHVsZS5jd3JhcChcIm1pbm1heFwiLCBcIm51bGxcIiwgW1wibnVtYmVyXCIsXCJudW1iZXJcIl0pO1xuICAgIGxldCBjb2x1bW5fbmFtZSA9IGNvbHVtbl90b19maWx0ZXIubmFtZSArJyBNaW4gTWF4IE5vcm1hbGl6ZWQnO1xuICAgIGxldCBmaWx0ZXJfYXJyYXkgPSBjb2x1bW5fdG9fZmlsdGVyLmdldFJhd0RhdGEoKTtcbiAgICBsZXQgbkRhdGFCeXRlcyA9IGZpbHRlcl9hcnJheS5sZW5ndGggKiBmaWx0ZXJfYXJyYXkuQllURVNfUEVSX0VMRU1FTlQ7XG4gICAgbGV0IGRhdGFQdHIgPSBNb2R1bGUuX21hbGxvYyhuRGF0YUJ5dGVzKTtcbiAgICBsZXQgZGF0YUhlYXAgPSBuZXcgVWludDhBcnJheShNb2R1bGUuSEVBUFU4LmJ1ZmZlciwgZGF0YVB0ciwgbkRhdGFCeXRlcyk7XG4gICAgZGF0YUhlYXAuc2V0KG5ldyBVaW50OEFycmF5KGZpbHRlcl9hcnJheS5idWZmZXIpKTtcbiAgICBqc193cmFwcGVkX3NtYShkYXRhUHRyLGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCByZXN1bHQgPSBuZXcgRmxvYXQzMkFycmF5KGRhdGFIZWFwLmJ1ZmZlciwgZGF0YUhlYXAuYnl0ZU9mZnNldCwgZmlsdGVyX2FycmF5Lmxlbmd0aCk7XG4gICAgbGV0IGNvbHVtbiA9IG5ldyBGbG9hdDMyQXJyYXkoZmlsdGVyX2FycmF5Lmxlbmd0aCk7XG4gICAgbGV0IGk9MDtcbiAgICBmb3IgKGkgPSAwOyBpIDwgZmlsdGVyX2FycmF5Lmxlbmd0aDsgaSsrKSB7Y29sdW1uW2ldPXJlc3VsdFtpXTt9XG4gICAgZGF0YS5jb2x1bW5zLmFkZChERy5Db2x1bW4uZnJvbUZsb2F0MzJBcnJheShjb2x1bW5fbmFtZSwgY29sdW1uKSk7XG4gICAgTW9kdWxlLl9mcmVlKGRhdGFIZWFwLmJ5dGVPZmZzZXQpO1xufVxuXG4vL25hbWU6IFotc2NvcmUgTm9ybWFsaXphdGlvblxuLy9pbnB1dDogZGF0YWZyYW1lIGRhdGEgW0lucHV0IGRhdGEgdGFibGVdXG4vL2lucHV0OiBjb2x1bW4gY29sdW1uX3RvX2ZpbHRlciBbU2lnbmFsIHRvIEZpbHRlcl1cbmV4cG9ydCBmdW5jdGlvbiBac2NvcmVfdHJhbnNmb3JtKGRhdGEsY29sdW1uX3RvX2ZpbHRlcikge1xuICAgIGxldCBqc193cmFwcGVkX3NtYSA9IE1vZHVsZS5jd3JhcChcInpzY29yZVwiLCBcIm51bGxcIiwgW1wibnVtYmVyXCIsXCJudW1iZXJcIl0pO1xuICAgIGxldCBjb2x1bW5fbmFtZSA9IGNvbHVtbl90b19maWx0ZXIubmFtZSArJyBaLXNjb3JlIE5vcm1hbGl6ZWQnO1xuICAgIGxldCBmaWx0ZXJfYXJyYXkgPSBjb2x1bW5fdG9fZmlsdGVyLmdldFJhd0RhdGEoKTtcbiAgICBsZXQgbkRhdGFCeXRlcyA9IGZpbHRlcl9hcnJheS5sZW5ndGggKiBmaWx0ZXJfYXJyYXkuQllURVNfUEVSX0VMRU1FTlQ7XG4gICAgbGV0IGRhdGFQdHIgPSBNb2R1bGUuX21hbGxvYyhuRGF0YUJ5dGVzKTtcbiAgICBsZXQgZGF0YUhlYXAgPSBuZXcgVWludDhBcnJheShNb2R1bGUuSEVBUFU4LmJ1ZmZlciwgZGF0YVB0ciwgbkRhdGFCeXRlcyk7XG4gICAgZGF0YUhlYXAuc2V0KG5ldyBVaW50OEFycmF5KGZpbHRlcl9hcnJheS5idWZmZXIpKTtcbiAgICBqc193cmFwcGVkX3NtYShkYXRhUHRyLGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCByZXN1bHQgPSBuZXcgRmxvYXQzMkFycmF5KGRhdGFIZWFwLmJ1ZmZlciwgZGF0YUhlYXAuYnl0ZU9mZnNldCwgZmlsdGVyX2FycmF5Lmxlbmd0aCk7XG4gICAgbGV0IGNvbHVtbiA9IG5ldyBGbG9hdDMyQXJyYXkoZmlsdGVyX2FycmF5Lmxlbmd0aCk7XG4gICAgbGV0IGk9MDtcbiAgICBmb3IgKGkgPSAwOyBpIDwgZmlsdGVyX2FycmF5Lmxlbmd0aDsgaSsrKSB7Y29sdW1uW2ldPXJlc3VsdFtpXTt9XG4gICAgZGF0YS5jb2x1bW5zLmFkZChERy5Db2x1bW4uZnJvbUZsb2F0MzJBcnJheShjb2x1bW5fbmFtZSwgY29sdW1uKSk7XG4gICAgTW9kdWxlLl9mcmVlKGRhdGFIZWFwLmJ5dGVPZmZzZXQpO1xufVxuXG4vL25hbWU6IEJveCBDb3ggVHJhbnNmb3JtXG4vL2lucHV0OiBkYXRhZnJhbWUgZGF0YSBbSW5wdXQgZGF0YSB0YWJsZV1cbi8vaW5wdXQ6IGNvbHVtbiBjb2x1bW5fdG9fZmlsdGVyIFtTaWduYWwgdG8gRmlsdGVyXVxuLy9pbnB1dDogZG91YmxlIGxhbWJkYVxuLy9pbnB1dDogZG91YmxlIG9mc2V0XG5leHBvcnQgZnVuY3Rpb24gYm94X2NveF90cmFuc2Zvcm0oZGF0YSxjb2x1bW5fdG9fZmlsdGVyLGxhbWJkYSwgb2ZzZXQpIHtcbiAgICBsZXQganNfd3JhcHBlZF9zbWEgPSBNb2R1bGUuY3dyYXAoXCJib3hjb3hcIiwgXCJudWxsXCIsIFtcIm51bWJlclwiLFwibnVtYmVyXCIsXCJudW1iZXJcIixcIm51bWJlclwiXSk7XG4gICAgbGV0IGNvbHVtbl9uYW1lID0gY29sdW1uX3RvX2ZpbHRlci5uYW1lICsnIEJveCBDb3ggVHJhbnNmb3JtZWQnO1xuICAgIGxldCBmaWx0ZXJfYXJyYXkgPSBjb2x1bW5fdG9fZmlsdGVyLmdldFJhd0RhdGEoKTtcbiAgICBsZXQgbkRhdGFCeXRlcyA9IGZpbHRlcl9hcnJheS5sZW5ndGggKiBmaWx0ZXJfYXJyYXkuQllURVNfUEVSX0VMRU1FTlQ7XG4gICAgbGV0IGRhdGFQdHIgPSBNb2R1bGUuX21hbGxvYyhuRGF0YUJ5dGVzKTtcbiAgICBsZXQgZGF0YUhlYXAgPSBuZXcgVWludDhBcnJheShNb2R1bGUuSEVBUFU4LmJ1ZmZlciwgZGF0YVB0ciwgbkRhdGFCeXRlcyk7XG4gICAgZGF0YUhlYXAuc2V0KG5ldyBVaW50OEFycmF5KGZpbHRlcl9hcnJheS5idWZmZXIpKTtcbiAgICBqc193cmFwcGVkX3NtYShkYXRhUHRyLGZpbHRlcl9hcnJheS5sZW5ndGgsbGFtYmRhLG9mc2V0KTtcbiAgICBsZXQgcmVzdWx0ID0gbmV3IEZsb2F0MzJBcnJheShkYXRhSGVhcC5idWZmZXIsIGRhdGFIZWFwLmJ5dGVPZmZzZXQsIGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBjb2x1bW4gPSBuZXcgRmxvYXQzMkFycmF5KGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBpPTA7XG4gICAgZm9yIChpID0gMDsgaSA8IGZpbHRlcl9hcnJheS5sZW5ndGg7IGkrKykge2NvbHVtbltpXT1yZXN1bHRbaV07fVxuICAgIGRhdGEuY29sdW1ucy5hZGQoREcuQ29sdW1uLmZyb21GbG9hdDMyQXJyYXkoY29sdW1uX25hbWUsIGNvbHVtbikpO1xuICAgIE1vZHVsZS5fZnJlZShkYXRhSGVhcC5ieXRlT2Zmc2V0KTtcbn1cblxuLy9uYW1lOiBGb3VyaWVyIEZpbHRlclxuLy9pbnB1dDogZGF0YWZyYW1lIGRhdGEgW0lucHV0IGRhdGEgdGFibGVdXG4vL2lucHV0OiBjb2x1bW4gY29sdW1uX3RvX2ZpbHRlciBbU2lnbmFsIHRvIEZpbHRlcl1cbi8vaW5wdXQ6IGRvdWJsZSBsb3djdXRcbi8vaW5wdXQ6IGRvdWJsZSBoaWN1dFxuZXhwb3J0IGZ1bmN0aW9uIGZvdXJpZXJfZmlsdGVyKGRhdGEsY29sdW1uX3RvX2ZpbHRlcixsb3djdXQsIGhpY3V0KSB7XG4gICAgbGV0IGpzX3dyYXBwZWRfZmYgPSBNb2R1bGUuY3dyYXAoXCJmZmlsdGVyXCIsIFwibnVsbFwiLCBbXCJudW1iZXJcIixcIm51bWJlclwiLFwibnVtYmVyXCIsXCJudW1iZXJcIl0pO1xuICAgIGxldCBjb2x1bW5fbmFtZSA9IGNvbHVtbl90b19maWx0ZXIubmFtZSArJyBGb3VyaWVyIEZpbHRlcmVkIChMOiAnK2xvd2N1dCArICc7IEg6ICcgKyBoaWN1dCArICcpJztcbiAgICBsZXQgZmlsdGVyX2FycmF5ID0gY29sdW1uX3RvX2ZpbHRlci5nZXRSYXdEYXRhKCk7XG4gICAgbGV0IG5EYXRhQnl0ZXMgPSBmaWx0ZXJfYXJyYXkubGVuZ3RoICogZmlsdGVyX2FycmF5LkJZVEVTX1BFUl9FTEVNRU5UO1xuICAgIGxldCBkYXRhUHRyID0gTW9kdWxlLl9tYWxsb2MobkRhdGFCeXRlcyk7XG4gICAgbGV0IGRhdGFIZWFwID0gbmV3IFVpbnQ4QXJyYXkoTW9kdWxlLkhFQVBVOC5idWZmZXIsIGRhdGFQdHIsIG5EYXRhQnl0ZXMpO1xuICAgIGRhdGFIZWFwLnNldChuZXcgVWludDhBcnJheShmaWx0ZXJfYXJyYXkuYnVmZmVyKSk7XG4gICAganNfd3JhcHBlZF9mZihkYXRhUHRyLGZpbHRlcl9hcnJheS5sZW5ndGgsbG93Y3V0LGhpY3V0KTtcbiAgICBsZXQgcmVzdWx0ID0gbmV3IEZsb2F0MzJBcnJheShkYXRhSGVhcC5idWZmZXIsIGRhdGFIZWFwLmJ5dGVPZmZzZXQsIGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBjb2x1bW4gPSBuZXcgRmxvYXQzMkFycmF5KGZpbHRlcl9hcnJheS5sZW5ndGgpO1xuICAgIGxldCBpPTA7XG4gICAgZm9yIChpID0gMDsgaSA8IGZpbHRlcl9hcnJheS5sZW5ndGg7IGkrKykge2NvbHVtbltpXT1yZXN1bHRbaV07fVxuICAgIGRhdGEuY29sdW1ucy5hZGQoREcuQ29sdW1uLmZyb21GbG9hdDMyQXJyYXkoY29sdW1uX25hbWUsIGNvbHVtbikpO1xuICAgIE1vZHVsZS5fZnJlZShkYXRhSGVhcC5ieXRlT2Zmc2V0KTtcbn0iLCJtb2R1bGUuZXhwb3J0cyA9IERHOyIsIm1vZHVsZS5leHBvcnRzID0gZ3JvazsiLCJtb2R1bGUuZXhwb3J0cyA9IHVpOyJdLCJzb3VyY2VSb290IjoiIn0=