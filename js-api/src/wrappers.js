define(['require', 'exports'], function (require, exports) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.toDart = exports.toJs = exports.paramsToJs = void 0;
  function paramsToJs(params) {
    return DG.paramsToJs(params);
  }
  exports.paramsToJs = paramsToJs;
  function toJs(d, check = false) {
    return DG.toJs(d);
  }
  exports.toJs = toJs;
  function toDart(x) {
    return DG.toDart(x);
  }
  exports.toDart = toDart;
});
