var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
  Object.defineProperty(o, 'default', { enumerable: true, value: v });
}) : function(o, v) {
  o['default'] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) if (k !== 'default' && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
  __setModuleDefault(result, mod);
  return result;
};
var __exportStar = (this && this.__exportStar) || function(m, exports) {
  for (var p in m) if (p !== 'default' && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
};
var __importDefault = (this && this.__importDefault) || function (mod) {
  return (mod && mod.__esModule) ? mod : { 'default': mod };
};
define(['require', 'exports', './src/chem', './src/ml', './src/utils', 'cash-dom', './src/const', './src/events', './src/dapi', './src/dataframe', './src/entities', './src/functions', './src/grid', './src/widgets', './src/view', './src/viewer', './src/docking', './src/wrappers_impl', './src/utils', './ui', './src/data'], function (require, exports, _chem, _ml, _utils, cash_dom_1, const_1, events_1, dapi_1, dataframe_1, entities_1, functions_1, grid_1, widgets_1, view_1, viewer_1, docking_1, wrappers_impl_1, utils_1, ui_1, data_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.utils = exports.ml = exports.chem = exports.Detector = exports.EntityMetaDartProxy = exports.ObjectHandler = exports.LruCache = exports.Utils = exports.timeAsync = exports.time = void 0;
  _chem = __importStar(_chem);
  _ml = __importStar(_ml);
  _utils = __importStar(_utils);
  cash_dom_1 = __importDefault(cash_dom_1);
  __exportStar(const_1, exports);
  __exportStar(events_1, exports);
  __exportStar(dapi_1, exports);
  __exportStar(dataframe_1, exports);
  __exportStar(entities_1, exports);
  __exportStar(functions_1, exports);
  __exportStar(grid_1, exports);
  __exportStar(widgets_1, exports);
  __exportStar(view_1, exports);
  __exportStar(viewer_1, exports);
  __exportStar(docking_1, exports);
  __exportStar(wrappers_impl_1, exports);
  Object.defineProperty(exports, 'time', { enumerable: true, get: function () { return utils_1.time; } });
  Object.defineProperty(exports, 'timeAsync', { enumerable: true, get: function () { return utils_1.timeAsync; } });
  Object.defineProperty(exports, 'Utils', { enumerable: true, get: function () { return utils_1.Utils; } });
  Object.defineProperty(exports, 'LruCache', { enumerable: true, get: function () { return utils_1.LruCache; } });
  Object.defineProperty(exports, 'ObjectHandler', { enumerable: true, get: function () { return ui_1.ObjectHandler; } });
  Object.defineProperty(exports, 'EntityMetaDartProxy', { enumerable: true, get: function () { return ui_1.EntityMetaDartProxy; } });
  Object.defineProperty(exports, 'Detector', { enumerable: true, get: function () { return data_1.Detector; } });
  exports.chem = _chem;
  exports.ml = _ml;
  exports.utils = _utils;
  cash_dom_1.default(function () {
    window.$ = cash_dom_1.default;
  });
});
