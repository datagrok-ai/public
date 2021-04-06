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
define(['require', 'exports', './src/chem', './src/ml', './src/dapi', './src/functions', './src/events', './src/shell', './src/data'], function (require, exports, _chem, _ml, dapi_1, functions_1, events_1, shell_1, data_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.ml = exports.chem = exports.data = exports.settings = exports.shell = exports.dapi = exports.events = exports.functions = void 0;
  _chem = __importStar(_chem);
  _ml = __importStar(_ml);
  exports.functions = new functions_1.Functions();
  exports.events = new events_1.Events();
  exports.dapi = new dapi_1.Dapi();
  exports.shell = new shell_1.Shell();
  exports.settings = new shell_1.Settings();
  exports.data = new data_1.Data();
  exports.chem = _chem;
  exports.ml = _ml;
});
