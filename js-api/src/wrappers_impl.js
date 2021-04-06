define(['require', 'exports', './entities', './const', './viewer'], function (require, exports, entities_1, const_1, viewer_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.toDart = exports.toJs = exports.paramsToJs = void 0;
  /** Converts list of Dart objects to JavaScript objects by calling {@link toJs}
     * @param {object[]} params
     * @returns {object[]} - list of JavaScript objects */
  function paramsToJs(params) {
    let result = [];
    for (let i = 0; i < params.length; i++) {
      let type = window.grok_GetType(params[i]);
      if (type !== null && (!const_1.TYPES_SCALAR.has(type) || type === const_1.TYPE.LIST || type === const_1.TYPE.MAP))
        result.push(toJs(params[i]));
      else
        result.push(params[i]);
    }
    return result;
  }
  exports.paramsToJs = paramsToJs;
  /**
     * Instantiates the corresponding JS handler for the Dart object [d]. See also {@link toDart}
     * @param d - Dart handle
     * @param {boolean} check - when true, throws an exception if the object can't be converted to JS.
     * @returns JavaScript wrapper for the Dart object
     * */
  function toJs(d, check = false) {
    let type = window.grok_GetType(d);
    if (type === const_1.TYPE.MAP) {
      let wrapper = window.grok_GetWrapper(d);
      for (let key in wrapper) {
        if (wrapper.hasOwnProperty(key)) {
          let type = window.grok_GetType(wrapper[key]);
          if (type !== null && (!const_1.TYPES_SCALAR.has(type) || type === const_1.TYPE.LIST || type === const_1.TYPE.MAP))
            wrapper[key] = toJs(wrapper[key]);
        }
      }
      return wrapper;
    }
    else if (type === const_1.TYPE.LIST) {
      return paramsToJs(d);
    }
    else if (type === 'DG.TypedEventArgs')
      return new viewer_1.TypedEventArgs(d);
    let wrapper = window.grok_GetWrapper(d);
    if (wrapper != null)
      return wrapper;
    if (type === const_1.TYPE.PROPERTY)
      return new entities_1.Property(d);
    if (check)
      throw `Not supported type: ${type}`;
    return d;
  }
  exports.toJs = toJs;
  /** Extracts a Dart handle from the JavaScript wrapper. See also {@link toJs} */
  function toDart(x) {
    if (x === undefined || x === null)
      return x;
    return (typeof x.d !== 'undefined') ? x.d : x;
  }
  exports.toDart = toDart;
});
