define(['require', 'exports', './wrappers', './widgets'], function (require, exports, wrappers_1, widgets_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.StepEditor = exports.callFuncWithDartParameters = exports.FuncCall = exports.Functions = void 0;
  let api = window;
  /** Grok functions */
  class Functions {
    register(func) {
      api.grok_RegisterFunc(func);
    }
    registerParamFunc(name, type, run, check = null, description = null) {
      api.grok_RegisterParamFunc(name, type, run, check, description);
    }
    call(name, parameters = {}, showProgress = false, progress = null) {
      return new Promise((resolve, reject) => api.grok_CallFunc(name, parameters, (out) => resolve(wrappers_1.toJs(out)), (err) => reject(err), showProgress, wrappers_1.toDart(progress)));
    }
    eval(name) {
      return new Promise((resolve, reject) => api.grok_EvalFunc(name, function (out) {
        return resolve(wrappers_1.toJs(out));
      }, (err) => reject(err)));
    }
    scriptSync(s) {
      return wrappers_1.toJs(api.grok_ScriptSync(s));
    }
  }
  exports.Functions = Functions;
  /** Represents a function call
     * {@link https://datagrok.ai/help/overview/functions/function-call*}
     * */
  class FuncCall {
    constructor(d) {
      this.d = d;
    }
    get func() { return wrappers_1.toJs(api.grok_FuncCall_Get_Func(this.d)); }
    set func(func) { api.grok_FuncCall_Get_Func(this.d, func.d); }
    /** Returns function call parameter value
         * @param {string} name
         * @returns {object} */
    getParamValue(name) {
      return wrappers_1.toJs(api.grok_FuncCall_Get_Param_Value(this.d, name));
    }
    setParamValue(name, value) {
      api.grok_FuncCall_Set_Param_Value(this.d, name, wrappers_1.toDart(value));
    }
    /** Executes the function call
         * @param {boolean} showProgress
         * @param {ProgressIndicator} progress
         * @returns {Promise<FuncCall>} */
    call(showProgress = false, progress = null) {
      return new Promise((resolve, reject) => api.grok_FuncCall_Call(this.d, (out) => resolve(wrappers_1.toJs(out)), (err) => reject(err), showProgress, wrappers_1.toDart(progress)));
    }
    getEditor(condensed) {
      return new Promise((resolve, reject) => api.grok_FuncCall_Get_Editor(this.d, condensed, (out) => resolve(out), (err) => reject(err)));
    }
  }
  exports.FuncCall = FuncCall;
  function callFuncWithDartParameters(f, params) {
    let jsParams = wrappers_1.paramsToJs(params);
    return f.apply(null, jsParams);
  }
  exports.callFuncWithDartParameters = callFuncWithDartParameters;
  class StepEditor extends widgets_1.DartWidget {
    constructor(d) {
      super(d);
    }
    static create() {
      return wrappers_1.toJs(api.grok_StepEditor_Create());
    }
    loadScript(script) {
      return new Promise((resolve, reject) => api.grok_StepEditor_LoadScript(this.d, script, (out) => resolve(wrappers_1.toJs(out)), (err) => reject(err)));
    }
    toScript() {
      return api.grok_StepEditor_ToScript(this.d);
    }
  }
  exports.StepEditor = StepEditor;
});
