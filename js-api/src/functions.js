import { paramsToJs, toDart, toJs } from './wrappers';
import { DartWidget } from './widgets';
let api = window;
/** Grok functions */
export class Functions {
  register(func) {
    api.grok_RegisterFunc(func);
  }
  registerParamFunc(name, type, run, check = null, description = null) {
    api.grok_RegisterParamFunc(name, type, run, check, description);
  }
  call(name, parameters = {}, showProgress = false, progress = null) {
    return new Promise((resolve, reject) => api.grok_CallFunc(name, parameters, (out) => resolve(toJs(out)), (err) => reject(err), showProgress, toDart(progress)));
  }
  eval(name) {
    return new Promise((resolve, reject) => api.grok_EvalFunc(name, function (out) {
      return resolve(toJs(out));
    }, (err) => reject(err)));
  }
  scriptSync(s) {
    return toJs(api.grok_ScriptSync(s));
  }
}
/** Represents a function call
 * {@link https://datagrok.ai/help/overview/functions/function-call*}
 * */
export class FuncCall {
  constructor(d) {
    this.d = d;
  }
  get func() { return toJs(api.grok_FuncCall_Get_Func(this.d)); }
  set func(func) { api.grok_FuncCall_Get_Func(this.d, func.d); }
  /** Returns function call parameter value
     * @param {string} name
     * @returns {object} */
  getParamValue(name) {
    return toJs(api.grok_FuncCall_Get_Param_Value(this.d, name));
  }
  setParamValue(name, value) {
    api.grok_FuncCall_Set_Param_Value(this.d, name, toDart(value));
  }
  /** Executes the function call
     * @param {boolean} showProgress
     * @param {ProgressIndicator} progress
     * @returns {Promise<FuncCall>} */
  call(showProgress = false, progress = null) {
    return new Promise((resolve, reject) => api.grok_FuncCall_Call(this.d, (out) => resolve(toJs(out)), (err) => reject(err), showProgress, toDart(progress)));
  }
}
export function callFuncWithDartParameters(f, params) {
  let jsParams = paramsToJs(params);
  return f.apply(null, jsParams);
}
export class StepEditor extends DartWidget {
  constructor(d) {
    super(d);
  }
  static create() {
    return toJs(api.grok_StepEditor_Create());
  }
  loadScript(script) {
    return new Promise((resolve, reject) => api.grok_StepEditor_LoadScript(this.d, script, (out) => resolve(toJs(out)), (err) => reject(err)));
  }
  toScript() {
    return api.grok_StepEditor_ToScript(this.d);
  }
}
