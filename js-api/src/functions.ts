import {paramsToJs, toDart, toJs} from "./wrappers";
import {Type} from "./const";
import {Func} from "./entities";
declare let grok: any;
declare let DG: any;
let api = <any>window;

/** Grok functions */
export class Functions {
  register(func: Func) {
    api.grok_RegisterFunc(func);
  }

  registerParamFunc(name: string, type: Type, run: Function, check: boolean | null = null, description: string | null = null) {
    api.grok_RegisterParamFunc(name, type, run, check, description);
  }

  call(name: string, parameters = {}, showProgress = false, progress = null) {
    return new Promise((resolve, reject) => api.grok_CallFunc(name, parameters, (out: any) => resolve(toJs(out)), (err: any) => reject(err), showProgress, toDart(progress)));
  }

  eval(name: string) {
    return new Promise((resolve, reject) => api.grok_EvalFunc(name, function (out: any) {
      return resolve(toJs(out));
    }, (err: any) => reject(err)));
  }

  scriptSync(s: string) {
    return toJs(api.grok_ScriptSync(s));
  }
}

/** Represents a function call
 * {@link https://datagrok.ai/help/overview/functions/function-call*}
 * */
export class FuncCall {
  private readonly d: any;
  constructor(d: any) {
    this.d = d;
  }

  /** Returns function call parameter value
   * @param {string} name
   * @returns {object} */
  getParamValue(name: string) {
    return toJs(api.grok_FuncCall_Get_Param_Value(this.d, name));
  }

  /** Executes the function call
   * @param {boolean} showProgress
   * @param {ProgressIndicator} progress
   * @returns {Promise<FuncCall>} */
  call(showProgress = false, progress = null) {
    return new Promise((resolve, reject) => api.grok_FuncCall_Call(this.d, (out: any) => resolve(toJs(out)), (err: any) => reject(err), showProgress, toDart(progress)));
  }
}

export function callFuncWithDartParameters(f: Function, params: object) {
  let jsParams = paramsToJs(params);
  return f.apply(null, jsParams);
}
