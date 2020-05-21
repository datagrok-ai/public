/** Grok functions */
import {paramsToJs} from "./utils";

export class Functions {
    register(func) { grok_RegisterFunc(func); }
    registerParamFunc(name, type, run, check = null, description = null) { grok_RegisterParamFunc(name, type, run, check, description); }
    call(name, parameters = {}, showProgress = false) { return new Promise((resolve, reject) => grok_CallFunc(name, parameters, (out) => resolve(out), showProgress)); }
    eval(name) { return new Promise((resolve, reject) => grok_EvalFunc(name, (out) => resolve(out))); }
}

export class FuncCall {
    constructor(d) { this.d = d; }

    getParamValue(name) { return grok_FuncCall_Get_Param_Value(this.d, name); }
}

export function callFuncWithDartParameters(f, params) {
    let jsParams = paramsToJs(params);
    return f.apply(null, jsParams);
}
