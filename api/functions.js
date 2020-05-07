
/** Grok functions */
class Functions {
    register(func) { grok_RegisterFunc(func); }
    registerParamFunc(name, type, run, check = null, description = null) { grok_RegisterParamFunc(name, type, run, check, description); }
}

class FuncCall {
    constructor(d) { this.d = d; }

    getParamValue(name) { return grok_FuncCall_Get_Param_Value(this.d, name); }
}

function callFuncWithDartParameters(f, params) {
    let jsParams = paramsToJs(params);
    return f.apply(null, jsParams);
}

