import {paramsToJs, toDart, toJs} from "./wrappers";
import {Type} from "./const";
import {Entity, Func} from "./entities";
import {DartWidget, ProgressIndicator,} from "./widgets";
import {MapProxy, _toIterable} from "./utils";
import {Observable} from "rxjs";
import {__obs} from "./events";
declare let grok: any;
declare let DG: any;
let api = <any>window;


const FuncCallParamMapProxy = new Proxy(class {
    d: any;
    input: boolean;
    constructor(d: any, input: boolean) {
      this.d = d;
      this.input = input;
    }
    keys(): Iterable<any> {
      return _toIterable(this.input ? api.grok_Func_InputParamMap_Keys(this.d) : api.grok_Func_OutputParamMap_Keys(this.d));
    }
    values() {
      return _toIterable(this.input ? api.grok_Func_InputParamMap_Values(this.d) : api.grok_Func_OutputParamMap_Values(this.d));
    }
    * [Symbol.iterator] () {
      for (let key of this.keys()) {
        const value = DG.toJs(this.input ? api.grok_Func_InputParamMap_Get(this.d, key) : api.grok_Func_OutputParamMap_Get(this.d, key));
        yield [key, value];
      }
    }
    entries() {
      return this;
    }
    forEach(callback: (key: string, value: any) => void) {
      for (const [key, value] of this) {
        callback(key, value);
      }
    }
    delete(key: string) {
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Delete(this.d, key) : api.grok_Func_OutputParamMap_Delete(this.d, key));
    }
    get(key: any) {
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Get(this.d, key) : api.grok_Func_OutputParamMap_Get(this.d, key));
    }
    has(key: string) {
      return this.input ? api.grok_Func_InputParamMap_Has(this.d, key) : api.grok_Func_OutputParamMap_Has(this.d, key);
    }
    set(key: string, value: any) {
      this.input ? api.grok_Func_InputParamMap_Set(this.d, key, DG.toDart(value)) : api.grok_Func_OutputParamMap_Set(this.d, key, DG.toDart(value));
      return this;
    }
    clear() {
      this.input ? api.grok_Func_InputParamMap_Clear(this.d) : api.grok_Func_OutputParamMap_Clear(this.d);
    }
    size() {
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Size(this.d) : api.grok_Func_OutputParamMap_Size(this.d));
    }
  }, {
    construct(target, args) {
      // @ts-ignore
      return new Proxy(new target(...args), {
        get: function (target: any, prop) {
          const val = target[prop];
          if (typeof val === 'function') {
            return function (...args :string[]) {
              return val.apply(target, args);
            };
          } else {
            return DG.toJs(target.input ? api.grok_Func_InputParamMap_Get(target.d, prop) : api.grok_Func_OutputParamMap_Get(target.d, prop));
          }
        },
        set: function (target, prop, value) {
          target.input ? api.grok_Func_InputParamMap_Set(target.d, prop, DG.toDart(value)) : api.grok_Func_OutputParamMap_Set(target.d, prop, DG.toDart(value));
          return true;
        },
        deleteProperty: function (target, prop) {
          target.input ? api.grok_Func_InputParamMap_Delete(target.d, DG.toDart(prop)) : api.grok_Func_OutputParamMap_Delete(target.d, DG.toDart(prop));
          return true;
        },
        has: function (target, prop) {
          return target.input ? api.grok_InputParamMap_Has(target.d, DG.toDart(prop)) : api.grok_OutputParamMap_Has(target.d, DG.toDart(prop));
        },
        getOwnPropertyDescriptor(target, prop) {
          return {
            enumerable: true,
            configurable: true
          };
        },
        ownKeys: function (target) {
          return Array.from(target.keys());
        }
      });
    }
  }
);


/** Grok functions */
export class Functions {
  register(func: Func): void {
    api.grok_RegisterFunc(func);
  }

  registerParamFunc(name: string, type: Type, run: Function, check: boolean | null = null, description: string | null = null): void {
    api.grok_RegisterParamFunc(name, type, run, check, description);
  }

  call(name: string, parameters: object = {}, showProgress: boolean = false, progress: ProgressIndicator | null = null): Promise<any> {
    return new Promise((resolve, reject) => api.grok_CallFunc(name, parameters, (out: any) => resolve(toJs(out)), (err: any) => reject(err), showProgress, toDart(progress)));
  }

  eval(name: string, context?: Context): Promise<any> {
    return new Promise((resolve, reject) => api.grok_EvalFunc(name, context?.d, function (out: any) {
      return resolve(toJs(out));
    }, (err: any) => reject(err)));
  }

  async find(name: string): Promise<any> {
    let f = await this.eval(name);
    if (f instanceof Func)
      return f;
    else
      return null;
  }

  scriptSync(s: string): any {
    return toJs(api.grok_ScriptSync(s));
  }

  getCurrentCall(): FuncCall {
    return toJs(api.grok_GetCurrentCall());
  }

  get onBeforeRunAction(): Observable<FuncCall> { return __obs('d4-before-run-action'); }
  get onAfterRunAction(): Observable<FuncCall> { return __obs('d4-after-run-action'); }
}

export class Context {
  readonly d: any;
  constructor(d: any) {
    this.d = d;
  }

  static create(): Context {
    return toJs(api.grok_Context_Create());
  }

  setVariable(name: string, value: any): void {
    api.grok_Context_Set_Variable(this.d, name, toDart(value));
  }

  getVariable(name: string): any {
    return toJs(api.grok_Context_Get_Variable(this.d, name));
  }

}

/** Represents a function call
 * {@link https://datagrok.ai/help/overview/functions/function-call*}
 * */
export class FuncCall {
  private readonly d: any;
  public inputs: any;
  public outputs: any;
  public aux: any;
  public options: any;

  constructor(d: any) {
    this.d = d;
    this.inputs = new FuncCallParamMapProxy(this.d, true);
    this.outputs = new FuncCallParamMapProxy(this.d, false);
    this.aux = new MapProxy(api.grok_FuncCall_Get_Aux(this.d));
    this.options = new MapProxy(api.grok_FuncCall_Get_Aux(this.d));
  }

  get func(): Func { return toJs(api.grok_FuncCall_Get_Func(this.d)); }
  set func(func: Func) {api.grok_FuncCall_Get_Func(this.d, func.d)}

  /** Returns function call parameter value
   * @param {string} name
   * @returns {object} */
  getParamValue(name: string): any {
    return toJs(api.grok_FuncCall_Get_Param_Value(this.d, name));
  }

  get context(): Context { return toJs(api.grok_FuncCall_Get_Context(this.d)); }
  set context(context: Context) { api.grok_FuncCall_Set_Context(this.d, context.d); }

  getOutputParamValue(): any {
    return toJs(api.grok_FuncCall_Get_Output_Param_Value(this.d));
  }

  setParamValue(name: string, value: any): void {
    api.grok_FuncCall_Set_Param_Value(this.d, name, toDart(value));
  }

  /** Executes the function call */
  call(showProgress: boolean = false, progress?: ProgressIndicator, options?: {processed: boolean}): Promise<FuncCall> {
    return new Promise((resolve, reject) => api.grok_FuncCall_Call(this.d, (out: any) => resolve(toJs(out)), (err: any) => reject(err), showProgress, toDart(progress), options?.processed));
  }

  edit() {
    api.grok_FuncCall_Edit(this.d);
  }

  getEditor(condensed?: boolean, showTableSelectors?: boolean): Promise<HTMLElement> {
    return new Promise((resolve, reject) => api.grok_FuncCall_Get_Editor(this.d, condensed, showTableSelectors, (out: any) => resolve(out), (err: any) => reject(err)));
  }
}

export function callFuncWithDartParameters<T>(f: (...params: any[]) => T, params: object): T {
  let jsParams = paramsToJs(params);
  return f.apply(null, jsParams);
}

export class StepEditor extends DartWidget {

  constructor(d: any) {
    super(d);
  }

  static create(): StepEditor {
    return toJs(api.grok_StepEditor_Create());
  }

  loadScript(script: String): Promise<void> {
    return new Promise((resolve, reject) => api.grok_StepEditor_LoadScript(this.d, script, (out: any) => resolve(toJs(out)), (err: any) => reject(err)));
  }

  toScript(): string {
    return api.grok_StepEditor_ToScript(this.d);
  }

}
