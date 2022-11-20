import {paramsToJs, toDart, toJs} from "./wrappers";
import {Type} from "./const";
import {Entity, Func, Property, User} from "./entities";
import {DartWidget, InputBase, ProgressIndicator,} from "./widgets";
import {MapProxy, _toIterable} from "./utils";
import {Observable} from "rxjs";
import {__obs, StreamSubscription} from "./events";
import * as rxjs from "rxjs";
import dayjs from "dayjs";
declare let grok: any;
declare let DG: any;
let api = <any>window;


const FuncCallParamMapProxy = new Proxy(class {
    dart: any;
    input: boolean;
    constructor(dart: any, input: boolean) {
      this.dart = dart;
      this.input = input;
    }
    keys(): Iterable<any> {
      return _toIterable(this.input ? api.grok_Func_InputParamMap_Keys(this.dart) : api.grok_Func_OutputParamMap_Keys(this.dart));
    }
    values() {
      return _toIterable(this.input ? api.grok_Func_InputParamMap_Values(this.dart) : api.grok_Func_OutputParamMap_Values(this.dart));
    }
    * [Symbol.iterator] () {
      for (let key of this.keys()) {
        const value = DG.toJs(this.input ? api.grok_Func_InputParamMap_Get(this.dart, key) : api.grok_Func_OutputParamMap_Get(this.dart, key));
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
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Delete(this.dart, key) : api.grok_Func_OutputParamMap_Delete(this.dart, key));
    }
    get(key: any) {
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Get(this.dart, key) : api.grok_Func_OutputParamMap_Get(this.dart, key));
    }
    has(key: string) {
      return this.input ? api.grok_Func_InputParamMap_Has(this.dart, key) : api.grok_Func_OutputParamMap_Has(this.dart, key);
    }
    set(key: string, value: any) {
      this.input ? api.grok_Func_InputParamMap_Set(this.dart, key, DG.toDart(value)) : api.grok_Func_OutputParamMap_Set(this.dart, key, DG.toDart(value));
      return this;
    }
    clear() {
      this.input ? api.grok_Func_InputParamMap_Clear(this.dart) : api.grok_Func_OutputParamMap_Clear(this.dart);
    }
    size() {
      return DG.toJs(this.input ? api.grok_Func_InputParamMap_Size(this.dart) : api.grok_Func_OutputParamMap_Size(this.dart));
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
            return DG.toJs(target.input ? api.grok_Func_InputParamMap_Get(target.dart, prop) : api.grok_Func_OutputParamMap_Get(target.dart, prop));
          }
        },
        set: function (target, prop, value) {
          target.input ? api.grok_Func_InputParamMap_Set(target.dart, prop, DG.toDart(value)) : api.grok_Func_OutputParamMap_Set(target.dart, prop, DG.toDart(value));
          return true;
        },
        deleteProperty: function (target, prop) {
          target.input ? api.grok_Func_InputParamMap_Delete(target.dart, DG.toDart(prop)) : api.grok_Func_OutputParamMap_Delete(target.dart, DG.toDart(prop));
          return true;
        },
        has: function (target, prop) {
          return target.input ? api.grok_InputParamMap_Has(target.dart, DG.toDart(prop)) : api.grok_OutputParamMap_Has(target.dart, DG.toDart(prop));
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

  async call(name: string, parameters: object = {}, showProgress: boolean = false, progress: ProgressIndicator | null = null): Promise<any> {
    return toJs(await api.grok_CallFunc(name, parameters, showProgress, toDart(progress)));
  }

  async eval(name: string, context?: Context): Promise<any> {
    return toJs(await api.grok_EvalFunc(name, context?.dart));
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


/** Represents a parameter of a function. */
export class FuncCallParam {
  readonly dart: any;

  /** Auxiliary data used for storing additional information associated with this parameter. */
  public aux: any;

  constructor(dart: any) {
    this.dart = dart;
    this.aux = new MapProxy(api.grok_FuncCallParam_Get_Aux(this.dart));
  }

  get name(): string { return this.property.name; }

  get value(): any { return toJs(api.grok_FuncCallParam_Get_Value(this.dart)); }

  /** A property that re*/
  get property(): Property { return toJs(api.grok_FuncCallParam_Get_Param(this.dart)); }

  processOutput(): void {
    api.grok_FuncCallParam_ProcessOutput(this.dart);
  }

  get onChanged(): Observable<any> {
    let object = this.dart;
    return rxjs.fromEventPattern(
      function (handler) {
        return api.grok_FuncCallParam_OnChanged(object, function (x: any) {
          handler(toJs(x));
        });
      },
      function (handler, dart) {
        new StreamSubscription(dart).cancel();
      }
    );
  }
}


export class Context {
  readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(): Context {
    return toJs(api.grok_Context_Create());
  }

  static cloneDefault(): Context {
    return toJs(api.grok_Context_CloneDefault());
  }

  setVariable(name: string, value: any): void {
    api.grok_Context_Set_Variable(this.dart, name, toDart(value));
  }

  getVariable(name: string): any {
    return toJs(api.grok_Context_Get_Variable(this.dart, name));
  }
}

/** Represents a function call
 * {@link https://datagrok.ai/help/datagrok/functions/function-call*}
 * */
export class FuncCall extends Entity {
  public readonly dart: any;
  public inputs: any;
  public outputs: any;
  public aux: any;
  public options: any;
  public inputParams: any;
  public outputParams: any;

  constructor(dart: any) {
    super(dart);
    this.inputs = new FuncCallParamMapProxy(this.dart, true);
    this.outputs = new FuncCallParamMapProxy(this.dart, false);
    this.inputParams = new MapProxy(api.grok_FuncCall_Get_Params(this.dart, true));
    this.outputParams = new MapProxy(api.grok_FuncCall_Get_Params(this.dart, false));
    this.aux = new MapProxy(api.grok_FuncCall_Get_Aux(this.dart));
    this.options = new MapProxy(api.grok_FuncCall_Get_Options(this.dart));
  }

  /** Function this call is associated with. */
  get func(): Func { return toJs(api.grok_FuncCall_Get_Func(this.dart)); }
  set func(func: Func) { api.grok_FuncCall_Set_Func(this.dart, func.dart) }

  get parentCall(): FuncCall { return toJs(api.grok_FuncCall_Get_ParentCall(this.dart)); }
  set parentCall(c: FuncCall) {api.grok_FuncCall_Set_ParentCall(this.dart, c.dart)}

  get started(): dayjs.Dayjs { return dayjs(api.grok_FuncCall_Get_Started(this.dart)); }
  get finished(): dayjs.Dayjs { return dayjs(api.grok_FuncCall_Get_Finished(this.dart)); }

  override get author(): User { return toJs(api.grok_FuncCall_Get_Author(this.dart)) }

  /** Returns function call parameter value
   * @param {string} name
   * @returns {object} */
  getParamValue(name: string): any {
    return toJs(api.grok_FuncCall_Get_Param_Value(this.dart, name));
  }

  /** Function call execution context. */
  get context(): Context { return toJs(api.grok_FuncCall_Get_Context(this.dart)); }
  set context(context: Context) { api.grok_FuncCall_Set_Context(this.dart, context.dart); }

  /** Error message, if this call resulted in an exception, or null. */
  get errorMessage(): string | null { return api.grok_FuncCall_Get_ErrorMessage(this.dart); }

  getOutputParamValue(): any {
    return toJs(api.grok_FuncCall_Get_Output_Param_Value(this.dart));
  }

  setParamValue(name: string, value: any): void {
    api.grok_FuncCall_Set_Param_Value(this.dart, name, toDart(value));
  }

  /** Executes the function call */
  call(showProgress: boolean = false, progress?: ProgressIndicator, options?: {processed?: boolean, report?: boolean}): Promise<FuncCall> {
    return new Promise((resolve, reject) => api.grok_FuncCall_Call(this.dart, (out: any) => resolve(toJs(out)), (err: any) => reject(err), showProgress, toDart(progress), options?.processed, options?.report));
  }

  /** Executes the function call synchronously*/
  callSync(options?: {processed?: boolean, report?: boolean}):FuncCall {
     return api.grok_FuncCall_Call_Sync(this.dart, options?.processed, options?.report);
  }

  /** Shows the corresponding dialog (or view). */
  edit() { api.grok_FuncCall_Edit(this.dart); }

  getEditor(condensed?: boolean, showTableSelectors?: boolean): Promise<HTMLDivElement> {
    return api.grok_FuncCall_Get_Editor(this.dart, condensed, showTableSelectors);
  }

  buildEditor(root: HTMLDivElement, options?: {condensed?: boolean, showTableSelectors?: boolean}): Promise<InputBase[]> {
    return api.grok_FuncCall_Build_Editor(this.dart, root, options?.condensed, options?.showTableSelectors);
  }

  /** Makes a shallow copy. */
  clone(): FuncCall { return api.grok_FuncCall_Clone(this.dart); }
}

export function callFuncWithDartParameters<T>(f: (...params: any[]) => T, params: object, dartResult: boolean): T {
  let jsParams = paramsToJs(params);
  if (dartResult)
    return toDart(f.apply(null, jsParams));
  return f.apply(null, jsParams);
}

export class StepEditor extends DartWidget {

  constructor(dart: any) {
    super(dart);
  }

  static create(): StepEditor {
    return toJs(api.grok_StepEditor_Create());
  }

  loadScript(script: String): Promise<void> {
    return new Promise((resolve, reject) => api.grok_StepEditor_LoadScript(this.dart, script, (out: any) => resolve(toJs(out)), (err: any) => reject(err)));
  }

  toScript(): string {
    return api.grok_StepEditor_ToScript(this.dart);
  }

}
