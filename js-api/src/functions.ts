import {paramsToJs, toDart, toJs} from "./wrappers";
import {Type} from "./const";
import {Entity, Func, Property, User} from "./entities";
import {DartWidget, InputBase, ProgressIndicator,} from "./widgets";
import {MapProxy, _toIterable} from "./utils";
import {Observable} from "rxjs";
import {__obs, StreamSubscription} from "./events";
import * as rxjs from "rxjs";
import dayjs from "dayjs";
import {IDartApi} from "./api/grok_api.g";
import {ViewBase} from "./views/view";
declare let grok: any;
declare let DG: any;
const api: IDartApi = <any>window;



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
        get: function (target: any, prop: string) {
          const val = target[prop];
          if (typeof val === 'function') {
            return function (...args :string[]) {
              return val.apply(target, args);
            };
          } else {
            return DG.toJs(target.input ? api.grok_Func_InputParamMap_Get(target.dart, prop) : api.grok_Func_OutputParamMap_Get(target.dart, prop));
          }
        },
        set: function (target, prop: string, value) {
          target.input ? api.grok_Func_InputParamMap_Set(target.dart, prop, DG.toDart(value)) : api.grok_Func_OutputParamMap_Set(target.dart, prop, DG.toDart(value));
          return true;
        },
        deleteProperty: function (target, prop: string) {
          target.input ? api.grok_Func_InputParamMap_Delete(target.dart, DG.toDart(prop)) : api.grok_Func_OutputParamMap_Delete(target.dart, DG.toDart(prop));
          return true;
        },
        has: function (target, prop: string) {
          return target.input ? api.grok_Func_InputParamMap_Has(target.dart, DG.toDart(prop)) : api.grok_Func_OutputParamMap_Has(target.dart, DG.toDart(prop));
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


export interface IFunctionRegistrationData {
  signature: string;    // int foo(string bar)
  run: Function;
  tags?: string;        // comma-separated tags
  isAsync?: boolean;    // whether is can be called synchronously
  namespace?: string;
  options?: {[key: string]: string}
}


export interface IFunctionCallOptions {
  /** Function call context. */
  context: Context;

  /** Specifies if this call could be cached, even if the parent function/connection is not. */
  cacheable?: boolean;

  /** Whether Func.beforeCommandExecuted and afterCommandExecuted events are fired. */
  reportable?: boolean;

  progress?: ProgressIndicator;
}


/** Grok functions */
export class Functions {

  /** Controls client caching. */
  get clientCache(): ClientCache { return new ClientCache(); }

  /** Registers a function globally. */
  register(func: IFunctionRegistrationData): Func {
    return new Func(api.grok_RegisterFunc(func));
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

  parse(command: string, safe: boolean = true): any {
    return toJs(api.grok_Parse_Command(command, safe));
  }

  handleOuterBracketsInColName(name: string, escape: boolean): string {
    return api.grok_ColumnName_HandleOuterBrackets(name, escape);
  }

 /** Returns a function with the specified name, or throws an error if
   * there is no such function. See also {@link find}. */
  async get(name: string): Promise<Func> {
    let f = await this.find(name);
    if (!f)
      throw `Function not found: "${name}"`;
    return f;
  }

  /** Returns a function with the specified name, or null if not found.
   * See also {@link find}. */
  async find(name: string): Promise<Func | null> {
    let f = await this.eval(name);
    return (f instanceof Func) ? f : null;
  }

  scriptSync(s: string): any {
    return toJs(api.grok_ScriptSync(s));
  }

  getCurrentCall(): FuncCall {
    return toJs(api.grok_GetCurrentCall());
  }

  get onBeforeRunAction(): Observable<FuncCall> { return __obs('d4-before-run-action'); }
  get onAfterRunAction(): Observable<FuncCall> { return __obs('d4-after-run-action'); }
  get onParamsUpdated(): Observable<FuncCall> { return __obs('d4-func-call-output-params-updated'); }
}


/** Client caching service that caches results of function invocations and stores
 * them in the IndexedDb. */
export class ClientCache {

  /** Clears cache content. */
  clear(metaId?:string): Promise<void> { return api.grok_ClientCache_Clear(metaId); }

  /** Starts client function caching service. */
  start(): Promise<void> { return api.grok_ClientCache_Start(); }

  /** Stops client function caching service. */
  stop(): void { api.grok_ClientCache_Stop(); }

  /** Removes expired records. Normally, Datagrok does it automatically when needed. */
  cleanup(): Promise<void> { return api.grok_ClientCache_Cleanup(); }

  /** Returns the number of */
  getRecordCount(): Promise<number> { return api.grok_ClientCache_GetRecordCount(); }

  /** Indicates whether the caching service is running. */
  get isRunning() { return api.grok_ClientCache_Get_IsRunning(); }
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

type FuncCallParams = {
  [name: string]: FuncCallParam;
} & { values(): FuncCallParam[]; }

/** Represents a function call
 * {@link https://datagrok.ai/help/datagrok/functions/function-call*}
 * */
export class FuncCall extends Entity {
  declare readonly dart: any;

  /** Named input values. See {@link inputParams} for parameter metadata. */
  public inputs: {[name: string]: any};

  /** Named output values. See {@link outputParams} for parameter metadata. */
  public outputs: {[name: string]: any};

  /** Input parameter metadata. See {@link inputs} for parameter values. */
  public inputParams: FuncCallParams;

  /** Output parameter metadata. See {@link outputs} for parameter values. */
  public outputParams: FuncCallParams;

  public aux: any;
  public options: any;

  constructor(dart: any) {
    super(dart);
    // @ts-ignore
    this.inputs = new FuncCallParamMapProxy(this.dart, true);
    // @ts-ignore
    this.outputs = new FuncCallParamMapProxy(this.dart, false);
    // @ts-ignore
    this.inputParams = new MapProxy(api.grok_FuncCall_Get_Params(this.dart, true));
    // @ts-ignore
    this.outputParams = new MapProxy(api.grok_FuncCall_Get_Params(this.dart, false));

    this.aux = new MapProxy(api.grok_FuncCall_Get_Aux(this.dart));
    this.options = new MapProxy(api.grok_FuncCall_Get_Options(this.dart));
  }

  /** Function this call is associated with. */
  get func(): Func { return toJs(api.grok_FuncCall_Get_Func(this.dart)); }
  set func(func: Func) { api.grok_FuncCall_Set_Func(this.dart, func.dart); }

  get parentCall(): FuncCall { return toJs(api.grok_FuncCall_Get_ParentCall(this.dart)); }
  set parentCall(c: FuncCall) { api.grok_FuncCall_Set_ParentCall(this.dart, c.dart); }

  get started(): dayjs.Dayjs { return dayjs(api.grok_FuncCall_Get_Started(this.dart)); }

  set started(value: dayjs.Dayjs) {
    if (!(dayjs.isDayjs(value) || value == null))
      value = dayjs(value);
    api.grok_FuncCall_Set_Started(this.dart, value?.valueOf());
  }

  get finished(): dayjs.Dayjs { return dayjs(api.grok_FuncCall_Get_Finished(this.dart)); }

  get status(): string { return api.grok_FuncCall_Get_Status(this.dart); }

  get adHoc(): boolean { return api.grok_FuncCall_Get_AdHoc(this.dart); }
  set adHoc(a: boolean) { api.grok_FuncCall_Set_AdHoc(this.dart, a); }

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

  /** Returns the first output parameter value, or null. */
  getOutputParamValue(): any {
    return toJs(api.grok_FuncCall_Get_Output_Param_Value(this.dart));
  }

  setAuxValue(name: string, value: any): void {
    api.grok_FuncCall_Set_Aux_Value(this.dart, name, toDart(value));
  }

  setParamValue(name: string, value: any): void {
    api.grok_FuncCall_Set_Param_Value(this.dart, name, toDart(value));
  }

  /** Executes the function call */
  async call(showProgress: boolean = false, progress?: ProgressIndicator, options?: {processed?: boolean, report?: boolean}): Promise<FuncCall> {
    await api.grok_FuncCall_Call(this.dart, showProgress, toDart(progress), options?.processed, options?.report);
    return this;
  }

  cancel(): Promise<void> {
    return api.grok_FuncCall_Cancel(this.dart);
  }

  /** Executes the function call synchronously*/
  callSync(options?: {processed?: boolean, report?: boolean}): FuncCall {
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

  /** Returns views with result. Should be called on succeeded FuncCall **/
  getResultViews(): ViewBase[] {
    return toJs(api.grok_FuncCall_GetOutputViews(this.dart));
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
