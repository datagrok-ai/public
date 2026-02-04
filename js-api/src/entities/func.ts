/**
 * Func, Script, and ScriptEnvironment classes.
 * @module entities/func
 */

import {ScriptingLanguage} from "../const";
import {toDart, toJs} from "../wrappers";
import {FuncCall} from "../functions";
import {MapProxy} from "../proxies";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import {Property} from "./property";

declare var grok: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

// Forward declaration
type Package = any;


/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/datagrok/functions/function}
 * */
export class Func extends Entity {
  public aux: any;
  public options: { [key: string]: any; };

  constructor(dart: any) {
    super(dart);
    this.aux = new MapProxy(api.grok_Func_Get_Aux(this.dart));
    // @ts-ignore
    this.options = new MapProxy(api.grok_Func_Get_Options(this.dart));
  }

  get description(): string { return api.grok_Func_Get_Description(this.dart); }

  get type(): string { return api.grok_Func_Get_Type(this.dart); }

  get path(): string { return api.grok_Func_Get_Path(this.dart); }

  /** Help URL. */
  get helpUrl(): string { return api.grok_Func_Get_HelpUrl(this.dart); }

  set helpUrl(url: string) { api.grok_Func_Set_HelpUrl(this.dart, url); }

  /** A package this function belongs to. */
  get package(): Package { return api.grok_Func_Get_Package(this.dart); }

  /** Indicates that the function (or script) is already vector, meaning it
   * accepts vector input (an entire column) and processes it in a single call,
   * rather than being executed separately for each scalar element (row) */
  get isVectorFunc(): boolean { return api.grok_Func_Get_IsVectorFunc(this.dart); }

  /** Returns {@link FuncCall} object in a stand-by state */
  prepare(parameters: {[name: string]: any} = {}): FuncCall {
    return toJs(api.grok_Func_Prepare(this.dart, parameters));
  };

  async prepareAsync(parameters: {[name: string]: any} = {}): Promise<FuncCall> {
    const call = await api.grok_Func_PrepareAsync(this.dart, parameters);
    return toJs(call);
  };

  /** Input parameters */
  get inputs(): Property[] {
    return toJs(api.grok_Func_Get_InputParams(this.dart));
  }

  /** Output parameters */
  get outputs(): Property[] {
    return toJs(api.grok_Func_Get_OutputParams(this.dart));
  }

  /**
   *  Executes the function with the specified {@link parameters}, and returns result.
   *  If necessary, the corresponding package will be loaded as part of the call.
   * */
  async apply(parameters: {[name: string]: any} | any[] = {}): Promise<any> {
    parameters ??= {};
    if (Array.isArray(parameters)) {
      if (parameters.length != this.inputs.length)
        throw `${this.name}: expected ${this.inputs.length} parameters, got ${parameters.length}`;
      const params: any = {};
      for (let i = 0; i < this.inputs.length; i++)
        params[this.inputs[i].name] = parameters[i];
      parameters = params;
    }

    return (await (this.prepare(parameters)).call()).getOutputParamValue();
  }

  /** Executes the function synchronously, and returns the result.
   *  If the function is asynchronous, throws an exception. */
  applySync(parameters: {[name: string]: any} = {}): any {
    return this.prepare(parameters).callSync().getOutputParamValue();
  }

  /** Returns functions with the specified attributes. */
  static find(params?: { package?: string, name?: string, tags?: string[], meta?: any, returnType?: string, returnSemType?: string}): Func[] {
    return api.grok_Func_Find(params?.package, params?.name, params?.tags, params?.meta, params?.returnType, params?.returnSemType);
  }

  /**
   * @deprecated Use find, it's the same now but does not make a server query and synchronous.
   */
  static async findAll(params?: { package?: string, name?: string, tags?: string[], meta?: any, returnType?: string, returnSemType?: string}): Promise<Func[]> {
    let functions = Func.find(params);
    let queries = await grok.dapi.queries.include('params,connection').filter(`name="${params?.name}"`).list();
    let scripts = await grok.dapi.scripts.include('params').filter(`name="${params?.name}"`).list();

    return [...functions, ...queries, ...scripts];
  }

  /** Returns a function with the specified name. */
  static byName(name: string): Func {
    return Func.find({ name: name})[0];
  }
}


/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
  public static readonly vecInputTableName = 'in_vec_table';
  public static readonly vecOutputTableName = 'out_vec_table';

  /** @constructs Script */
  constructor(dart: any) {
    super(dart);
  }

  static create(script: string): Script { return new Script(api.grok_Script_Create(script)); }

  static fromParams(inputs: Property[], outputs: Property[], script: string = ''): Script {
    return new Script(api.grok_Script_FromParams(inputs.map((i) => toDart(i)), outputs.map((i) => toDart(i)), script));
  }

  /** Script */
  get script(): string { return api.grok_Script_GetScript(this.dart); }
  set script(s: string) { api.grok_Script_SetScript(this.dart, s); }

  /** Script */
  get clientCode(): string { return api.grok_Script_ClientCode(this.dart); }

  /** Script language. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get language(): ScriptingLanguage { return api.grok_Script_GetLanguage(this.dart); }
  set language(s: ScriptingLanguage) { api.grok_Script_SetLanguage(this.dart, s); }

  /** Indicates that the script is already vector, meaning it accepts vector
   * input (an entire column) and processes it in a single call, rather than
   * being executed separately for each scalar element (row) */
  get isVectorFunc(): boolean { return api.grok_Script_Get_IsVectorFunc(this.dart); }

  /** Environment name. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get environment(): string { return api.grok_Script_Get_Environment(this.dart); }
  set environment(s: string) { api.grok_Script_Set_Environment(this.dart, s); }

  /** Reference header parameter. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get reference(): string { return api.grok_Script_Get_Reference(this.dart); }
  set reference(s: string) { api.grok_Script_Set_Reference(this.dart, s); }

  /** Sample table. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get sample(): string { return api.grok_Script_Get_Sample(this.dart); }
  set sample(s: string) { api.grok_Script_Set_Sample(this.dart, s); }

  /** Script tags. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get tags(): string[] { return api.grok_Script_Get_Tags(this.dart); }
  set tags(tags: string[]) { api.grok_Script_Set_Tags(this.dart, tags); }
}


/** Represents a script environment */
export class ScriptEnvironment extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Create instance of ScriptEnvironment */
  static create(name: string): ScriptEnvironment {
    return new ScriptEnvironment(api.grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content */
  get environment(): string { return api.grok_ScriptEnvironment_Environment(this.dart); }
}
