/**
 * Miscellaneous entity classes: Model, Notebook, Package, DockerContainer, ProgressIndicator.
 * @module entities/misc
 */

import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {Observable} from "rxjs";
import {observeStream} from "../events";
import {PackageLogger} from "../logger";
import type {FilesDataSource} from "../dapi";
import {Entity} from "./entity";
import {Group} from "./user";
import {Credentials} from "./data-connection";
import {DockerContainerStatus} from "./types";

declare var DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/-info}
 * */
export class Model extends Entity {
  /** @constructs Model */
  constructor(dart: any) {
    super(dart);
  }
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/compute/jupyter-notebook}
 * */
export class Notebook extends Entity {

  /** @constructs Notebook */
  constructor(dart: any) {
    super(dart);
  }

  /** Environment name */
  get environment(): string { return api.grok_Notebook_Get_Environment(this.dart); }
  set environment(e: string) { api.grok_Notebook_Set_Environment(this.dart, e); }

  /** Description */
  get description(): string { return api.grok_Notebook_Get_Description(this.dart); }
  set description(e: string) { api.grok_Notebook_Set_Description(this.dart, e); }

  get notebook(): {[_:string]: any} { return api.grok_Notebook_Get_Notebook_Content(this.dart); }
  set notebook(data: {[_:string]: any}) { api.grok_Notebook_Set_Notebook_Content(this.dart, data); }
}


/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package extends Entity {
  _webRoot: string | undefined;
  public _version: string = '';

  constructor(dart: any | undefined = undefined) {
    super(dart);

    if (typeof dart === 'string') {
      this._webRoot = dart;
      this.dart = null;
    }
  }

  /** Override init() method to provide package-specific initialization.
   * It is guaranteed to get called exactly once before the execution of any function below.
   */
  init(): Promise<null> { return Promise.resolve(null); }

  private _name: string = '';

  get webRoot(): string {
    if (this._webRoot === undefined)
      return api.grok_Package_Get_WebRoot(this.dart);
    else
      return this._webRoot;
  }

  set webRoot(x) {
    this._webRoot = x;
  }

  get packageOwner(): string {
    return api.grok_Package_Get_Package_Author(this.dart);
  }

  get version(): string {
    if (this.dart != null)
      return api.grok_Package_Get_Version(this.dart);
    else
      return this._version;
  }

  set version(x) {
    if (this.dart != null)
      api.grok_Package_Set_Version(this.dart, x);
    else
      this._version = x;
  }

  /** Package short name */
  get name(): string {
    if (this.dart != null)
      return api.grok_Entity_Get_Name(this.dart);
    else
      return this._name;
  }

  set name(x) {
    if (this.dart != null)
      api.grok_Entity_Set_Name(this.dart, x);
    else
      this._name = x;
  }

  getModuleName(file: string): string {
    if (this.dart != null)
      return api.grok_Package_GetModuleName(this.dart, file);
    else
      return '';
  }

  getIconUrl(): string {
    return api.grok_Package_GetIconUrl(this.dart);
  }

  /** Returns a JavaScript module for this package. */
  getModule(file: string): any {
    if (this.dart != null)
      return api.grok_Package_GetModule(this.dart, file)();
    else
      return null;
  }

  /** Returns metadata associated with the package.
   * The metadata gets generated when the package is built.
   * It is a concatenation of JSON files located under the /meta folder.
   * See example: /packages/PowerPack. */
  get meta(): {[key: string]: any} | null {
    return (this.dart == null) ? null : toJs(api.grok_Package_Get_Meta(this.dart));
  }

  /** Loads package. */
  async load(options?: {file: string}): Promise<Package> {
    return api.grok_Dapi_Packages_Load(this.dart, options?.file);
  }

  private _logger?: PackageLogger;

  get logger(): PackageLogger {
    if (this._logger)
      return this._logger;
    return this._logger = new PackageLogger(this);
  }

  /** Returns credentials for package. */
  getCredentials(): Promise<Credentials> {
    return api.grok_Package_Get_Credentials(this.name);
  }

  /**
   * @deprecated The {@link getProperties} should not be used. Use {@link settings} instead
   */
  getProperties(): Promise<any> {
    return this.getSettings();
  }

  /**
   * @deprecated The {@link getSettings} should not be used. Use {@link settings} instead
   */
  getSettings(): Promise<Map<string, any>> {
    return api.grok_Package_Get_Settings(this.name);
  }

  /** Returns settings for a package. */
  get settings(): {[index: string]: any} {
    return api.grok_Package_Get_Settings_Sync(this.name);
  }

  /** Updates settings for a package. */
  setSettings(props: Map<string, any>, group: Group): Promise<void> {
    return api.grok_Package_Set_Settings(this.name, props, group?.dart);
  }

  /** Global application data */
  get files(): FilesDataSource {
    return new DG.FilesDataSource(`System:AppData/${this.name}`);
  }

  public async getTests(core: boolean = false) {
    try {
      await this.load({ file: 'package-test.js' });
      let module = this.getModule('package-test.js');
      if (core && module.initAutoTests)
        await module.initAutoTests();
      return module.tests;
    } catch (e: any) {
      this.logger.error(e?.msg ?? 'get module error')
      return undefined;
    }
  }
}


// export class DockerImage extends Entity {
//   constructor(dart: any) {
//     super(dart);
//   }
// }


export class DockerContainer extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get status(): DockerContainerStatus {
    return api.grok_DockerContainer_Status(this.dart);
  }
}


export class ProgressIndicator {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create() {
    return toJs(api.grok_ProgressIndicator_Create());
  }

  get percent(): number {
    return api.grok_ProgressIndicator_Get_Percent(this.dart);
  }

  /** Flag indicating whether the operation was canceled by the user. */
  get canceled(): boolean { return api.grok_ProgressIndicator_Get_Canceled(this.dart); }

  get description(): string { return api.grok_ProgressIndicator_Get_Description(this.dart); }
  set description(s: string) { api.grok_ProgressIndicator_Set_Description(this.dart, s); }

  update(percent: number, description: string): void {
    api.grok_ProgressIndicator_Update(this.dart, percent, description);
  }

  log(line: string): void {
    api.grok_ProgressIndicator_Log(this.dart, line);
  }

  get onProgressUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Updated(this.dart));
  }

  get onLogUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Log_Updated(this.dart));
  }

  get onCanceled(): Observable<any> {
    return observeStream(api.grok_Progress_Canceled(this.dart));
  }
}
