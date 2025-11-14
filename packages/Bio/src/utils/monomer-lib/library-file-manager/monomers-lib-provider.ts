/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {JSONSchemaType} from 'ajv';
import {MonomerLibFileValidator} from './file-validator';
import {DEFAULT_FILES_LIB_PROVIDER_NAME, IMonomerLibHelper, IMonomerLibProvider, IMonomerLinkData,
  IMonomerSet, Monomer, MonomerLibData, MonomerSet, MonomerSetPlaceholder}
  from '@datagrok-libraries/bio/src/types/monomer-library';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  HELM_REQUIRED_FIELD as REQ,
} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {Observable, Subject} from 'rxjs';

export class MonomerLibFromFilesProvider implements IMonomerLibProvider {
  protected get LIB_PATH(): string { return 'System:AppData/Bio/monomer-libraries/'; }
  protected get SETS_PATH(): string { return 'System:AppData/Bio/monomer-sets/'; }
  public static readonly HELM_JSON_SCHEMA_PATH = 'System:AppData/Bio/tests/libraries/HELMmonomerSchema.json';
  private static _instance: MonomerLibFromFilesProvider | null = null;
  public filesPromise: Promise<void> = Promise.resolve();
  name: string = DEFAULT_FILES_LIB_PROVIDER_NAME;
  private _onChanged = new Subject<void>();
  get onChanged(): Observable<void> {
    return this._onChanged;
  }
  private constructor(
        private readonly fileValidator: MonomerLibFileValidator,
        private readonly libHelper: IMonomerLibHelper,
        private readonly logger: ILogger,
  ) {

  }

  public static async getInstance(
    libHelper: IMonomerLibHelper, logger: ILogger,
  ): Promise<MonomerLibFromFilesProvider> {
    if (this._instance == null) {
      const helmSchemaRaw = await grok.dapi.files.readAsText(MonomerLibFromFilesProvider.HELM_JSON_SCHEMA_PATH);
      const helmSchema = JSON.parse(helmSchemaRaw) as JSONSchemaType<any>;

      const fileValidator = new MonomerLibFileValidator(helmSchema);
      MonomerLibFromFilesProvider._instance = new MonomerLibFromFilesProvider(fileValidator, libHelper, logger);
      MonomerLibFromFilesProvider._instance.refreshLists();
    }
    return MonomerLibFromFilesProvider._instance!;
  }

  protected toLog(): string {
    return `MonomerLibFileManager`;
  }
  private async validateAgainstHELM(fileContent: string, fileName: string): Promise<void> {
    const isValid = this.isValidHELMLibrary(fileContent, fileName);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private isValidHELMLibrary(fileContent: string, fileName: string): boolean {
    return this.fileValidator.validateFile(fileContent, fileName);
  }

  async refreshLists(): Promise<void> {
    await Promise.all([
      this.updateValidLibList(),
      this.updateValidSetList(),
    ]);
  }


  async addOrUpdateLibraryString(fileName: string, contentString: string): Promise<void> {
    fileName = fileName.endsWith('.json') ? fileName : `${fileName}.json`;
    try {
      await this.validateAgainstHELM(contentString, fileName);
      await grok.dapi.files.writeAsText(this.LIB_PATH + `${fileName}`, contentString);
      await this.updateValidLibList();
      const fileExists = await grok.dapi.files.exists(this.LIB_PATH + `${fileName}`);
      if (!fileExists)
        grok.shell.error(`Failed to add ${fileName} library`);
      else
        grok.shell.info(`Added ${fileName} HELM library`);
    } catch (e) {
      this.logger.error(e);
      grok.shell.error(`Failed to add ${fileName} library`);
    }
  }

  async addOrUpdateLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    const contentString = JSON.stringify([...(monomers.map((m) => ({...m, wem: undefined, lib: undefined})))], null, 2);
    await this.addOrUpdateLibraryString(libraryName, contentString);
  }

  get loadPromise(): Promise<void> {
    return this.filesPromise;
  }

  async listLibraries(): Promise<string[]> {
    return this.validLibList;
  }

  async listSets(): Promise<string[]> {
    return this.validSetList;
  }

  async loadLibraries(names: string[], isFullPath: boolean = false): Promise<MonomerLibData[]> {
    let rawLibData: any[] = [];
    names = names.map((n) => {
      const normalized = n.endsWith('.json') ? n : `${n}.json`;
      return normalized;
    });
    const paths = isFullPath ? names : names.map((n) => this.LIB_PATH + n);
    const libs: MonomerLibData[] = [];

    for (let i = 0; i < paths.length; i++) {
      const formatedLib: MonomerLibData = {};
      try {
        const file = await grok.dapi.files.readAsText(paths[i]);
        rawLibData = JSON.parse(file);
        rawLibData.forEach((monomer) => {
          const polymerType = monomer[REQ.POLYMER_TYPE];
          const monomerSymbol = monomer[REQ.SYMBOL];
          if (!formatedLib[polymerType])
            formatedLib[polymerType] = {};
          formatedLib[polymerType][monomerSymbol] = monomer as Monomer;
        });
      } catch (err) {
        this.logger.error(`Failed to load monomer library file: ${paths[i]}`);
        this.logger.error(err);
        continue;
      } finally {
        libs.push(formatedLib);
      }
    }
    return libs;
  }

  async loadSets(names: string[]): Promise<IMonomerSet[]> {
    const sets: IMonomerSet[] = [];
    const lib = this.libHelper.getMonomerLib();
    for (const name of names) {
      const correctedName = name.endsWith('.json') ? name : `${name}.json`;
      const path = this.SETS_PATH + correctedName;
      try {
        const file = await grok.dapi.files.readAsText(path);
        const raw = JSON.parse(file);
        const description = raw['description'];
        const placeholders = Object.entries(raw['placeholders']).map(([k, v]: [string, any]) => {
          const placeholderSymbol = k;
          const polymerType = v['polymerType'] as PolymerType;
          const monomerType = v['monomerType'] as MonomerType;
          const monomerLinks = v['set'] as IMonomerLinkData[];

          return new MonomerSetPlaceholder(lib, placeholderSymbol, polymerType, monomerType, monomerLinks);
        });

        sets.push(new MonomerSet(description, placeholders));
      } catch (err) {
        this.logger.error(`Failed to load monomer set file: ${path}`);
        this.logger.error(err);
        continue;
      }
    }
    return sets;
  }

  async deleteLibrary(name: string): Promise<void> {
    const correctedName = name.endsWith('.json') ? name : `${name}.json`;
    await grok.dapi.files.delete(this.LIB_PATH + correctedName);
    await this.updateValidLibList();
  }

  async deleteSet(name: string): Promise<void> {
    const correctedName = name.endsWith('.json') ? name : `${name}.json`;
    await grok.dapi.files.delete(this.SETS_PATH + correctedName);
    await this.updateValidSetList();
  }

  async addOrUpdateSetString(name: string, contentString: string): Promise<void> {
    const correctedName = name.endsWith('.json') ? name : `${name}.json`;
    await grok.dapi.files.writeAsText(this.SETS_PATH + correctedName, contentString);
    await this.updateValidSetList();
  }

  async addOrUpdateSet(setName: string, monomerSet: IMonomerSet): Promise<void> {
    const contentString = JSON.stringify(monomerSet, null, 2);
    await this.addOrUpdateSetString(setName, contentString);
  }

  async deleteMonomersFromLibrary(libraryName: string, monomers: ({ polymerType: PolymerType; symbol: string; }[])) {
    const correctedName = libraryName.endsWith('.json') ? libraryName : `${libraryName}.json`;
    const path = this.LIB_PATH + correctedName;
    const file = await grok.dapi.files.readAsText(path);
    let rawLibData: Monomer[] = JSON.parse(file);
    rawLibData = rawLibData.filter((m) => !monomers.some((toDelete) =>
      m.polymerType === toDelete.polymerType && m.symbol === toDelete.symbol));
    const contentString = JSON.stringify(rawLibData, null, 2);
    await grok.dapi.files.writeAsText(path, contentString);
    await this.updateValidLibList();
  }


  private validLibList: string[] = [];
  private libListHasChanged(newList: string[]): boolean {
    const currentList = this.validLibList;
    return newList.length !== currentList.length || newList.some((el, i) => el !== currentList[i]);
  }

  private validSetList: string[] = [];
  private setListHasChanged(newList: string[]): boolean {
    const currentList = this.validSetList;
    return newList.length !== currentList.length || newList.some((el, i) => el !== currentList[i]);
  }

  async getLibraryAsString(libName: string): Promise<string> {
    const correctedName = libName.endsWith('.json') ? libName : `${libName}.json`;
    return grok.dapi.files.readAsText(this.LIB_PATH + correctedName);
  }

  async getSingleLibrary(name: string): Promise<MonomerLibData | null> {
    return (await this.loadLibraries([name]))?.[0];
  }

  async getSingleLibraryWithFullPath(path: string): Promise<MonomerLibData | null> {
    return (await this.loadLibraries([path], true))?.[0];
  }

  async updateOrAddMonomersInLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    const correctedName = libraryName.endsWith('.json') ? libraryName : `${libraryName}.json`;
    const path = this.LIB_PATH + correctedName;
    const file = await grok.dapi.files.readAsText(path);
    const rawLibData: Monomer[] = JSON.parse(file);
    for (const monomer of monomers) {
      const existingIdx = rawLibData.findIndex((m) =>
        m.polymerType === monomer.polymerType && m.symbol === monomer.symbol);
      if (existingIdx !== -1)
        rawLibData[existingIdx] = {...monomer, lib: undefined, wem: undefined};
      else
        rawLibData.push({...monomer, lib: undefined, wem: undefined});
    }
    const contentString = JSON.stringify(rawLibData, null, 2);
    await grok.dapi.files.writeAsText(path, contentString);
    await this.updateValidLibList();
  }

  private async updateValidLibList(): Promise<void> {
    const logPrefix: string = `${this.toLog()}.updateValidLibList()`;
    this.logger.debug(`${logPrefix}, start`);
    this.filesPromise = this.filesPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const invalidFiles = [] as string[];
      // console.log(`files before validation:`, this.libraryEventManager.getValidFilesPathList());
      const libPathList = await this.getLibFileListAtLocation();

      if (!this.libListHasChanged(libPathList)) {
        this.logger.debug(`${logPrefix}, end, not changed`);
        return;
      }

      for (const path of libPathList) {
        if (!path.endsWith('.json')) {
          invalidFiles.push(path);
          continue;
        }

        const fileContent = await grok.dapi.files.readAsText(this.LIB_PATH + `${path}`);
        if (!this.isValidHELMLibrary(fileContent, path))
          invalidFiles.push(path);
      }

      const validLibPathList = libPathList.filter((path) => !invalidFiles.includes(path));

      if (this.libListHasChanged(validLibPathList)) {
        this.validLibList = validLibPathList;
        await this.libHelper.loadMonomerLib(true);
      }
      // console.log(`files after validation:`, this.libraryEventManager.getValidFilesPathList());

      if (validLibPathList.some((el) => !el.endsWith('.json')))
        this.logger.warning(`Wrong validation: ${validLibPathList}`);

      if (invalidFiles.length > 0) {
        const message = `Invalid monomer library files in ${this.LIB_PATH}` +
            `, consider fixing or removing them: ${invalidFiles.join(', ')}`;

        this.logger.warning(message);
        // grok.shell.warning(message);
      }
      this.logger.debug(`${logPrefix}, OUT`);
      this._onChanged.next();
    });
    this.logger.debug(`${logPrefix}, end`);
    return this.filesPromise;
  }

  private async updateValidSetList(): Promise<void> {
    const logPrefix: string = `${this.toLog()}.updateValidSetList()`;
    this.logger.debug(`${logPrefix}, start`);
    this.filesPromise = this.filesPromise.then(async () => {
      this.logger.debug(`${logPrefix}, IN`);
      const invalidFiles = [] as string[];
      const setPathList: string[] = await this.getSetFileListAtLocation();

      if (!this.setListHasChanged(setPathList)) {
        this.logger.debug(`${logPrefix}, end, not changed`);
        return;
      }

      for (const path of setPathList) {
        if (!path.endsWith('.json')) {
          invalidFiles.push(path);
          continue;
        }

        //const fileContent = await grok.dapi.files.readAsText(this.SETS_PATH + `${path}`);
        // TODO: Validate monomer set
        // if (!this.isValidMonomerSet(fileContent, path))
        //   invalidFiles.push(path);
      }

      const validSetPathList = setPathList.filter((path) => !invalidFiles.includes(path));
      if (this.setListHasChanged(validSetPathList)) {
        this.validSetList = validSetPathList;
        this.libHelper.loadMonomerSets(true);
      }

      this.logger.debug(`${logPrefix}, OUT`);
      this._onChanged.next();
    });
    this.logger.debug(`${logPrefix}, end`);
    return this.filesPromise;
  }

  /** Get relative paths for files in LIB_PATH, SET_PATH  */
  private async getLibFileListAtLocation(): Promise<string[]> {
    const logPrefix = `${this.toLog()}.getLibFileListAtLocation()`;
    this.logger.debug(`${logPrefix}, start`);

    const libPaths = await grok.dapi.files.list(this.LIB_PATH)
      .then((l) => l.filter((f) => f.isFile).map((fi) => fi.fullPath));

    const checkForUi = false;
    const existingLibPaths: string[] = [];
    if (checkForUi) {
      // WARNING: an extra sanity check,
      // caused by unexpected behavior of grok.dapi.files.list() when it returns non-existent paths
      for (const path of libPaths) {
        const exists = await grok.dapi.files.exists(path);
        if (exists)
          existingLibPaths.push(path);
      }
    } else
      existingLibPaths.push(...libPaths);

    return existingLibPaths.map((p) => /* relative to LIB_PATH */ p.substring(this.LIB_PATH.length));
  }

  private async getSetFileListAtLocation(): Promise<string[]> {
    const logPrefix = `${this.toLog()}.getSetFileListAtLocation()`;
    this.logger.debug(`${logPrefix}, start`);

    const setPaths = await grok.dapi.files.list(this.SETS_PATH)
      .then((l) => l.map((fi) => fi.fullPath));

    const checkForUi = false;
    const existingSetPaths: string[] = [];
    if (checkForUi) {
      for (const path of setPaths) {
        const exists = await grok.dapi.files.exists(path);
        if (exists)
          existingSetPaths.push(path);
      }
    } else
      existingSetPaths.push(...setPaths);

    return existingSetPaths.map((p) => /* relative to SET_PATH */ p.substring(this.SETS_PATH.length));
  }
}
