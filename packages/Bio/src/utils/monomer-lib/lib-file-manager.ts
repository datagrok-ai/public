/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {MonomerLib} from './monomer-lib';
import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {PolyToolMonomerLibHandler} from '../poly-tool/monomer-lib-handler';
import {
  createJsonMonomerLibFromSdf,
} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {HELM_REQUIRED_FIELDS as REQ} from '@datagrok-libraries/bio/src/utils/const';


/** Singleton for adding, validation and reading of monomer library files.
 * All files stored under LIB_PATH directory must be in .json format and satisfy HELM standard.
 * All libraries  that are not in .json format or do not satisfy HELM standard are considered as custom libraries.
 * They must be brought to the standard format before adding to LIB_PATH. */
export class MonomerLibFileManager {
  private constructor() { }

  private static instance: MonomerLibFileManager | undefined;

  static async getInstance(): Promise<MonomerLibFileManager> {
    if (MonomerLibFileManager.instance === undefined) {
      MonomerLibFileManager.instance = new MonomerLibFileManager();
      await MonomerLibFileManager.instance.init();
    }
    return MonomerLibFileManager.instance;
  }

  private async init(): Promise<void> {
    this.validateAllFiles();
  }

  /** Add standard .json monomer library  */
  async addStandardLibFile(fileContent: string, fileName: string): Promise<void> {
    await this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Added ${fileName} HELM library`);
  }

  /** Transform non-standad monomer librarieies to standard format */
  async addCustomLibFile(fileContent: string, fileName: string): Promise<void> {
    await this.validate(fileContent, fileName);
    await grok.dapi.files.writeAsText(LIB_PATH + `${fileName}`, fileContent);
    grok.shell.info(`Added ${fileName} HELM library`);
  }

  async deleteLibFile(fileName: string): Promise<void> {
    grok.dapi.files.delete(LIB_PATH + `${fileName}`);
    grok.shell.warning(`Deleted ${fileName} library`);
    await this.validateAllFiles();
  }

  async readLibraryFile(path: string, fileName: string): Promise<IMonomerLib> {
    let rawLibData: any[] = [];
    let file;
    let dfSdf;
    const fileSource = new DG.FileSource(path);
    if (fileName.endsWith('.sdf')) {
      const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'importSdf'});
      if (funcList.length === 1) {
        file = await fileSource.readAsBytes(fileName);
        dfSdf = await grok.functions.call('Chem:importSdf', {bytes: file});
        rawLibData = createJsonMonomerLibFromSdf(dfSdf[0]);
      } else {
        grok.shell.warning('Chem package is not installed');
      }
    } else if (fileName.endsWith('.json')) {
      const file = await fileSource.readAsText(fileName);
      rawLibData = JSON.parse(file);
    } else if (fileName.endsWith('.csv')) {
      // todo: replace by DataFrame's method after update of js-api
      function toJson(df: DG.DataFrame): any[] {
        return Array.from({length: df.rowCount}, (_, idx) =>
          df.columns.names().reduce((entry: { [key: string]: any }, colName) => {
            entry[colName] = df.get(colName, idx);
            return entry;
          }, {})
        );
      }

      const df = await fileSource.readCsv(fileName);
      const json = toJson(df);
      const polyToolMonomerLib = new PolyToolMonomerLibHandler(json);
      if (polyToolMonomerLib.isValid())
        rawLibData = polyToolMonomerLib.getJsonMonomerLib();
      else
        throw new Error('Invalid format of CSV monomer lib');
    } else {
      throw new Error('Monomer library of unknown file format, supported formats: SDF, JSON, CSV');
    }

    const monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
    const polymerTypes: string[] = [];
    rawLibData.forEach((monomer) => {
      if (!polymerTypes.includes(monomer[REQ.POLYMER_TYPE])) {
        monomers[monomer[REQ.POLYMER_TYPE]] = {};
        polymerTypes.push(monomer[REQ.POLYMER_TYPE]);
      }
      monomers[monomer[REQ.POLYMER_TYPE]][monomer[REQ.SYMBOL]] = monomer as Monomer;
    });

    return new MonomerLib(monomers);
  }

  private async validate(fileContent: string, fileName: string): Promise<void> {
    await this.validateAllFiles();
    const isValid = this.isValid(fileContent);
    if (!isValid)
      throw new Error(`File ${fileName} does not satisfy HELM standard`);
  }

  private async validateAllFiles(): Promise<void> {
    const list = await grok.dapi.files.list(LIB_PATH);
    const invalidFiles: string[] = [];
    for (const file of list) {
      const fileContent = await grok.dapi.files.readAsText(LIB_PATH + `${file.name}`);
      if (!this.isValid(fileContent))
        invalidFiles.push(file.name);
    }

    if (invalidFiles.length > 0)
      grok.shell.warning(`The following files in ${LIB_PATH} do not satisfy HELM standard}, delete or fix them: ${invalidFiles.join(', ')}`);
  }

  private isValid(fileContent: string): boolean {
    if (fileContent.length > 0)
      return true;
    return false;
  }
}
