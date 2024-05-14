import * as DG from 'datagrok-api/dg';

import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';

import {APP_NAME} from '../view/const';
import {DEFAULT_LIB_FILENAME, LIB_PATH} from './data-loader/const';
import {tryCatch} from './helpers';

export class OligoToolkitPackage extends DG.Package {
  private _monomerLib?: IMonomerLib;

  get monomerLib(): IMonomerLib {
    if (!this._monomerLib)
      throw new Error('Monomer lib not loaded');
    return this._monomerLib!;
  }

  async initMonomerLib(): Promise<void> {
    if (this._monomerLib !== undefined)
      return;

    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
      `Initializing ${APP_NAME.COMBINED} monomer library ...`);
    await tryCatch(async () => {
      const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
      this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
    }, () => pi.close());
  }
}
