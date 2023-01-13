import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PdbHelper} from '../utils/pdb-helper';
import {IPdbHelper} from '@datagrok-libraries/bio';
import {getPdbHelper} from '../package';

export class PdbApp {

  private _funcName: string = '';

  constructor() {}

  async init(df?: DG.DataFrame, funcName: string = 'pdbApp'): Promise<void> {
    this._funcName = funcName;
    await this.loadData(df);
  }

  async loadData(df?: DG.DataFrame): Promise<void> {
    const ph = new PdbHelper();
    if (!df) {
      const ph: IPdbHelper = getPdbHelper();

    }
  }
}