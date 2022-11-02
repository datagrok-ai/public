import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer, MonomerType} from './types';
import {capTheMonomer} from './utils/to-atomic-level';


export class MonomerWorks {
  monomerLib: IMonomerLib;

  constructor(monomerLib: IMonomerLib) {
    this.monomerLib = monomerLib;
  }

  //types according to Monomer possible
  public getCappedMonomer(monomerType: MonomerType, monomerName: string): Monomer | null {
    if (!this.monomerLib.types.includes(monomerType))
      throw new Error(`MonomerWorksError: monomer type \`${monomerType}\` is incorrect`);

    //return capTheMonomer(this.monomerLib!.get(type, name));
    return this.monomerLib.get(monomerType, monomerName);
  }
}
