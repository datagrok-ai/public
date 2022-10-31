import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from './types';
import {capTheMonomer} from './utils/to-atomic-level';


export class MonomerWorks {
  monomerLib: IMonomerLib;

  constructor(monomerLib: IMonomerLib) {
    this.monomerLib = monomerLib;
  }

  //types according to Monomer possible
  public getCappedMonomer(monomerType: string, monomerName: string): Monomer | null {
    const types = Object.keys(this.monomerLib!);
    if (!types.includes(monomerType))
      throw '';

    //return capTheMonomer(this.monomerLib!.get(type, name));
    return this.monomerLib.get(monomerType, monomerName);
  }
}
