import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '../types';



export class MonomerWorks {
  monomerLib: IMonomerLib;

  constructor(monomerLib: IMonomerLib) {
    this.monomerLib = monomerLib;
  }

  //types according to Monomer possible
  public getCappedRotatedMonomer(monomerType: string, monomerName: string): string | null {
    const monomer = this.monomerLib.getMonomer(monomerType, monomerName);
    if (monomer)
      return monomer.molfile; //TODO cap 

    return null;
  }
}
