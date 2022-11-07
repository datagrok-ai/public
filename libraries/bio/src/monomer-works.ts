import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from './types';



export class MonomerWorks {
  monomerLib: IMonomerLib;

  constructor(monomerLib: IMonomerLib) {
    this.monomerLib = monomerLib;
  }

  //types according to Monomer possible
  public getCappedMonomer(monomerType: string, monomerName: string): Monomer | null {
    //TODO
    return null;
  }
}
