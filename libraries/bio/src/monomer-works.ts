import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib} from './types/index';
import {capTheMonomer} from './utils/to-atomic-level';
import {getMonomerLib} from './utils/monomer-lib';


export class MonomerWorks {
  monomerLib: IMonomerLib | null = null;
  //Forbid non sequences
  sequenceCol: DG.Column | null = null;

  public async init(sequenceCol: DG.Column | null): Promise<void> {
    const funcList: DG.Func[] = DG.Func.find({package: 'Helm', name: 'getAllLibsData'});
    if (funcList.length === 0)
      await grok.functions.call('Helm:initHelm');

    this.monomerLib = await getMonomerLib();
    this.monomerLib.onChanged.subscribe(() => {
      // invalidate something
    });
    this.sequenceCol = sequenceCol;
  }

  //types according to Monomer possible
  public getCappedMonomer(name: string, type: string): string {
    const types = Object.keys(this.monomerLib!);
    if (!types.includes(type))
      throw '';

    //return capTheMonomer(this.monomerLib!.get(type, name));
    return this.monomerLib!.get(type, name)!.m;
  }
}
