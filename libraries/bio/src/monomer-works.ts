import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { MonomerLib } from './types/index';
import { capTheMonomer } from './utils/to-atomic-level';


export class MonomerWorks {
  monomerLib: MonomerLib | null = null;
  //Forbid non sequences
  sequenceCol: DG.Column | null = null;

  public async init(sequenceCol: DG.Column | null): Promise<void> {
    const funcList: DG.Func[] = DG.Func.find({package: 'Helm', name: 'getAllLibsData'});
    if(funcList.length === 0)
      await grok.functions.call('Helm:initHelm')

    await this.refreshLib();
    this.sequenceCol = sequenceCol;
    
    grok.events.onCustomEvent('monomersLibChanged').subscribe((_) => {
      this.refreshLib();
    });
  }

  private async refreshLib(): Promise<void> {
    this.monomerLib = await grok.functions.call('Helm:getAllLibsData')
  }

  //types according to Monomer possible 
  public getCappedMonomer(name: string, type: string) : string {
    const types = Object.keys(this.monomerLib!);
    if(!types.includes(type))
      throw '';

    return capTheMonomer(this.monomerLib![type][name]);
  }
}
