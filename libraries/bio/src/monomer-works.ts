import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SeqPalette} from './seq-palettes';

export class MonomerWorks {
  monomerLibs: string[] | null = null;
  acceptedMonomers: Set<string> | null = null;

  constructor() {

  }

  /*
  getColors(seqCol: DG.Column<string>, method?: string): SeqPalette {
    if (!method) {
      // return Unknown
      if () {} //monomer not in acceptedMonomers color it red

    } else {
      switch (method) {
        case :
          break;
        default:
          throw new Error('Not implemented');
      }
    }

    return ;
  }
  */
}
