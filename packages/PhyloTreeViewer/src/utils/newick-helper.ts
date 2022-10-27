import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';
import {newickToDf} from './index';


export class NewickHelper implements bio.INewickHelper {

  newickToDf(newick: string, name: string): DG.DataFrame {
    return newickToDf(newick, name);
  }
}
