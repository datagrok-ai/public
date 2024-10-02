import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export async function getRdKitModule(): Promise<RDModule> {
  const funcList = DG.Func.find({package: 'Chem', name: 'getRdKitModule'});
  if (funcList.length === 0)
    throw new Error('Package "Chem" must be installed for getRdKitModule.');

  const res: RDModule = (await funcList[0].prepare().call()).getOutputParamValue() as RDModule;
  return res;
}
