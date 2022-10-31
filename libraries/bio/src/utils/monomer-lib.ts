import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Monomer, IMonomerLib} from '../types';

/** Handles Helm package presence and initialization  */
export async function getMonomerLib(): Promise<IMonomerLib> {
  // All checks for Helm package presence and is initialized
  if (DG.Func.find({package: 'Helm', name: 'getMonomerLibObj'}).length === 0)
    throw new Error('Package "Helm" must be installed for monomer libraries.');

  return (await grok.functions.call('Helm:getMonomerLibObj')) as IMonomerLib;
}