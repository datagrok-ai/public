import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {expectObject} from '@datagrok-libraries/utils/src/test';

import {IMonomerLib, MonomerLibSummaryType} from '../types/index';

export function expectMonomerLib(lib: IMonomerLib, summary: { [polymerType: string]: number }): void {
  const exp: MonomerLibSummaryType = summary;
  const act: MonomerLibSummaryType = lib.getSummaryObj();

  for (const pt in exp)
    if (!exp[pt]) delete exp[pt];
  for (const pt in act)
    if (!act[pt]) delete act[pt];

  if (Object.keys(exp).length == 0 && Object.keys(act).length != 0)
    throw new Error('Expected empty monomer lib, actual is not.');
  else if (Object.keys(exp).length != 0 && Object.keys(act).length == 0)
    throw new Error('Expected non-empty monomer lib, actual is empty.');
  else {
    try {
      expectObject(act, exp);
      expectObject(exp, act);
    } catch (err) {
      throw new Error(`Expected monomer lib ${JSON.stringify(exp)} does not match actual ${JSON.stringify(act)}.`);
    }
  }
}
