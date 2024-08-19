import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {BuildDataFunc, DataBase} from '../utils/units-handler-base';

export class MoleculeBase extends DataBase {
  constructor(name?: string) {
    super(name);
  }
}

export type MoleculeBuildDataFunc = BuildDataFunc<string, MoleculeBase>;
