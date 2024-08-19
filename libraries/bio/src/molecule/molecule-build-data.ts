import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {MoleculeBuildDataFunc, MoleculeBase} from './types';

export const buildDataMolV2000: MoleculeBuildDataFunc = (src: string): MoleculeBase => {
  const ma = src.match(/\@\<(?<name>)\>MOLECULE/);
  const name: string | undefined = ma ? ma.groups?.['name'] : undefined;
  return new MoleculeBase(name);
};

export const buildDataMolV3000: MoleculeBuildDataFunc = (src: string): MoleculeBase => {
  const name: string | undefined = undefined;
  return new MoleculeBase(name);
};

