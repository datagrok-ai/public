import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Molecule3DBase, Molecule3DBuildDataFunc} from './types';
import {Molecule3DData} from '../viewers/molecule3d';

export const buildDataPdb: Molecule3DBuildDataFunc = (src: any): Molecule3DBase => {
  return new Molecule3DBase(undefined);
};

export const buildDataPdbqt: Molecule3DBuildDataFunc = (src: any): Molecule3DBase => {
  return new Molecule3DBase(undefined);
};

export const buildDataMmcif: Molecule3DBuildDataFunc = (src: any): Molecule3DBase => {
  return new Molecule3DBase(undefined);
};
