import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {BuildDataFunc, UnitsHandlerBase} from '../utils/units-handler-base';
import {Molecule3DBase, Molecule3DBuildDataFunc} from './types';
import {buildDataMmcif, buildDataPdb, buildDataPdbqt} from './molecule-3d-build-data';

export const Temps = new class {
  uh = `units-handler.${DG.SEMTYPE.MOLECULE3D}`;
  data = `data.${DG.SEMTYPE.MOLECULE3D}`;
}();

export enum Molecule3DUnits {
  pdb = 'pdb',
  pdbqt = 'pdbqt',
  mmcif = 'mmcif',
}

export const Molecule3DFileExts: { [units: string]: string } = {
  [Molecule3DUnits.pdb]: 'pdb',
  [Molecule3DUnits.pdbqt]: 'pdbqt',
  [Molecule3DUnits.mmcif]: 'mmcif',
};

export const Molecule3DBuildDataFuncs: { [units: string]: Molecule3DBuildDataFunc } = {
  [Molecule3DUnits.pdb]: buildDataPdb,
  [Molecule3DUnits.pdbqt]: buildDataPdbqt,
  [Molecule3DUnits.mmcif]: buildDataMmcif,
};

export const Molecule3DNameColumnNames: string[] = ['name', 'symbol'];


export class Molecule3DUnitsHandler extends UnitsHandlerBase<any, Molecule3DBase> {
  constructor(col: DG.Column) {
    super(col, DG.SEMTYPE.MOLECULE3D, Molecule3DNameColumnNames);
  }

  protected getFileExt(): string {
    return Molecule3DFileExts[this.units];
  }

  protected getBuildDataFunc(): BuildDataFunc<any, Molecule3DBase> {
    return Molecule3DBuildDataFuncs[this.units];
  }

  public static getOrCreate(col: DG.Column): Molecule3DUnitsHandler {
    let res = col.temp[Temps.uh];
    if (!res)
      res = col.temp[Temps.uh] = new Molecule3DUnitsHandler(col);
    return res;
  }
}

