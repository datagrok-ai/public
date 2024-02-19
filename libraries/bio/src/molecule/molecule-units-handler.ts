import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {BuildDataFunc, UnitsHandlerBase} from '../utils/units-handler-base';
import {buildDataMolV2000, buildDataMolV3000} from './molecule-build-data';
import {MoleculeBuildDataFunc, MoleculeBase} from './types';
import {IPdbHelper} from '../pdb/pdb-helper';

export const Temps = new class {
  uh = `units-handler.${DG.SEMTYPE.MOLECULE}`;
  data = `data.${DG.SEMTYPE.MOLECULE}`;
}();

export enum MoleculeUnits {
  molV2000 = 'molV2000',
  molV3000 = 'molV3000',
}

export const MoleculeFileExts: { [units: string]: string } = {
  [MoleculeUnits.molV2000]: 'mol',
  [MoleculeUnits.molV3000]: 'mol',
};

export const MoleculeBuildDataFuncs: { [units: string]: MoleculeBuildDataFunc } = {
  [MoleculeUnits.molV2000]: buildDataMolV2000,
  [MoleculeUnits.molV3000]: buildDataMolV3000,
};

export const MoleculeNameColumnNames: string[] = ['name', 'symbol'];

export class MoleculeUnitsHandler extends UnitsHandlerBase<string, MoleculeBase> {
  constructor(col: DG.Column) {
    super(col, DG.SEMTYPE.MOLECULE, MoleculeNameColumnNames);
  }

  protected override getFileExt(): string {
    return MoleculeFileExts[this.units];
  }

  protected override getBuildDataFunc(): BuildDataFunc<string, MoleculeBase> {
    return MoleculeBuildDataFuncs[this.units];
  }

  public async getAsPdb(ph: IPdbHelper): Promise<DG.Column<string>> {
    if (this.units !== 'mol')
      throw new Error(`Unsupported units '${this.units}' of the molecule column, must be 'mol'.`);
    const colLength: number = this.column.length;
    const resCol: DG.Column = DG.Column.fromType(DG.TYPE.STRING, this.column.name, colLength);
    for (let rowI = 0; rowI < colLength; ++rowI) {
      const molVal = this.column.get(rowI);
      const pdbVal = molVal === null ? null : await ph.molToPdb(molVal);
      resCol.set(rowI, pdbVal);
    }
    return resCol;
  }

  public static getOrCreate(col: DG.Column): MoleculeUnitsHandler {
    let res = col.temp[Temps.uh];
    if (!res)
      res = col.temp[Temps.uh] = new MoleculeUnitsHandler(col);
    return res;
  }
}
