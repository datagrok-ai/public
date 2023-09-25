/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MonomerLibHelper} from './monomer-lib';

import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {Monomer} from '@datagrok-libraries/bio/src/types';

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string> ) {
    this.helmColumn = helmColumn;
  }

  async convert(): Promise<DG.Column<string>> {
    const helmCoordinatesColumn: DG.Column<string> = await this.getHelmCoordinatesColumn();
    helmCoordinatesColumn.toList().forEach(
      (pseudoMolfile: string) => {
        this.getPolymerMolfile(pseudoMolfile);
      });
    return helmCoordinatesColumn;
  }

  async getHelmCoordinatesColumn(): Promise<DG.Column<string>> {
    const polymerMolfileCol: DG.Column<string> = await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerMolfileCol;
  }

  private getPolymerMolfile(helmCoordinates: string): string {
    const meta = new MonomerShiftMetadata(helmCoordinates);
    console.log('symbols:', meta.monomerSymbols);
    meta.monomerSymbols.forEach((symbol: string, idx: number) => {
      const shift = meta.getMonomerShifts(idx);
      console.log('shift:', shift);
    });
    return '';
  }
}

class MonomerShiftMetadata {
  constructor(helmCoordinates: string) {
    this.molfileHandler = MolfileHandler.getInstance(helmCoordinates);
  }

  private molfileHandler: MolfileHandlerBase;

  get monomerSymbols(): string[] {
    return this.molfileHandler.atomTypes;
  }

  getMonomerShifts(monomerIdx: number): {x: number, y: number} {
    const x = this.molfileHandler.x[monomerIdx];
    const y = this.molfileHandler.y[monomerIdx];
    return {x, y};
  }
}

class MonomerWrapper {
  constructor(monomerSymbol:string) {
    const monomerLib = MonomerLibHelper.instance.getBioLib();
    const monomer = monomerLib.getMonomer(HELM_POLYMER_TYPE.PEPTIDE, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.monomer = monomer;
  }

  private monomer: Monomer;
}

class MolfileCompiler {
  constructor() { }
}
