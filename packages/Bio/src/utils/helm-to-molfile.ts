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
    meta.monomerSymbols.forEach((monomerSymbol: string, idx: number) => {
      const monomer = new MonomerWrapper(monomerSymbol, meta.getMonomerShifts(idx));
      const molfile = monomer.getMonomerMolfile(meta.getMonomerShifts(idx), 0);
      console.log(`molfile for ${monomerSymbol}:`, molfile);
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
  constructor(monomerSymbol:string, private shift: {x: number, y: number}) {
    const monomerLib = MonomerLibHelper.instance.getBioLib();
    const monomer = monomerLib.getMonomer(HELM_POLYMER_TYPE.PEPTIDE, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.monomer = monomer;
  }

  private monomer: Monomer;

  getMonomerMolfile(shift: {x: number, y: number}, angle: number): string {
    const molfile = this.monomer.molfile;
    return molfile;
  }
}

class MolfileWrapper {
  constructor(private molfileV2K: string) {
    this.molfileV2K = molfileV2K;
    const lines = this.molfileV2K.split('\n');
    const atomCountIdx = {begin: 0, end: 3};
    const bondCountIdx = {begin: 3, end: 6};
    const countsLineIdx = 3;
    const atomBlockIdx = 4;
    const atomCount = parseInt(lines[countsLineIdx].substring(atomCountIdx.begin, atomCountIdx.end));
    const bondCount = parseInt(lines[countsLineIdx].substring(bondCountIdx.begin, bondCountIdx.end));
    this.atomLines = lines.slice(atomBlockIdx, atomBlockIdx + atomCount);
    this.bondLines = lines.slice(atomBlockIdx + atomCount, atomBlockIdx + atomCount + bondCount);
  }

  atomLines: string[];
  bondLines: string[];

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.atomLines = this.atomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      const newX = x + shift.x;
      const newY = y + shift.y;
      const newLine = line.replace(x.toString(), newX.toString()).replace(y.toString(), newY.toString());
      return newLine;
    });
  }

  rotateCoordinates(angle: number): void {
    this.atomLines = this.atomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      const newX = x * Math.cos(angle) - y * Math.sin(angle);
      const newY = x * Math.sin(angle) + y * Math.cos(angle);
      const newLine = line.replace(x.toString(), newX.toString()).replace(y.toString(), newY.toString());
      return newLine;
    });
  }

  getData(): {atomLines: string[], bondLines: string[]} {
    return {atomLines: this.atomLines, bondLines: this.bondLines};
  }
}

class PolymerWrapper {
  constructor() { }

  compileToMolfile(): string {
    const molfile = '';
    return molfile;
  }
}
