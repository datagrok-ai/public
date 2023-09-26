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

  async convertToMolfile(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> = await this.getPolymerGraphColumn();
    const molfileList = polymerGraphColumn.toList().map(
      (pseudoMolfile: string, idx: number) => {
        const sourceHelm = this.helmColumn.get(idx);
        if (!sourceHelm)
          throw new Error(`HELM column is empty at row ${idx}`);
        return this.getPolymerMolfile(sourceHelm, pseudoMolfile);
      });
    console.log('molfile list:', molfileList);
    const molfileColumn = DG.Column.fromList('string', 'molfiles', molfileList);
    return molfileColumn;
  }

  async getPolymerGraphColumn(): Promise<DG.Column<string>> {
    const polymerPseudoMolfileCol: DG.Column<string> =
      await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerPseudoMolfileCol;
  }

  private getPolymerMolfile(sourceHelm: string, polymerGraph: string): string {
    const meta = new MonomerShiftMetadata(polymerGraph);
    const polymerWrapper = new PolymerWrapper(sourceHelm, polymerGraph);
    meta.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = meta.getMonomerShifts(monomerIdx);
      polymerWrapper.addShiftedMonomer(monomerSymbol, shift);
    });
    const polymerMolfile = polymerWrapper.compileToMolfile();
    return polymerMolfile;
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
    this.molfileWrapper = new MolfileWrapper(monomer.molfile);
  }

  private monomer: Monomer;
  private molfileWrapper: MolfileWrapper;

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.molfileWrapper.shiftCoordinates(shift);
  }

  getMolfileData(bondShift: number): {atomLines: string[], bondLines: string[]} {
    return this.molfileWrapper.getData(bondShift);
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
      // get new coords with fixed precision 4 and padded with zeros
      const newX = (x + shift.x).toFixed(4);
      const newY = (y + shift.y).toFixed(4);
      const newLine = `${newX.toString().padStart(10, ' ')}${newY.toString().padStart(10, ' ')}${line.substring(20)}`;
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

  getRGroupIndices(): number[] {
    return (MolfileHandler.getInstance(this.molfileV2K))
      .atomTypes.map((atomType, idx) => {
        if (atomType === 'R#')
          return idx;
      }).filter((idx) => idx !== undefined) as number[];
  }

  /** Get R-group ids from Atom Alias block and RGP block  */
  // getRGroupIndexById(id: number): number {
  // }

  shiftBondNumbers(shift: number): void {
    this.bondLines = this.bondLines.map((line: string) => {
      const bond1 = parseInt(line.substring(0, 3));
      const bond2 = parseInt(line.substring(3, 6));
      const newBond1 = bond1 + shift;
      const newBond2 = bond2 + shift;
      return `${newBond1.toString().padStart(3, ' ')}${newBond2.toString().padStart(3, ' ')}${line.substring(6)}`;
    });
  }

  getData(bondShift: number): {atomLines: string[], bondLines: string[]} {
    this.shiftBondNumbers(bondShift);
    return {atomLines: this.atomLines, bondLines: this.bondLines};
  }
}

class PolymerWrapper {
  constructor(private sourceHelm: string, private polymerGraph: string) { }

  private monomerWrappers: MonomerWrapper[] = [];

  addShiftedMonomer(monomerSymbol: string, shift: {x: number, y: number}): void {
    const monomerWrapper = new MonomerWrapper(monomerSymbol);
    monomerWrapper.shiftCoordinates(shift);
    this.monomerWrappers.push(monomerWrapper);
  }

  compileToMolfile(): string {
    const molfileHeader = '\nDatagrok\n';
    const atomLines: string[] = [];
    const bondLines: string[] = [];
    this.monomerWrappers.forEach((monomerWrapper: MonomerWrapper) => {
      const bondShift = atomLines.length;
      const {atomLines: atomLinesPart, bondLines: bondLinesPart} = monomerWrapper.getMolfileData(bondShift);
      atomLines.push(...atomLinesPart);
      bondLines.push(...bondLinesPart);
    });
    const atomCount = atomLines.length;
    if (atomCount > 999)
      throw new Error(`Atom count is ${atomCount} and exceeds 999`);
    const bondCount = bondLines.length;
    const countsLine = `${
      atomCount.toString().padStart(3, ' ')
    }${
      bondCount.toString().padStart(3, ' ')
    }  0  0  1  0              0 V2000`;
    const molfileEnd = 'M  END\n';
    const newLineChar = '\n';
    const blockList = [molfileHeader, countsLine, atomLines.join(newLineChar), bondLines.join(newLineChar), molfileEnd];
    const molfile = blockList.join(newLineChar);
    return molfile;
  }
}
