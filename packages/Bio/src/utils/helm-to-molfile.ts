/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MonomerLibHelper} from './monomer-lib';

import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';

const HELM_SECTION_SEPARATOR = '$';
const HELM_ITEM_SEPARATOR = '|';
const enum HELM_MONOMER_TYPE {
  BACKBONE,
  BRANCH,
}

type Bond = {
  /** Global (for complex polymer) or local (for simple polymer) monomer id  */
  monomerId: number,
  rGroupId: number
}

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string>) {
    this.helmColumn = helmColumn;
  }

  async convertToMolfileColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> = await this.getPolymerGraphColumn();
    const molfileList = polymerGraphColumn.toList().map(
      (pseudoMolfile: string, idx: number) => {
        const helm = this.helmColumn.get(idx);
        if (!helm)
          return '';
        return this.getPolymerMolfile(helm, pseudoMolfile);
      });
    const molfileColName = this.helmColumn.dataFrame.columns.getUnusedName(`molfile(${this.helmColumn.name})`);
    const molfileColumn = DG.Column.fromList('string', molfileColName, molfileList);
    return molfileColumn;
  }

  private async getPolymerGraphColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> =
      await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerGraphColumn;
  }

  private getPolymerMolfile(helm: string, polymerGraph: string): string {
    const meta = new MonomerPositionDataManager(polymerGraph);
    const polymer = new Polymer(helm);
    meta.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = meta.getMonomerShifts(monomerIdx);
      polymer.addShiftedMonomer(monomerSymbol, shift);
    });
    const polymerMolfile = polymer.compileToMolfile();
    return polymerMolfile;
  }
}

class MonomerPositionDataManager {
  constructor(helmCoordinatesPseudoMolfile: string) {
    this.molfileHandler = MolfileHandler.getInstance(helmCoordinatesPseudoMolfile);
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
  constructor(
    private monomerSymbol: string,
    private polymerType: HELM_POLYMER_TYPE,
    // id: {polymer: number, monomer: number}
  ) {
    const monomerLib = MonomerLibHelper.instance.getBioLib();
    const monomer = monomerLib.getMonomer(polymerType, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.molfileWrapper = new MolfileWrapper(monomer.molfile, monomerSymbol);
  }

  private molfileWrapper: MolfileWrapper;

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.molfileWrapper.shiftCoordinates(shift);
  }

  getMolfileAtomAndBondLines(bondShift: number): {atomLines: string[], bondLines: string[]} {
    return this.molfileWrapper.getAtomAndBondLines(bondShift);
  }
}

class MolfileWrapper {
  constructor(private molfileV2K: string, private monomerSymbol: string) {
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

    this.rgpIdToAtomicIdxMap = this.getRGroupIdToAtomicIdxMap(lines);
    this.shiftMonomerToDefaultPosition();
  }

  private atomLines: string[];
  private bondLines: string[];
  private rgpIdToAtomicIdxMap: {[key: number]: number} = {};

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.atomLines = this.atomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      const newX = x + shift.x;
      const newY = y + shift.y;
      if (isNaN(newX) || isNaN(newY))
        throw new Error(`Cannot shift coordinates by ${shift.x}, ${shift.y}`);
      const newLine = `${newX.toFixed(4).padStart(10, ' ')}${newY.toFixed(4).padStart(10, ' ')}${line.substring(20)}`;
      return newLine;
    });
  }

  rotateCoordinates(angle: number): void {
    this.atomLines = this.atomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      const newX = x * Math.cos(angle) - y * Math.sin(angle);
      const newY = x * Math.sin(angle) + y * Math.cos(angle);
      const newLine = `${newX.toFixed(4).padStart(10, ' ')}${newY.toFixed(4).padStart(10, ' ')}${line.substring(20)}`;
      return newLine;
    });
  }

  shiftBondNumbers(shift: number): void {
    this.bondLines = this.bondLines.map((line: string) => {
      const bond1 = parseInt(line.substring(0, 3));
      const bond2 = parseInt(line.substring(3, 6));
      const newBond1 = bond1 + shift;
      const newBond2 = bond2 + shift;
      return `${newBond1.toString().padStart(3, ' ')}${newBond2.toString().padStart(3, ' ')}${line.substring(6)}`;
    });
  }

  getAtomAndBondLines(bondShift: number): {atomLines: string[], bondLines: string[]} {
    this.shiftBondNumbers(bondShift);
    return {atomLines: this.atomLines, bondLines: this.bondLines};
  }

  getAtomCoordinatesByIdx(atomicIdx: number): {x: number, y: number} {
    const molfileHandler = MolfileHandler.getInstance(this.molfileV2K);
    const x = molfileHandler.x[atomicIdx];
    const y = molfileHandler.y[atomicIdx];
    return {x, y};
  }

  private getRGroupAtomicIndices(): number[] {
    return (MolfileHandler.getInstance(this.molfileV2K))
      .atomTypes.map((atomType, idx) => {
        if (atomType === 'R#')
          return idx;
      }).filter((idx) => idx !== undefined) as number[];
  }

  private getRGroupIdToAtomicIdxMap(lines: string[]): {[key: number]: number} {
    const rgroupLines = lines.filter((line: string) => line.startsWith('M  RGP'));

    const map: {[key: number]: number} = {};
    rgroupLines.forEach((line: string) => {
      const indices = line.split(/\s+/).filter((item) => item)
        .slice(3).map((item) => parseInt(item));
      const atomIdxToRgpIdxList = new Array<[number, number]>(indices.length / 2);
      for (let i = 0; i < indices.length; i += 2)
        atomIdxToRgpIdxList[i / 2] = [indices[i + 1], indices[i] - 1];
      const mapPart = Object.fromEntries(atomIdxToRgpIdxList);
      Object.assign(map, mapPart);
    });

    const atomAliasLinesIndices = lines.map((line: string, idx: number) => {
      if (line.startsWith('A  '))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
    const atomAliasLines = atomAliasLinesIndices.map((idx) => lines[idx]);
    const atomAliasTextLines = atomAliasLinesIndices.map((idx) => lines[idx + 1]);
    atomAliasLines.forEach((line: string, idx: number) => {
      const rgpAtomIdx = parseInt(line.split(/\s+/)[1]) - 1;
      const rgpId = parseInt(atomAliasTextLines[idx].substring(1));
      map[rgpId] = rgpAtomIdx;
    });

    const rGroupAtomicIndices = this.getRGroupAtomicIndices();
    const unaccounted = rGroupAtomicIndices.filter((idx) => !Object.values(map).includes(idx));
    if (unaccounted.length !== 0)
      throw new Error(`Unaccounted R group indices: ${unaccounted}`);

    return map;
  }

  getRGroupAtomicIdxById(rgroupId: number): number {
    return this.rgpIdToAtomicIdxMap[rgroupId];
  }

  getAttachmentAtomByRGroupId(rgroupId: number): number {
    const rGroupAtomicIdx = this.getRGroupAtomicIdxById(rgroupId);
    const molfileHandler = MolfileHandler.getInstance(this.molfileV2K);
    const allBondedPairs = molfileHandler.pairsOfBondedAtoms;
    const bondedPairWithRGP = allBondedPairs.find((pair) => pair.includes(rGroupAtomicIdx));
    if (!bondedPairWithRGP)
      throw new Error(`Cannot find bonded pair with RGP ${rgroupId} for monomer ${this.monomerSymbol}`);
    const attachmentIdx = bondedPairWithRGP.find((idx) => idx !== rGroupAtomicIdx);
    if (!attachmentIdx)
      throw new Error(`Cannot find attachment atom for RGP ${rgroupId} for monomer ${this.monomerSymbol}`);
    return attachmentIdx;
  }

  private shiftR1GroupToOrigin(): void {
    const r1Idx = this.getRGroupAtomicIdxById(1);
    const {x, y} = this.getAtomCoordinatesByIdx(r1Idx);
    this.shiftCoordinates({x: -x, y: -y});
  }

  private alignR2AlongX(): void {
    const r2Coordinates = this.getAtomCoordinatesByIdx(this.getRGroupAtomicIdxById(2));
    const tan = r2Coordinates.y / r2Coordinates.x;
    const angle = Math.atan(tan);
    if (isNaN(angle))
      throw new Error(`Cannot calculate angle for R2 group for monomer ${this.monomerSymbol}`);
    this.rotateCoordinates(-angle);
  }

  private shiftMonomerToDefaultPosition(): void {
    this.shiftR1GroupToOrigin();
    if (this.rgpIdToAtomicIdxMap[2])
      this.alignR2AlongX();
  }
}

class SimplePolymer {
  constructor(private simplePolymer: string) {
    this.polymerType = this.getPolymerType();
    this.id = this.getId();
    const {monomers, monomerTypes} = this.extractMonomerSymbolsAndTypes();
    this.monomers = monomers;
    this.monomerTypes = monomerTypes;
  }

  polymerType: string;
  id: number;
  monomers: string[];
  monomerTypes: HELM_MONOMER_TYPE[];

  get fullId(): string {
    return this.polymerType + this.id.toString();
  }

  private getPolymerType(): string {
    const regex = new RegExp(
      `(${HELM_POLYMER_TYPE.PEPTIDE}|${HELM_POLYMER_TYPE.RNA})[0-9]+{`
    );
    const match = this.simplePolymer.match(regex);
    if (!match)
      throw new Error(`Unsupported polymer type in ${this.simplePolymer}`);
    const polymerType = match[1];
    return polymerType;
  }

  private getId(): number {
    const regex = new RegExp(`${this.polymerType}([0-9]+){`);
    const match = this.simplePolymer.match(regex);
    if (!match)
      throw new Error(`Cannot parse simple polymer id from ${this.simplePolymer}`);
    const id = parseInt(match[1]);
    return id;
  }

  private extractMonomerSymbolsAndTypes(): {monomers: string[], monomerTypes: HELM_MONOMER_TYPE[]} {
    const helmWrapperRegex = new RegExp(`${this.polymerType}${this.id}{|}`, 'g');
    const monomerGroups = this.simplePolymer.replace(helmWrapperRegex, '').split('.');
    const monomerList: string[] = [];
    const monomerTypeList: HELM_MONOMER_TYPE[] = [];
    monomerGroups.forEach((monomerGroup) => {
      const splitted = monomerGroup.split(/\(|\)/)
        .map((el) => el.replace(/[\[\]]/g, ''));
      monomerList.push(...splitted);
      // WARNING: only the groups of the form r(A)p, as in RNA, are supported
      const monomerTypes = splitted.map(
        (_, idx) => (idx % 2 === 0) ? HELM_MONOMER_TYPE.BACKBONE : HELM_MONOMER_TYPE.BRANCH
      );
      monomerTypeList.push(...monomerTypes);
    });
    return {monomers: monomerList, monomerTypes: monomerTypeList};
  }

  getBondData(): Bond[][] {
    const result: Bond[][] = [];
    const backboneMonomerIndices = this.monomerTypes.map((type, idx) => {
      if (type === HELM_MONOMER_TYPE.BACKBONE)
        return idx;
    }
    ).filter((idx) => idx !== undefined) as number[];
    const branchMonomerIndices = this.monomerTypes.map((type, idx) => {
      if (type === HELM_MONOMER_TYPE.BRANCH)
        return idx;
    }
    ).filter((idx) => idx !== undefined) as number[];
    for (let i = 0; i < backboneMonomerIndices.length - 1; i++) {
      const backboneIdx = backboneMonomerIndices[i];
      const nextBackboneIdx = backboneMonomerIndices[i + 1];
      result.push([{monomerId: backboneIdx, rGroupId: 2}, {monomerId: nextBackboneIdx, rGroupId: 1}]);
    }
    for (let i = 0; i < branchMonomerIndices.length; i++) {
      const branchIdx = branchMonomerIndices[i];
      const backboneIdx = branchIdx - 1;
      result.push([{monomerId: backboneIdx, rGroupId: 3}, {monomerId: branchIdx, rGroupId: 1}]);
    }
    return result;
  }
}

class ConnectionList {
  constructor(connectionList: string) {
    const splitted = connectionList.split(HELM_ITEM_SEPARATOR);
    splitted.forEach((connectionItem: string) => this.validateConnectionItem(connectionItem));
    this.connectionItems = splitted;
    console.log(`${connectionList}:`, this.getConnectionData());
  }

  private connectionItems: string[];

  private validateConnectionItem(connectionItem: string): void {
    const allowedType = `(${HELM_POLYMER_TYPE.PEPTIDE}|${HELM_POLYMER_TYPE.RNA})`;
    const regex = new RegExp(`${allowedType}[0-9]+,${allowedType}[0-9]+,[0-9]+:R[0-9]+-[0-9]+:R[0-9]+`, 'g');
    if (!connectionItem.match(regex))
      throw new Error(`Cannot parse connection item from ${connectionItem}`);
  }

  getConnectionData(): {polymerId: string, bond: Bond}[][] {
    const result: {polymerId: string, bond: Bond}[][] = [];
    this.connectionItems.forEach((connectionItem: string) => {
      const pair: {polymerId: string, bond: Bond}[] = [];
      const splitted = connectionItem.split(',');
      splitted[2].split('-').forEach((item, idx) => {
        const polymerId = splitted[idx];
        const data = item.split(':');
        const monomerId = parseInt(data[0]);
        const rGroupId = parseInt(data[1].slice(1));
        const bondData = {monomerId, rGroupId};
        pair.push({polymerId, bond: bondData});
      });
      result.push(pair);
    });
    return result;
  }
}

class Helm {
  constructor(private helm: string) {
    const helmSections = this.helm.split(HELM_SECTION_SEPARATOR);
    const simplePolymers = helmSections[0].split(HELM_ITEM_SEPARATOR);
    this.simplePolymers = simplePolymers
      .map((item) => new SimplePolymer(item));
    if (helmSections[1] !== '')
      this.connectionList = new ConnectionList(helmSections[1]);
    this.bondData = this.getBondData();
    console.log(`${this.helm}:`, this.bondData);
  }

  private bondData: Bond[][];

  toString() {
    return this.helm;
  }

  get polymerType(): HELM_POLYMER_TYPE {
    const polymerType = this.simplePolymers[0].polymerType;
    return polymerType as HELM_POLYMER_TYPE;
  }

  private simplePolymers: SimplePolymer[];
  private connectionList?: ConnectionList;

  private shiftBondMonomerIds(shift: number, bonds: Bond[][]): void {
    bonds.forEach((bond) => {
      bond.forEach((bondPart) => {
        bondPart.monomerId += shift;
      });
    });
  }

  private getMonomerIdShifts(): {[simplePolymerId: string]: number} {
    const result: {[simplePolymerId: string]: number} = {};
    let shift = 0;
    this.simplePolymers.forEach((simplePolymer) => {
      result[simplePolymer.fullId] = shift;
      shift += simplePolymer.monomers.length;
    });
    console.log(`${this.helm}:`, result);
    return result;
  }

  private getBondData(): Bond[][] {
    const shifts = this.getMonomerIdShifts();
    const result: Bond[][] = [];
    this.simplePolymers.forEach((simplePolymer) => {
      const bondData = simplePolymer.getBondData();
      const shift = shifts[simplePolymer.fullId];
      this.shiftBondMonomerIds(shift, bondData);
      result.push(...bondData);
    });
    if (this.connectionList) {
      const connectionData = this.connectionList.getConnectionData();
      connectionData.forEach((connection) => {
        const data: Bond[] = [];
        connection.forEach((connectionItem) => {
          const shift = shifts[connectionItem.polymerId];
          const bond = connectionItem.bond;
          bond.monomerId += shift;
          data.push(bond);
        });
        result.push(data);
        console.log(`${connection} added pair:`, data);
      });
    }
    return result;
  }
}

class Polymer {
  constructor(helm: string) {
    this.helm = new Helm(helm);
    this.polymerType = this.helm.polymerType;
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;
  private polymerType: HELM_POLYMER_TYPE;

  addShiftedMonomer(
    monomerSymbol: string,
    shift: {x: number, y: number},
    // id: {polymer: number, monomer: number}
  ): void {
    const monomerWrapper = new MonomerWrapper(monomerSymbol, this.polymerType);
    monomerWrapper.shiftCoordinates(shift);
    this.monomerWrappers.push(monomerWrapper);
  }

  compileToMolfile(): string {
    const molfileHeader = '\nDatagrok\n';
    const atomLines: string[] = [];
    const bondLines: string[] = [];
    this.monomerWrappers.forEach((monomerWrapper: MonomerWrapper) => {
      const bondShift = atomLines.length;
      const {atomLines: atomLinesPart, bondLines: bondLinesPart} = monomerWrapper.getMolfileAtomAndBondLines(bondShift);
      atomLines.push(...atomLinesPart);
      bondLines.push(...bondLinesPart);
    });
    const atomCount = atomLines.length;
    if (atomCount > 999)
      throw new Error(`Atom count in polymer ${this.helm.toString()} is ${atomCount} and exceeds 999`);
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
