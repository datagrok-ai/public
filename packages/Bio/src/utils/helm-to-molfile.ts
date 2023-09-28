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

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string>) {
    this.helmColumn = helmColumn;
  }

  async convertToMolfile(): Promise<DG.Column<string>> {
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

  async getPolymerGraphColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> =
      await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerGraphColumn;
  }

  private getPolymerMolfile(helm: string, polymerGraph: string): string {
    const meta = new MonomerMetadataManager(polymerGraph);
    const polymerWrapper = new PolymerWrapper(helm);
    meta.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = meta.getMonomerShifts(monomerIdx);
      polymerWrapper.addShiftedMonomer(monomerSymbol, shift);
    });
    const polymerMolfile = polymerWrapper.compileToMolfile();
    return polymerMolfile;
  }
}

class MonomerMetadataManager {
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
  constructor(
    monomerSymbol: string,
    polymerType: HELM_POLYMER_TYPE,
    // id: {polymer: number, monomer: number}
  ) {
    const monomerLib = MonomerLibHelper.instance.getBioLib();
    const monomer = monomerLib.getMonomer(polymerType, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.molfileWrapper = new MolfileWrapper(monomer.molfile);
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

    this.rgpIdToAtomicIdxMap = this.getRGroupIdToAtomicIdxMap();

    this.shiftR1GroupToOrigin();
  }

  private atomLines: string[];
  private bondLines: string[];
  private rgpIdToAtomicIdxMap: {[key: number]: number} = {};

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.atomLines = this.atomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
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

  private getRGroupIdToAtomicIdxMap(): {[key: number]: number} {
    const rgroupLines = this.molfileV2K.split('\n').filter((line: string) => line.startsWith('M  RGP'));

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
      throw new Error(`Cannot find bonded pair with RGP ${rgroupId}`);
    const attachmentIdx = bondedPairWithRGP.find((idx) => idx !== rGroupAtomicIdx);
    if (!attachmentIdx)
      throw new Error(`Cannot find attachment atom for RGP ${rgroupId}`);
    return attachmentIdx;
  }

  shiftR1GroupToOrigin(): void {
    const r1Idx = this.getRGroupAtomicIdxById(1);
    const {x, y} = this.getAtomCoordinatesByIdx(r1Idx);
    this.shiftCoordinates({x: -x, y: -y});
  }
}

class SimplePolymer {
  constructor(private simplePolymer: string) {
    this.polymerType = this.getPolymerType();
    this.id = this.getId();
  }

  polymerType: string;
  id: number;

  get fullId(): string {
    return this.polymerType + this.id.toString();
  }

  // getMonomerList(simplePolymer: string): string[] {
  //   return simplePolymer.replace(/PEPTIDE[0-9]+{|}/g, '').split('.');
  // }

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

  // list of structs: monomer #, monomer symbol, monomer type (backbone/branch)
}

class ConnectionList {
  constructor(private connectionList: string) {
    const splitted = connectionList.split(HELM_ITEM_SEPARATOR);
    splitted.forEach((connectionItem: string) => this.validateConnectionItem(connectionItem));
  }

  private validateConnectionItem(connectionItem: string): void {
    const allowedType = `(${HELM_POLYMER_TYPE.PEPTIDE}|${HELM_POLYMER_TYPE.RNA})`;
    const regex = new RegExp(`${allowedType}[0-9]+,${allowedType}[0-9]+,[0-9]+:R[0-9]+-[0-9]+:R[0-9]+`, 'g');
    if (!connectionItem.match(regex))
      throw new Error(`Cannot parse connection item from ${connectionItem}`);
  }

  // getConnectionMap() {
  //   const connectionDataList = this.connectionList.map(
  //     (connectionItem: string) => this.getConnectionData(connectionItem)
  //   );
  // }

  private getConnectionData(connectionItem: string): {[key: string]: {[key: string]: number}} {
    const splitted = connectionItem.split(',');
    const keys = ['source', 'target'];
    const data: {[key: string]: {[key: string]: number}} = {};
    splitted.forEach((item: string, idx: number) => {
      const polymerId = parseInt(item.replace('PEPTIDE', ''));
      data[keys[idx]]['polymerId'] = polymerId;
    });
    const attachmentData = splitted[2].split('-');
    keys.forEach((key: string, idx: number) => {
      const parts = attachmentData[idx].split(':');
      data[key]['monomerPosition'] = parseInt(parts[0]);
      data[key]['attachment'] = parseInt(parts[1].replace('R', ''));
    });
    return data;
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
  }

  get polymerType(): HELM_POLYMER_TYPE {
    const polymerType = this.simplePolymers[0].polymerType;
    return polymerType as HELM_POLYMER_TYPE;
  }

  private simplePolymers: SimplePolymer[];
  private connectionList?: ConnectionList;
}

class PolymerWrapper {
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
