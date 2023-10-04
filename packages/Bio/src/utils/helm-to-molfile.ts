/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MonomerLibHelper} from './monomer-lib';

import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';

const MAX_ATOM_COUNT = 999;

const HELM_SECTION_SEPARATOR = '$';
const HELM_ITEM_SEPARATOR = '|';
const enum HELM_MONOMER_TYPE {
  BACKBONE,
  BRANCH,
}

type Bond = {
  /** Global (for complex polymer) or local (for simple polymer) monomer index, starting from 0 */
  monomerIdx: number,
  /** RGroup id, starting from 1  */
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
      polymer.addShiftedMonomer(monomerSymbol, monomerIdx, shift);
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
    private atomNumberShift: number,
    polymerType: HELM_POLYMER_TYPE,
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

  get atomLines(): string[] {
    return this.molfileWrapper.getAtomAndBondLines().atomLines;
  }

  get bondLines(): string[] {
    const bonds = this.molfileWrapper.getShiftedBonds(this.atomNumberShift);
    console.log(`new bond numbers`, bonds);
    return this.molfileWrapper.getAtomAndBondLines(bonds).bondLines;
  }

  deleteRGroups(rGroupIds: number[]): void {
    this.molfileWrapper.removeRGroupsFromAtomBlock(rGroupIds);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, attachmentAtomIdx: number): void {
    this.molfileWrapper.replaceRGroupWithAttachmentAtom(rGroupId, attachmentAtomIdx);
  };

  getAttachmentAtomByRGroupId(rGroupId: number): number {
    const attachmentAtom = this.molfileWrapper.getAttachmentAtomByRGroupId(rGroupId);
    return attachmentAtom + this.atomNumberShift;
  }

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    this.molfileWrapper.deleteBondLineWithSpecifiedRGroup(rGroupId);
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
    this.bondedAtoms = this.getBondedAtoms();
  }

  private atomLines: string[];
  private bondLines: string[];
  private bondedAtoms: number[][] = [];
  private rgpIdToAtomicIdxMap: {[key: number]: number} = {};

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    const bondLineIdxToBeDeleted = this.bondedAtoms.map((bondedPair, idx) => {
      if (bondedPair.some((atomId) => atomId === -rGroupId))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
    this.deleteSpecifiedBondLines(bondLineIdxToBeDeleted);
  }

  private deleteSpecifiedBondLines(bondLineIdxToBeDeleted: number[]): void {
    this.bondLines = this.bondLines.filter((_, idx) => !bondLineIdxToBeDeleted.includes(idx));
    this.bondedAtoms = this.bondedAtoms.filter((_, idx) => !bondLineIdxToBeDeleted.includes(idx));
  }

  private getBondedAtoms(): number[][] {
    const bonds = MolfileHandler.getInstance(this.molfileV2K).pairsOfBondedAtoms;
    return bonds.map((pair) => Array.from(pair));
  }

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

  getShiftedBonds(shift: number): number[][] {
    return this.bondedAtoms.map((bondedPair) => {
      return bondedPair.map((atomId) => atomId + shift);
    });
  }

  getAtomAndBondLines(newBondedAtoms?: number[][]): {atomLines: string[], bondLines: string[]} {
    const bondedAtoms = newBondedAtoms ? newBondedAtoms : this.bondedAtoms;
    this.bondLines = this.bondLines.map((line: string, idx) => {
      const newBond1 = bondedAtoms[idx][0];
      const newBond2 = bondedAtoms[idx][1];
      return `${newBond1.toString().padStart(3, ' ')}${newBond2.toString().padStart(3, ' ')}${line.substring(6)}`;
    });
    return {atomLines: this.atomLines, bondLines: this.bondLines};
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    const bondIdx = this.bondedAtoms.findIndex((bondedPair) => {
      return bondedPair.includes(-rGroupId);
    });
    const attachmentAtom = this.bondedAtoms[bondIdx].find((atomId) => atomId !== -rGroupId)!;
    this.bondedAtoms[bondIdx] = [externalAtom, attachmentAtom];
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

  private getRGroupAtomicIdxById(rgroupId: number): number {
    return this.rgpIdToAtomicIdxMap[rgroupId];
  }

  getAttachmentAtomByRGroupId(rgroupId: number): number {
    let result = -1;
    for (const bondedPair of this.bondedAtoms) {
      if (bondedPair.includes(-rgroupId)) {
        result = bondedPair.find((idx) => idx !== -rgroupId)!;
        break;
      }
    }
    if (result === -1)
      throw new Error(`Cannot find attachment atom for RGP ${rgroupId} for monomer ${this.monomerSymbol}`);
    return result;
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

  removeRGroupsFromAtomBlock(rGroupIds: number[]): void {
    const rGroupIdsToAtomicIdx = new Map<number, number>();
    rGroupIds.forEach((rgroupId) => {
      const atomicIdx = this.getRGroupAtomicIdxById(rgroupId);
      if (atomicIdx !== undefined)
        rGroupIdsToAtomicIdx.set(rgroupId, atomicIdx);
      else
        throw new Error(`Cannot find atomic index for R group ${rgroupId} for monomer ${this.monomerSymbol}`);
    });
    const rGroupAtomicIndices = Array.from(rGroupIdsToAtomicIdx.values());
    this.atomLines = this.atomLines.filter((_, idx) => {
      return !rGroupAtomicIndices.includes(idx);
    });
    rGroupIds.forEach((rGroupId) => {
      let isAdded = false;
      const rGroupAtomicIdx = rGroupIdsToAtomicIdx.get(rGroupId)!;
      this.bondedAtoms = this.bondedAtoms.map((bondedPair, idx) => {
        return bondedPair.map((atomId) => {
          if (atomId - 1 > rGroupAtomicIdx)
            return atomId - 1;
          if (atomId - 1 === rGroupAtomicIdx) {
            if (isAdded)
              throw new Error(`R group ${rGroupId} is already handled`);
            // WARNING: atomId is 1-based
            isAdded = true;
            return -rGroupId;
          }
          return atomId;
        });
      });
      for (const [key, value] of rGroupIdsToAtomicIdx.entries()) {
        if (value > rGroupAtomicIdx)
          rGroupIdsToAtomicIdx.set(key, value - 1);
      };
    });
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
      result.push([{monomerIdx: backboneIdx, rGroupId: 2}, {monomerIdx: nextBackboneIdx, rGroupId: 1}]);
    }
    for (let i = 0; i < branchMonomerIndices.length; i++) {
      const branchIdx = branchMonomerIndices[i];
      const backboneIdx = branchIdx - 1;
      result.push([{monomerIdx: backboneIdx, rGroupId: 3}, {monomerIdx: branchIdx, rGroupId: 1}]);
    }
    return result;
  }
}

class ConnectionList {
  constructor(connectionList: string) {
    const splitted = connectionList.split(HELM_ITEM_SEPARATOR);
    splitted.forEach((connectionItem: string) => this.validateConnectionItem(connectionItem));
    this.connectionItems = splitted;
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
        // WARNING: monomer idx starts from 0
        const monomerIdx = parseInt(data[0]) - 1;
        const rGroupId = parseInt(data[1].slice(1));
        const bondData = {monomerIdx, rGroupId};
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
  }

  readonly bondData: Bond[][];

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
        bondPart.monomerIdx += shift;
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
          bond.monomerIdx += shift;
          data.push(bond);
        });
        result.push(data);
      });
    }
    return result;
  }
}

class Polymer {
  constructor(helm: string) {
    this.helm = new Helm(helm);
    this.polymerType = this.helm.polymerType;

    this.rGroupsToRemove = new Map<number, number[]>();
    this.helm.bondData.forEach((bond) => {
      bond.forEach((bondPart) => {
        const monomerIdx = bondPart.monomerIdx;
        const rGroupId = bondPart.rGroupId;
        if (!this.rGroupsToRemove.get(monomerIdx))
          this.rGroupsToRemove.set(monomerIdx, []);
        this.rGroupsToRemove.get(monomerIdx)!.push(rGroupId);
      });
    });
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private atomNumberShift = 0;
  private helm: Helm;
  private polymerType: HELM_POLYMER_TYPE;
  private rGroupsToRemove: Map<number, number[]>;

  addShiftedMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const monomerWrapper = new MonomerWrapper(monomerSymbol, this.atomNumberShift, this.polymerType);
    monomerWrapper.shiftCoordinates(shift);

    if (this.rGroupsToRemove.has(monomerIdx))
      monomerWrapper.deleteRGroups(this.rGroupsToRemove.get(monomerIdx)!);

    this.monomerWrappers.push(monomerWrapper);

    this.atomNumberShift += monomerWrapper.atomLines.length;
  }

  restoreBondsBetweenMonomers(): void {
    this.helm.bondData.forEach((bond) => {
      const bondPart1 = bond[0];
      const bondPart2 = bond[1];
      const monomerIdx1 = bondPart1.monomerIdx;
      const monomerIdx2 = bondPart2.monomerIdx;
      const rGroupId1 = bondPart1.rGroupId;
      const rGroupId2 = bondPart2.rGroupId;
      const monomerWrapper1 = this.monomerWrappers[monomerIdx1];
      const monomerWrapper2 = this.monomerWrappers[monomerIdx2];
      // const attachmentAtom1 = monomerWrapper1.getAttachmentAtomByRGroupId(rGroupId1);
      const attachmentAtom2 = monomerWrapper2.getAttachmentAtomByRGroupId(rGroupId2);
      monomerWrapper1.replaceRGroupWithAttachmentAtom(rGroupId1, attachmentAtom2);
      monomerWrapper2.deleteBondLineWithSpecifiedRGroup(rGroupId2);
    });
  }

  compileToMolfile(): string {
    const molfileHeader = '\nDatagrok\n';
    const atomLines: string[] = [];
    const bondLines: string[] = [];
    this.restoreBondsBetweenMonomers();
    console.log('restoring for:', this.helm.toString());
    this.monomerWrappers.forEach((monomerWrapper: MonomerWrapper) => {
      console.log(`atom lines`, monomerWrapper.atomLines);
      atomLines.push(...monomerWrapper.atomLines);
      console.log(`bond lines`, monomerWrapper.bondLines);
      bondLines.push(...monomerWrapper.bondLines);
    });
    const atomCount = atomLines.length;
    if (atomCount > MAX_ATOM_COUNT) {
      throw new Error(
        `Atom count in polymer ${this.helm.toString()} is ${atomCount} and exceeds ${MAX_ATOM_COUNT}`
      );
    }
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
