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

type RGroupBondPosition = {
  bondLineIdx: number,
  rgpIdx: number,
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
      polymer.addMonomer(monomerSymbol, monomerIdx, shift);
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
    monomerSymbol: string,
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
    return this.molfileWrapper.getAtomLines();
  }

  get bondLines(): string[] {
    return this.molfileWrapper.getBondLines(this.atomNumberShift);
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

class RGroupHandler {
  constructor(rGroupLines: string[], private atoms: MolfileAtoms, private bonds: MolfileBonds) {
    this.rGroupToAtomicIdxMap = this.getRGroupIdToAtomicIdxMap(rGroupLines);
  }

  /** Relates R group Id (starting from 1) to its atomic index within the
   * molfile  */
  rGroupToAtomicIdxMap: Map<number, number>;

  /** Maps R group id (starting from 1) to its position in the bond block  */
  private rGroupBondPositionMap = new Map<number, RGroupBondPosition>();

  getAtomicIdx(rGroupId: number): number | null {
    const atomicIdx = this.rGroupToAtomicIdxMap.get(rGroupId);
    return atomicIdx === undefined ? null : atomicIdx;
  }

  removeRGroups(rGroupIds: number[]): void {
    rGroupIds.forEach((rgroupId) => {
      const atomicIdx = this.rGroupToAtomicIdxMap.get(rgroupId);
      if (atomicIdx === undefined)
        throw new Error(`Cannot find atomic index for R group ${rgroupId}`);
    });

    const rGroupAtomicIndices = Array.from(this.rGroupToAtomicIdxMap.values());
    rGroupAtomicIndices.forEach((idx) => this.atoms.deleteAtom(idx));

    rGroupIds.forEach((rGroupId) => {
      const rGroupAtomicIdx = this.rGroupToAtomicIdxMap.get(rGroupId)!;

      if (this.rGroupBondPositionMap.has(rGroupId))
        throw new Error(`R group ${rGroupId} is already handled`);

      const rGroupPosition = this.bonds.replaceAtomByDummy(rGroupAtomicIdx + 1);
      if (rGroupPosition.bondLineIdx === -1)
        throw new Error(`Cannot delete atom for R group ${rGroupId}`);
      this.rGroupBondPositionMap.set(rGroupId, rGroupPosition);

      this.rGroupToAtomicIdxMap.delete(rGroupId);
      for (const [rGroupId, atomIdx] of this.rGroupToAtomicIdxMap) {
        if (atomIdx > rGroupAtomicIdx)
          this.rGroupToAtomicIdxMap.set(rGroupId, atomIdx - 1);
      };
    });
  }

  private getRGroupIdToAtomicIdxMap(lines: string[]): Map<number, number> {
    const rgroupLines = lines.filter((line: string) => line.startsWith('M  RGP'));

    const map = new Map<number, number>();
    rgroupLines.forEach((line: string) => {
      const indices = line.split(/\s+/).filter((item) => item)
        .slice(3).map((item) => parseInt(item));
      const atomIdxToRgpIdxList = new Array<[number, number]>(indices.length / 2);
      for (let i = 0; i < indices.length; i += 2)
        atomIdxToRgpIdxList[i / 2] = [indices[i + 1], indices[i] - 1];
      for (const [key, value] of atomIdxToRgpIdxList) {
        if (map.has(key))
          throw new Error(`R group ${key} is already in the map`);
        map.set(key, value);
      }
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
      if (map.has(rgpId))
        throw new Error(`R group ${rgpId} is already in the map`);
      map.set(rgpId, rgpAtomIdx);
    });

    const rGroupAtomicIndices = this.atoms.getRGroupAtomicIndices();
    const unaccounted = rGroupAtomicIndices.filter((idx) => !Array.from(map.values()).includes(idx));
    if (unaccounted.length !== 0)
      throw new Error(`Unaccounted R group indices: ${unaccounted}`);
    return map;
  }

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    const position = this.rGroupBondPositionMap.get(rGroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    const {bondLineIdx} = position;
    this.bonds.deleteBondLines([bondLineIdx]);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    const position = this.rGroupBondPositionMap.get(rGroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    const {bondLineIdx, rgpIdx} = position;
    this.bonds.bondedAtoms[bondLineIdx][rgpIdx] = externalAtom;
  }

  /** Atom id is molfile id starting from 1  */
  getAttachmentAtomIdByRGroupId(rgroupId: number): number {
    const position = this.rGroupBondPositionMap.get(rgroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rgroupId}`);
    const {bondLineIdx, rgpIdx} = position;
    return this.bonds.bondedAtoms[bondLineIdx][(rgpIdx + 1) % 2];
  }
}

class MolfileBonds {
  constructor(bondLines: string[]) {
    this.rawBondLines = bondLines;
    this.bondedPairs = this.rawBondLines.map((line: string) => {
      const firstAtom = parseInt(line.substring(0, 3));
      const secondAtom = parseInt(line.substring(3, 6));
      return [firstAtom, secondAtom];
    });
  }

  private bondedPairs: number[][] = [];
  private rawBondLines: string[] = [];

  /** Get bond lines with new values for bonded atoms  */
  getBondLines(atomNumberShift: number): string[] {
    return this.bondedPairs.map((bondedPair, idx) => {
      const newBond1 = bondedPair[0] + atomNumberShift;
      const newBond2 = bondedPair[1] + atomNumberShift;
      return `${ newBond1.toString().padStart(3, ' ')}${
        newBond2.toString().padStart(3, ' ')}${
        this.rawBondLines[idx].substring(6)}`;
    });
  }

  get bondedAtoms(): number[][] {
    return this.bondedPairs;
  }

  deleteBondLines(indices: number[]): void {
    this.rawBondLines = this.rawBondLines.filter((_, idx) => !indices.includes(idx));
    this.bondedPairs = this.bondedPairs.filter((_, idx) => !indices.includes(idx));
  }

  /** Replaces atom id (starts from 1) with -1 in the bond block  */
  replaceAtomByDummy(id: number): RGroupBondPosition {
    const replacedAtomPosition = {bondLineIdx: -1, rgpIdx: -1};
    this.bondedPairs = this.bondedPairs.map((bondedPair, bondLineIdx) => {
      return bondedPair.map((atomId, idx) => {
        if (atomId > id)
          return atomId - 1;
        if (atomId === id) {
          replacedAtomPosition.bondLineIdx = bondLineIdx;
          replacedAtomPosition.rgpIdx = idx;
          return -1;
        }
        return atomId;
      });
    });
    return replacedAtomPosition;
  }
}

class MolfileAtoms {
  constructor(atomLines: string[]) {
    this.rawAtomLines = atomLines;
    this.coordinates = this.rawAtomLines.map((line: string) => {
      const x = parseFloat(line.substring(0, 10));
      const y = parseFloat(line.substring(10, 20));
      return {x, y};
    });
  }

  private coordinates: {x: number, y: number}[] = [];
  private rawAtomLines: string[] = [];

  get atomCoordinates(): {x: number, y: number}[] {
    return this.coordinates;
  }

  get atomLines(): string[] {
    return this.rawAtomLines.map((line: string, idx: number) => {
      const coordinates = this.coordinates[idx];
      const x = coordinates.x.toFixed(4).padStart(10, ' ');
      const y = coordinates.y.toFixed(4).padStart(10, ' ');
      return `${x}${y}${line.substring(20)}`;
    });
  }

  deleteAtom(index: number): void {
    this.coordinates.splice(index, 1);
    this.rawAtomLines.splice(index, 1);
  }

  shift(shift: {x: number, y: number}): void {
    this.coordinates = this.coordinates.map((coordinates) => {
      const newX = coordinates.x + shift.x;
      const newY = coordinates.y + shift.y;
      if (isNaN(newX) || isNaN(newY))
        throw new Error(`Cannot shift coordinates by ${shift.x}, ${shift.y}`);
      return {x: newX, y: newY};
    });
  }

  rotate(angle: number): void {
    this.coordinates = this.coordinates.map((coordinates) => {
      const x = coordinates.x;
      const y = coordinates.y;
      const newX = x * Math.cos(angle) - y * Math.sin(angle);
      const newY = x * Math.sin(angle) + y * Math.cos(angle);
      if (isNaN(newX) || isNaN(newY))
        throw new Error(`Cannot rotate coordinates by ${angle}`);
      return {x: newX, y: newY};
    });
  }

  getRGroupAtomicIndices(): number[] {
    return this.rawAtomLines.map((line: string, idx: number) => {
      if (line.includes('R#'))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
  }
}

class MolfileWrapper {
  constructor(molfileV2K: string, private monomerSymbol: string) {
    const lines = molfileV2K.split('\n');

    const atomCountIdx = {begin: 0, end: 3};
    const bondCountIdx = {begin: 3, end: 6};
    const countsLineIdx = 3;
    const atomBlockIdx = 4;

    const atomCount = parseInt(lines[countsLineIdx].substring(atomCountIdx.begin, atomCountIdx.end));
    const bondCount = parseInt(lines[countsLineIdx].substring(bondCountIdx.begin, bondCountIdx.end));

    const atomLines = lines.slice(atomBlockIdx, atomBlockIdx + atomCount);
    this.atoms = new MolfileAtoms(atomLines);

    const bondLines = lines.slice(atomBlockIdx + atomCount, atomBlockIdx + atomCount + bondCount);
    this.bonds = new MolfileBonds(bondLines);

    this.rGroups = new RGroupHandler(lines, this.atoms, this.bonds);

    this.shiftMonomerToDefaultPosition();
  }

  private atoms: MolfileAtoms;
  private bonds: MolfileBonds;
  private rGroups: RGroupHandler;

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    this.rGroups.deleteBondLineWithSpecifiedRGroup(rGroupId);
  }

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.atoms.shift(shift);
  }

  rotateCoordinates(angle: number): void {
    this.atoms.rotate(angle);
  }

  getBondLines(atomNumberShift: number): string[] {
    return this.bonds.getBondLines(atomNumberShift);
  }

  getAtomLines(): string[] {
    return this.atoms.atomLines;
  }

  removeRGroupsFromAtomBlock(rGroupIds: number[]): void {
    this.rGroups.removeRGroups(rGroupIds);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    this.rGroups.replaceRGroupWithAttachmentAtom(rGroupId, externalAtom);
  }

  getAttachmentAtomByRGroupId(rgroupId: number): number {
    return this.rGroups.getAttachmentAtomIdByRGroupId(rgroupId);
  }

  private shiftR1GroupToOrigin(): void {
    const r1Idx = this.rGroups.getAtomicIdx(1);
    if (r1Idx === null)
      throw new Error(`Cannot find R1 group for monomer ${this.monomerSymbol}`);
    const {x, y} = this.atoms.atomCoordinates[r1Idx];
    this.atoms.shift({x: -x, y: -y});
  }

  private alignR2AlongX(): void {
    const r2Idx = this.rGroups.getAtomicIdx(2);
    if (r2Idx === null)
      throw new Error(`Cannot find R2 group for monomer ${this.monomerSymbol}`);
    const r2Coordinates = this.atoms.atomCoordinates[r2Idx];
    const tan = r2Coordinates.y / r2Coordinates.x;
    const angle = Math.atan(tan);
    if (isNaN(angle))
      throw new Error(`Cannot calculate angle for R2 group for monomer ${this.monomerSymbol}`);
    this.rotateCoordinates(-angle);
  }

  private shiftMonomerToDefaultPosition(): void {
    this.shiftR1GroupToOrigin();
    const r2Idx = this.rGroups.getAtomicIdx(2);
    if (r2Idx !== null)
      this.alignR2AlongX();
  }
}

/** Wrapper over simple polymer substring of HELM, like RNA1{d(A)p}  */
class SimplePolymer {
  constructor(private simplePolymer: string) {
    this.polymerType = this.getPolymerType();
    this.idx = this.getIdx();
    const {monomers, monomerTypes} = this.getMonomerSymbolsAndTypes();
    this.monomers = monomers;
    this.monomerTypes = monomerTypes;
  }

  readonly polymerType: string;
  readonly monomers: string[];
  private idx: number;
  private monomerTypes: HELM_MONOMER_TYPE[];

  /** Simple polymer id in the form 'polymer type' + 'index'  */
  get id(): string {
    return this.polymerType + this.idx.toString();
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

  private getIdx(): number {
    const regex = new RegExp(`${this.polymerType}([0-9]+){`);
    const match = this.simplePolymer.match(regex);
    if (!match)
      throw new Error(`Cannot parse simple polymer id from ${this.simplePolymer}`);
    const id = parseInt(match[1]);
    return id;
  }

  private getMonomerSymbolsAndTypes(): {monomers: string[], monomerTypes: HELM_MONOMER_TYPE[]} {
    const helmWrapperRegex = new RegExp(`${this.polymerType}${this.idx}{|}`, 'g');
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

  /** Get list of pairs for bonded monomers, monomers indexed locally
   * (within the simple polymer)  */
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

  /** List of pairs for bonded monomers, monomers indexed globally (withing the
   * complex polymer) */
  readonly bondData: Bond[][];

  private simplePolymers: SimplePolymer[];
  private connectionList?: ConnectionList;

  toString() {
    return this.helm;
  }

  getPolymerTypeByMonomerIdx(monomerGlobalIdx: number): HELM_POLYMER_TYPE {
    const simplePolymer = this.getSimplePolymerByMonomerIdx(monomerGlobalIdx);
    const polymerType = simplePolymer.polymerType;
    return polymerType as HELM_POLYMER_TYPE;
  }

  private getSimplePolymerByMonomerIdx(monomerGlobalIdx: number): SimplePolymer {
    if (this.simplePolymers.length === 1)
      return this.simplePolymers[0];
    const shifts = this.getMonomerIdxShifts();
    const shiftValues = Object.values(shifts);
    const lowerBound = shiftValues.sort((a, b) => a - b).find(
      (shift) => monomerGlobalIdx >= shift
    );
    if (lowerBound === undefined)
      throw new Error(`Cannot find simple polymer for monomer ${monomerGlobalIdx}`);
    const simplePolymerId = Object.keys(shifts).find((simplePolymerId) => shifts[simplePolymerId] === lowerBound)!;
    const simplePolymer = this.simplePolymers.find((simplePolymer) => simplePolymer.id === simplePolymerId)!;
    return simplePolymer;
  }

  private shiftBondMonomerIds(shift: number, bonds: Bond[][]): void {
    bonds.forEach((bond) => {
      bond.forEach((bondPart) => {
        bondPart.monomerIdx += shift;
      });
    });
  }

  private getMonomerIdxShifts(): {[simplePolymerId: string]: number} {
    const result: {[simplePolymerId: string]: number} = {};
    let shift = 0;
    this.simplePolymers.forEach((simplePolymer) => {
      result[simplePolymer.id] = shift;
      shift += simplePolymer.monomers.length;
    });
    return result;
  }

  private getBondData(): Bond[][] {
    const shifts = this.getMonomerIdxShifts();
    const result: Bond[][] = [];
    this.simplePolymers.forEach((simplePolymer) => {
      const bondData = simplePolymer.getBondData();
      const shift = shifts[simplePolymer.id];
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

    this.monomerToRGroupsMap = new Map<number, number[]>();
    this.helm.bondData.forEach((bond) => {
      bond.forEach((bondPart) => {
        const monomerIdx = bondPart.monomerIdx;
        const rGroupId = bondPart.rGroupId;
        if (!this.monomerToRGroupsMap.get(monomerIdx))
          this.monomerToRGroupsMap.set(monomerIdx, []);
        this.monomerToRGroupsMap.get(monomerIdx)!.push(rGroupId);
      });
    });
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private atomCountShift = 0;
  private helm: Helm;
  /** Maps global monomer index to r-group ids (starting from 1) participating
   * in connection */
  private monomerToRGroupsMap: Map<number, number[]>;

  addMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const polymerType = this.helm.getPolymerTypeByMonomerIdx(monomerIdx);
    const monomerWrapper = new MonomerWrapper(monomerSymbol, this.atomCountShift, polymerType);
    monomerWrapper.shiftCoordinates(shift);

    // if (this.monomerToRGroupsMap.has(monomerIdx))
    //   monomerWrapper.deleteRGroups(this.monomerToRGroupsMap.get(monomerIdx)!);

    this.monomerWrappers.push(monomerWrapper);
    this.atomCountShift += monomerWrapper.atomLines.length;
  }

  private restoreBondsBetweenMonomers(): void {
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

    // this.restoreBondsBetweenMonomers();

    this.monomerWrappers.forEach((monomerWrapper: MonomerWrapper) => {
      atomLines.push(...monomerWrapper.atomLines);
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
