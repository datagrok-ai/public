/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {HELM_POLYMER_TYPE, HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import {MonomerLibManager} from './monomer-lib/lib-manager';

import {_package} from '../package';

const enum V2K_CONST {
  MAX_ATOM_COUNT = 999,
  RGP_LINE_START = 'M  RGP',
  ATOM_ALIAS_LINE_START = 'A  ',
}

const HELM_SECTION_SEPARATOR = '$';
const HELM_ITEM_SEPARATOR = '|';
const R_GROUP_ELEMENT_SYMBOL = 'R#';
const HYDROGEN_SYMBOL = 'H';
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

/** Position of a node in the connection list / bond block  */
type PositionInBonds = {
  bondLineIdx: number,
  nodeIdx: number,
}

/** Translate HELM column into molfile column and append to the dataframe */
export async function helm2mol(df: DG.DataFrame, helmCol: DG.Column<string>): Promise<void> {
  const molCol = await getMolColumnFromHelm(df, helmCol);
  df.columns.add(molCol, true);
  await grok.data.detectSemanticTypes(df);
}


/** Translate HELM column into molfile column and append to the dataframe */
export async function getMolColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>
): Promise<DG.Column<string>> {
  const converter = new HelmToMolfileConverter(helmCol, df);
  const molCol = await converter.convertToRdKitBeautifiedMolfileColumn();
  molCol.semType = DG.SEMTYPE.MOLECULE;
  return molCol;
}

export async function getSmilesColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>
): Promise<DG.Column<string>> {
  const converter = new HelmToMolfileConverter(helmCol, df);
  const smilesCol = await converter.convertToSmiles();
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return smilesCol;
}

export class HelmToMolfileConverter {
  constructor(private helmColumn: DG.Column<string>, private df: DG.DataFrame) {
    this.helmColumn = helmColumn;
  }

  async convertToSmiles(): Promise<DG.Column<string>> {
    const smiles = await this.getSmilesList();
    const columnName = this.df.columns.getUnusedName(`smiles(${this.helmColumn.name})`);
    return DG.Column.fromStrings(columnName, smiles.map((molecule) => {
      if (molecule === null)
        return '';
      return molecule;
    }));
  }

  private async getSmilesList(): Promise<string[]> {
    const molfilesV2K = (await this.convertToMolfileV2KColumn()).toList();
    const smiles = molfilesV2K.map((mol) => DG.chem.convert(mol, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles));
    return smiles;
  }

  async convertToRdKitBeautifiedMolfileColumn(): Promise<DG.Column<string>> {
    const smiles = await this.getSmilesList();
    const rdKitModule: RDModule = await grok.functions.call('Chem:getRdKitModule');
    const beautifiedMols = smiles.map((item) =>{
      if (item === '')
        return null;
      const mol = rdKitModule.get_mol(item);
      if (!mol)
        return null;
      mol.normalize_depiction(1);
      mol.straighten_depiction(true);
      return mol;
    });
    const columnName = this.df.columns.getUnusedName(`molfile(${this.helmColumn.name})`);
    return DG.Column.fromStrings(columnName, beautifiedMols.map((mol) => {
      if (mol === null)
        return '';
      return mol.get_molblock();
    }));
  }

  async convertToMolfileV2KColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> = await this.getPolymerGraphColumn();
    const molfileList = polymerGraphColumn.toList().map(
      (pseudoMolfile: string, idx: number) => {
        const helm = this.helmColumn.get(idx);
        if (!helm)
          return '';
        let result = '';
        try {
          result = this.getPolymerMolfile(helm, pseudoMolfile);
        } catch (err: any) {
          const [errMsg, errStack] = errInfo(err);
          _package.logger.error(errMsg, undefined, errStack);
        } finally {
          return result;
        }
      });
    const molfileColName = this.df.columns.getUnusedName(`molfileV2K(${this.helmColumn.name})`);
    const molfileColumn = DG.Column.fromList('string', molfileColName, molfileList);
    return molfileColumn;
  }

  private async getPolymerGraphColumn(): Promise<DG.Column<string>> {
    const polymerGraphColumn: DG.Column<string> =
      await grok.functions.call('HELM:getMolfiles', {col: this.helmColumn});
    return polymerGraphColumn;
  }

  private getPolymerMolfile(helm: string, polymerGraph: string): string {
    const globalPositionHandler = new GlobalMonomerPositionHandler(polymerGraph);
    const polymer = new Polymer(helm);
    globalPositionHandler.monomerSymbols.forEach((monomerSymbol: string, monomerIdx: number) => {
      const shift = globalPositionHandler.getMonomerShifts(monomerIdx);
      polymer.addMonomer(monomerSymbol, monomerIdx, shift);
    });
    const polymerMolfile = polymer.compileToMolfile();
    return polymerMolfile;
  }
}

class GlobalMonomerPositionHandler {
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
    polymerType: HELM_POLYMER_TYPE,
  ) {
    const monomerLib = MonomerLibManager.instance.getBioLib();
    const monomer = monomerLib.getMonomer(polymerType, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.molfileWrapper = new MolfileWrapper(monomer.molfile, monomerSymbol);
    this.capGroupElements = monomer.rgroups.map((rgroup) => {
      const smiles = rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES] ||
        // WARNING: ignore because both key variants coexist in HELM Core Library!
        // @ts-ignore
        rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
      // extract the element symbol
      return smiles.replace(/(\[|\]|\*|:|\d)/g, '');
    });
  }

  private molfileWrapper: MolfileWrapper;
  private capGroupElements: string[] = [];

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.molfileWrapper.shiftCoordinates(shift);
  }

  getAtomLines(): string[] {
    return this.molfileWrapper.getAtomLines();
  }

  getBondLines(): string[] {
    return this.molfileWrapper.getBondLines();
  }

  removeBondedRGroups(rGroupIds: number[]): void {
    this.molfileWrapper.removeRGroups(rGroupIds);
  }

  capTrailingRGroups(): void {
    this.molfileWrapper.capRGroups(this.capGroupElements);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, attachmentAtomIdx: number): void {
    this.molfileWrapper.replaceRGroupWithAttachmentAtom(rGroupId, attachmentAtomIdx);
  };

  getAttachmentAtomByRGroupId(rGroupId: number): number {
    const attachmentAtom = this.molfileWrapper.getAttachmentAtomByRGroupId(rGroupId);
    return attachmentAtom;
  }

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    this.molfileWrapper.deleteBondLineWithSpecifiedRGroup(rGroupId);
  }

  shiftBonds(shift: number): void {
    this.molfileWrapper.shiftBonds(shift);
  }
}

class RGroupHandler {
  constructor(rGroupLines: string[], private atoms: MolfileAtoms, private bonds: MolfileBonds) {
    this.rGroupIdToAtomicIndexMap = this.getRGroupIdToAtomicIdxMap(rGroupLines);
  }

  /** Relates R group id (starting from 1) to its atomic index within the
   * molfile  */
  rGroupIdToAtomicIndexMap: Map<number, number>;

  /** Maps R group id (starting from 1) to its position in the bond block  */
  private rGroupBondPositionMap = new Map<number, PositionInBonds>();

  getAtomicIdx(rGroupId: number): number | null {
    const atomicIdx = this.rGroupIdToAtomicIndexMap.get(rGroupId);
    return atomicIdx === undefined ? null : atomicIdx;
  }

  private removeRGroupsFromAtomBlock(rGroupIds: number[]): void {
    rGroupIds.forEach((rgroupId) => {
      const atomicIdx = this.rGroupIdToAtomicIndexMap.get(rgroupId);
      if (atomicIdx === undefined)
        throw new Error(`Cannot find atomic index for R group ${rgroupId}`);
    });

    const rGroupAtomicIndices = Array.from(this.rGroupIdToAtomicIndexMap.entries()).filter(
      ([rGroupId, _]) => rGroupIds.includes(rGroupId)
    ).map(([_, atomicIdx]) => atomicIdx);
    this.atoms.deleteAtoms(rGroupAtomicIndices);
  }

  removeRGroups(rGroupIds: number[]): void {
    this.removeRGroupsFromAtomBlock(rGroupIds);

    rGroupIds.forEach((rGroupId) => {
      const dummyPosition = this.replaceRGroupInBondsByDummy(rGroupId);
      this.rGroupBondPositionMap.set(rGroupId, dummyPosition);
    });
  }

  /** Replace RGroups by -1, update associated maps, and return the position in
   * bond block */
  private replaceRGroupInBondsByDummy(rGroupId: number): PositionInBonds {
    const rGroupAtomicIdx = this.rGroupIdToAtomicIndexMap.get(rGroupId)!;

    if (this.rGroupBondPositionMap.has(rGroupId))
      throw new Error(`R group ${rGroupId} is already handled`);

    const positions = this.bonds.getPositionsInBonds(rGroupAtomicIdx + 1);
    if (positions.length === 0)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    if (positions.length > 1)
      throw new Error(`More than one position for R group ${rGroupId}`);

    const rGroupPosition = positions[0];

    this.bonds.replacePositionsInBondsByDummy([rGroupPosition]);
    this.bonds.removeAtomIdFromBonds(rGroupAtomicIdx + 1);
    this.removeRGroupFromAtomicIdxMap(rGroupId, rGroupAtomicIdx);

    return rGroupPosition;
  }

  private removeRGroupFromAtomicIdxMap(deletedId: number, deletedAtomicIdx: number): void {
    this.rGroupIdToAtomicIndexMap.delete(deletedId);
    for (const [rGroupId, rGroupAtomicIdx] of this.rGroupIdToAtomicIndexMap) {
      if (rGroupAtomicIdx > deletedAtomicIdx)
        this.rGroupIdToAtomicIndexMap.set(rGroupId, rGroupAtomicIdx - 1);
    }
  }

  private getRGroupIdToAtomicIdxMap(lines: string[]): Map<number, number> {
    function getAtomIdxToRgpIdxList(rgpLine: string): [number, number][] {
      const indices = rgpLine.split(/\s+/).filter((item) => item)
        .slice(3).map((item) => parseInt(item));
      const atomIdxToRgpIdxList = new Array<[number, number]>(indices.length / 2);
      for (let i = 0; i < indices.length; i += 2)
        atomIdxToRgpIdxList[i / 2] = [indices[i + 1], indices[i] - 1];
      return atomIdxToRgpIdxList;
    }

    const map = new Map<number, number>();

    const rgroupLines = lines.filter((line: string) => line.startsWith(V2K_CONST.RGP_LINE_START));
    rgroupLines.forEach((line: string) => {
      const atomIdxToRgpIdxList = getAtomIdxToRgpIdxList(line);
      for (const [key, value] of atomIdxToRgpIdxList) {
        if (map.has(key))
          throw new Error(`R group ${key} is already in the map`);
        map.set(key, value);
      }
    });

    const atomAliasLinesIndices = lines.map((line: string, idx: number) => {
      if (line.startsWith(V2K_CONST.ATOM_ALIAS_LINE_START))
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
    this.rGroupBondPositionMap.delete(rGroupId);
    this.rGroupIdToAtomicIndexMap.delete(rGroupId);
    // update values of other positions
    this.rGroupBondPositionMap.forEach((position) => {
      if (position.bondLineIdx > bondLineIdx)
        position.bondLineIdx -= 1;
    });
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    const position = this.rGroupBondPositionMap.get(rGroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    const {bondLineIdx, nodeIdx} = position;
    this.bonds.bondedAtoms[bondLineIdx][nodeIdx] = externalAtom;
  }

  /** Atom id is molfile id starting from 1  */
  getAttachmentAtomIdByRGroupId(rgroupId: number): number {
    const position = this.rGroupBondPositionMap.get(rgroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rgroupId}`);
    const {bondLineIdx, nodeIdx} = position;
    return this.bonds.bondedAtoms[bondLineIdx][(nodeIdx + 1) % 2];
  }

  /** WARNING: capping RGRoups and deletion of the bonded ones don't commute */
  capRGroups(capGroupElements: string[]): void {
    this.rGroupIdToAtomicIndexMap.forEach((atomicIdx, rGroupId) => {
      const element = capGroupElements[rGroupId - 1];
      if (element === HYDROGEN_SYMBOL) {
        this.removeRGroups([rGroupId]);
        this.deleteBondLineWithSpecifiedRGroup(rGroupId);
      } else {
        this.atoms.replaceElementSymbol(atomicIdx, element);
      }
    });
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
  getBondLines(): string[] {
    return this.bondedPairs.map((bondedPair, idx) => {
      if (bondedPair.some((atom) => atom === -1))
        throw new Error(`Bonded pair ${bondedPair} contains -1`);
      return `${bondedPair[0].toString().padStart(3, ' ')}${
        bondedPair[1].toString().padStart(3, ' ')
      }${this.rawBondLines[idx].substring(6)}`;
    });
  }

  get bondedAtoms(): number[][] {
    return this.bondedPairs;
  }

  deleteBondLines(indices: number[]): void {
    this.rawBondLines = this.rawBondLines.filter((_, idx) => !indices.includes(idx));
    this.bondedPairs = this.bondedPairs.filter((_, idx) => !indices.includes(idx));
  }

  /** Atom id starts from 1  */
  getPositionsInBonds(atomId: number): PositionInBonds[] {
    const positions: PositionInBonds[] = [];
    this.bondedPairs.forEach((bondedPair, bondLineIdx) => {
      bondedPair.forEach((atom, nodeIdx) => {
        if (atom === atomId)
          positions.push({bondLineIdx, nodeIdx});
      });
    });
    return positions;
  }

  replacePositionsInBondsByDummy(positions: PositionInBonds[], dummy?: number): void {
    if (dummy === undefined)
      dummy = -1;
    positions.forEach((position) => {
      const {bondLineIdx, nodeIdx} = position;
      this.bondedPairs[bondLineIdx][nodeIdx] = dummy!;
    });
  }

  removeAtomIdFromBonds(atomId: number): void {
    this.bondedPairs = this.bondedPairs.map((bondedPair) => {
      return bondedPair.map((id) => {
        if (id > atomId)
          return id - 1;
        return id;
      });
    });
  }

  shift(shift: number): void {
    this.bondedPairs = this.bondedPairs.map((bondedPair) => {
      return bondedPair.map((id) => id + shift);
    });
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

  replaceElementSymbol(atomIdx: number, newElementSymbol: string): void {
    this.rawAtomLines[atomIdx] = this.rawAtomLines[atomIdx].replace(R_GROUP_ELEMENT_SYMBOL, newElementSymbol);
  }

  deleteAtoms(indices: number[]): void {
    this.coordinates = this.coordinates.filter((_, idx) => !indices.includes(idx));
    this.rawAtomLines = this.rawAtomLines.filter((_, idx) => !indices.includes(idx));
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
      if (line.includes(R_GROUP_ELEMENT_SYMBOL))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
  }
}

class MolfileWrapper {
  constructor(molfileV2K: string, private monomerSymbol: string) {
    const lines = molfileV2K.split('\n');

    // TODO: port to consts
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

  getBondLines(): string[] {
    return this.bonds.getBondLines();
  }

  getAtomLines(): string[] {
    return this.atoms.atomLines;
  }

  removeRGroups(rGroupIds: number[]): void {
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

  shiftBonds(shift: number): void {
    this.bonds.shift(shift);
  }

  capRGroups(capGroupElements: string[]): void {
    this.rGroups.capRGroups(capGroupElements);
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
   * complex polymer scope) */
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

    this.bondedRGroupsMap = new Map<number, number[]>();
    this.helm.bondData.forEach((bond) => {
      bond.forEach((bondPart) => {
        const monomerIdx = bondPart.monomerIdx;
        const rGroupId = bondPart.rGroupId;
        if (!this.bondedRGroupsMap.get(monomerIdx))
          this.bondedRGroupsMap.set(monomerIdx, []);
        this.bondedRGroupsMap.get(monomerIdx)!.push(rGroupId);
      });
    });
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;
  /** Maps global monomer index to r-group ids (starting from 1) participating
   * in connection */
  private bondedRGroupsMap: Map<number, number[]>;

  addMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const polymerType = this.helm.getPolymerTypeByMonomerIdx(monomerIdx);
    const monomerWrapper = new MonomerWrapper(monomerSymbol, polymerType);
    monomerWrapper.shiftCoordinates(shift);

    this.monomerWrappers.push(monomerWrapper);
  }

  private removeRGroups(): void {
    this.monomerWrappers.forEach((monomerWrapper, monomerIdx) => {
      if (this.bondedRGroupsMap.has(monomerIdx))
        monomerWrapper.removeBondedRGroups(this.bondedRGroupsMap.get(monomerIdx)!);
      monomerWrapper.capTrailingRGroups();
    });
  }

  private getAtomNumberShifts(): number[] {
    const atomNumberShifts: number[] = [];
    let shift = 0;
    this.monomerWrappers.forEach((monomerWrapper) => {
      atomNumberShifts.push(shift);
      shift += monomerWrapper.getAtomLines().length;
    });
    return atomNumberShifts;
  }

  private restoreBondsBetweenMonomers(): void {
    this.helm.bondData.forEach((bond) => {
      const monomerIdx = bond.map((bondPart) => bondPart.monomerIdx);
      const rGroupId = bond.map((bondPart) => bondPart.rGroupId);
      const monomer = monomerIdx.map((idx) => this.monomerWrappers[idx]);

      const attachmentAtom = monomer[1].getAttachmentAtomByRGroupId(rGroupId[1]);
      monomer[0].replaceRGroupWithAttachmentAtom(rGroupId[0], attachmentAtom);
      monomer[1].deleteBondLineWithSpecifiedRGroup(rGroupId[1]);
    });
  }

  compileToMolfile(): string {
    const molfileHeader = '\nDatagrok\n';
    const atomLines: string[] = [];
    const bondLines: string[] = [];

    this.removeRGroups();

    const atomNumberShifts = this.getAtomNumberShifts();
    this.monomerWrappers.forEach((monomerWrapper, idx) => {
      monomerWrapper.shiftBonds(atomNumberShifts[idx]);
    });

    this.restoreBondsBetweenMonomers();

    this.monomerWrappers.forEach((monomerWrapper) => {
      atomLines.push(...monomerWrapper.getAtomLines());
      bondLines.push(...monomerWrapper.getBondLines());
    });

    const atomCount = atomLines.length;
    if (atomCount > V2K_CONST.MAX_ATOM_COUNT) {
      throw new Error(
        `Atom count in polymer ${this.helm.toString()} is ${atomCount} and exceeds ${V2K_CONST.MAX_ATOM_COUNT}`
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
