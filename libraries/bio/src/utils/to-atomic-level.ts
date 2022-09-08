/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '../../src/viewers/web-logo';
import {HELM_CORE_FIELDS, HELM_CORE_LIB_MONOMER_SYMBOL} from './monomer-utils';

type BondIndices = {
  // The node to which the radical is attached from the "left" (the order is specified by Macromolecule sequence):
  leftNode: number,
  rightNode: number,
  leftRemovedNode: number, // "leftmost" r-group node of the monomer
  rightRemovedNode: number,
  leftRemovedBond: number, // idx of the bond between leftNode and leftRemovedNode
  rightRemovedBond: number,
}

type MonomerData = {
  molfileV3K: string, // molfile corresponding to a rotated monomer
  bondIndices: BondIndices,
}

type AtomData = {
  atomIdx: number,
  atomType: string,
  x: number,
  y: number,
}

export async function _toAtomicLevel(
  df: DG.DataFrame,
  macroMolCol: DG.Column<string>,
  monomersLibList: any[]
): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }

  // jagged array of monomer symbols
  const monomerSequencesArray: string[][] = getMonomerSequencesArray(macroMolCol);
  const monomersDict = await getMonomersDict(monomerSequencesArray, monomersLibList);
  const reconstructed: string[] = new Array(macroMolCol.length);
  for (let row = 0; row < monomerSequencesArray.length; ++row) {
    const monomerSeq = monomerSequencesArray[row];
    reconstructed[row] = bind(monomerSeq, monomersDict);
  }
  const newCol = DG.Column.fromStrings('reconstructed', reconstructed);
  // todo: name properly
  newCol.semType = DG.SEMTYPE.MOLECULE;
  newCol.tags[DG.TAGS.UNITS] = 'molblock';
  df.columns.add(newCol, true);
  await grok.data.detectSemanticTypes(df);
}

function getFormattedMonomerLib(monomersLibList: any[]) {
  const map = new Map<string, any>();
  monomersLibList.forEach(
    (it) => {
      // todo: generalize for the case of nucleotides
      if (it['polymerType'] === 'PEPTIDE') {
        // const monomerMap = new Map<string, any>();
        const monomerObject: { [key: string]: any } = {};
        HELM_CORE_FIELDS.forEach((field) => {
          monomerObject[field] = it[field];
        });
        map.set(it[HELM_CORE_LIB_MONOMER_SYMBOL], monomerObject);
      }
    });
  return map;
}

function getMonomerSequencesArray(macroMolCol: DG.Column<string>): string[][] {
  const result: string[][] = new Array(macroMolCol.length);

  const colUnits = macroMolCol.tags[DG.TAGS.UNITS];
  const separator = macroMolCol.getTag('separator');
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(colUnits, separator);
  for (let row = 0; row < macroMolCol.length; ++row) {
    const macroMolecule = macroMolCol.get(row);
    // todo: handle the exception case when macroMolecule can be null
    result[row] = macroMolecule ? splitterFunc(macroMolecule) : [];
  }
  return result;
}

async function getMonomersDict(
  monomerSequencesArray: string[][],
  monomersLibList: any[]
): Promise<Map<string, MonomerData>> {
  const formattedMonomerLib = getFormattedMonomerLib(monomersLibList);
  const monomersDict = new Map<string, MonomerData>();

  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');

  for (let row = 0; row < monomerSequencesArray.length; ++row) {
    const monomerSeq: string[] = monomerSequencesArray[row];
    for (const sym of monomerSeq) {
      if (!monomersDict.has(sym)) {
        const monomerData: MonomerData | null = getMonomerData(sym, formattedMonomerLib, moduleRdkit);
        if (monomerData)
          monomersDict.set(sym, monomerData);
        // todo: handle exception when there is no monomer with symbol sym in
        // monomersLibList
      }
    }
  }
  return monomersDict;
}

function getMonomerData(
  monomerSymbol: string,
  formattedMonomerLib: Map<string, any>,
  moduleRdkit: any // todo: specify type
): MonomerData | null {
  if (!formattedMonomerLib.has(monomerSymbol)) {
    return null;
  } else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    // todo: create a constant for HELM molfile field
    // todo: field names to constants
    const rgroups = parseRGroups(libObject['rgroups']);
    const molfileV2K = substituteRGroups(libObject['molfile'], rgroups);
    let molfileV3K = convertToV3K(molfileV2K, moduleRdkit);
    // todo: field name to constant
    const bondIndices = parseBondIndices(molfileV2K, molfileV3K);
    molfileV3K = getFromattedMolfileV3K(molfileV3K, bondIndices);
    return {molfileV3K: molfileV3K, bondIndices: bondIndices};
  }
}

function parseRGroups(rgroupObjList: any[]): string[] {
  // specifically for HELMCoreLibrary
  // considered only monoatomic rgroups
  // supposing that elements in rgroupObjList are sorted w.r.t. the rgroups idx
  // todo: possible generalizations
  const rgroupsArray: string[] = [];
  for (const obj of rgroupObjList) {
    let rgroup: string = obj['capGroupSmiles'];
    // todo: verify that there are no multi-element rgroups
    rgroup = rgroup.replace(/(\[|\]|\*|:|\d)/g, '');
    rgroupsArray.push(rgroup);
  }
  return rgroupsArray;
}

function substituteRGroups(molfileV2K: string, rGroups: string[]) {
  let modifiedMolfile = molfileV2K;
  for (const value of rGroups)
    modifiedMolfile = modifiedMolfile.replace('R#', value);
  return modifiedMolfile;
}

// todo: type of moduleRdkit
// todo: consider the use of unified converter (relies creation of moduleRdkit
// on each iteration, though)
function convertToV3K(molfileV2K: string, moduleRdkit: any) {
  const molObj = moduleRdkit.get_mol(molfileV2K);
  const molfileV3K = molObj.get_v3Kmolblock();
  molObj.delete();
  return molfileV3K;
}

function parseBondIndices(molfileV2K: string, molfileV3K: string): BondIndices {
  // todo: consider the case when there is no simple leftmost/rightmost choice
  // todo: consider the case when there are multiple consequent M  RGP lines,
  // like in HELMCoreLibrary nucleotides

  const removedNodes = getRemovedNodes(molfileV2K);

  const indices = getLRNodesAndBonds(molfileV3K, removedNodes);

  return {
    leftNode: indices.leftNode,
    rightNode: indices.rightNode,
    leftRemovedNode: removedNodes.left,
    rightRemovedNode: removedNodes.right,
    leftRemovedBond: indices.leftRemovedBond,
    rightRemovedBond: indices.rightRemovedBond,
  };
}

function getRemovedNodes(molfileV2K: string): {left: number, right: number} {
  const rGroupIndices = parseRGroupIndices(molfileV2K);

  const leftRemovedNode = rGroupIndices[0];
  const rightRemovedNode = rGroupIndices[rGroupIndices.length - 2];
  return {left: leftRemovedNode, right: rightRemovedNode};
}

function parseRGroupIndices(molfileV2K: string): number[] {
  // todo: handle the exceptional case when there is not enough rgroups
  // todo: handle the exceptional case when the order is different
  const V2K_RGP_SHIFT = 8;

  const begin = molfileV2K.indexOf('M  RGP', 0) + V2K_RGP_SHIFT;
  const end = molfileV2K.indexOf('\n', begin);

  const rgpStringParsed = molfileV2K.substring(begin, end).replaceAll('  ', ' ').replaceAll('  ', ' ').split(' ');
  const rgpIndices = rgpStringParsed.map((el) => parseInt(el));
  rgpIndices.shift(); // rgpIndices[0] is the number of R-groups
  return rgpIndices;
}

function getLRNodesAndBonds(
  v3KMolblock: string,
  removedNodes: {left: number, right: number}
): {leftNode: number, rightNode: number, leftRemovedBond: number, rightRemovedBond: number} {
  let leftNode = 0;
  let rightNode = 0;
  let leftRemovedBond = 0;
  let rightRemovedBond = 0;
  const counts = getAtomAndBondCountsV3K(v3KMolblock);

  let idxV3KBondBlock = v3KMolblock.indexOf('M  V30 BEGIN BOND');
  idxV3KBondBlock = v3KMolblock.indexOf('\n', idxV3KBondBlock);
  let begin = idxV3KBondBlock;
  let end = idxV3KBondBlock;

  const V3K_BOND_DATA_SHIFT = 4;

  let i = 0;

  while ((i < counts.bondCount) && (leftNode === 0 || rightNode === 0)) {
    begin = v3KMolblock.indexOf('V30', begin) + V3K_BOND_DATA_SHIFT;
    end = v3KMolblock.indexOf('\n', begin);
    const bondStringParsed = v3KMolblock.substring(begin, end).replaceAll('  ', ' ').replaceAll('  ', ' ').split(' ');
    const bondData = bondStringParsed.map((el) => parseInt(el));

    if (bondData[2] === removedNodes.left) { // bondData[2] is the 1st node/atom of the bond
      leftNode = bondData[3]; // bondData[3] is the 2nd node/atom of the bond
      leftRemovedBond = bondData[0]; // bondData[0] is the idx of the associated bond/edge
    } else if (bondData[3] === removedNodes.left) {
      leftNode = bondData[2];
      leftRemovedBond = bondData[0];
    } else if (bondData[2] === removedNodes.right) {
      rightNode = bondData[3];
      rightRemovedBond = bondData[0];
    } else if (bondData[3] === removedNodes.right) {
      rightNode = bondData[2];
      rightRemovedBond = bondData[0];
    }
    ++i;
  }

  return {
    leftNode: leftNode,
    rightNode: rightNode,
    leftRemovedBond: leftRemovedBond,
    rightRemovedBond: rightRemovedBond
  };
}

function getAtomAndBondCountsV3K(v3KMolblock: string): {atomCount: number, bondCount: number} {
  v3KMolblock = v3KMolblock.replaceAll('\r', ''); // to handle old and new sdf standards

  const V3K_COUNTS_SHIFT = 7;

  // parse atom count
  let idxBegin = v3KMolblock.indexOf('COUNTS') + V3K_COUNTS_SHIFT;
  let idxEnd = v3KMolblock.indexOf(' ', idxBegin);
  const numOfAtoms = parseInt(v3KMolblock.substring(idxBegin, idxEnd));

  // parse bond count
  idxBegin = idxEnd + 1;
  idxEnd = v3KMolblock.indexOf(' ', idxBegin);
  const numOfBonds = parseInt(v3KMolblock.substring(idxBegin, idxEnd));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

function getFromattedMolfileV3K(molfileV3K: string, bondIndices: BondIndices): string {
  const atom = parseAtomData(molfileV3K);

  // let formattedMolfileV3K = rotateV3KBackbone(molfileV3K, bondIndices);
}

function parseAtomData(molfileV3K: string): AtomData[] {
  const atomsArray = [];

  
  const numbers = extractAtomAndBondCountsV3K(v3KMolblock);
  let begin = v3KMolblock.indexOf('M  V30 BEGIN ATOM'); // V3000 block for atom coordinates
  begin = v3KMolblock.indexOf('\n', begin);
  let end = begin;

  const atomIndices: number[] = Array(numbers.atomCount);
  const atomTypes: string[] = Array(numbers.atomCount);
  const x: number[] = Array(numbers.atomCount);
  const y: number[] = Array(numbers.atomCount);

  for (let i = 0; i < numbers.atomCount; i++) {
    begin = v3KMolblock.indexOf('V30', begin) + 4;
    end = v3KMolblock.indexOf(' ', begin);
    atomIndices[i] = parseInt(v3KMolblock.substring(begin, end));

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    atomTypes[i] = v3KMolblock.substring(begin, end);

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    x[i] = parseFloat(v3KMolblock.substring(begin, end));

    begin = end + 1;
    end = v3KMolblock.indexOf(' ', begin);
    y[i] = parseFloat(v3KMolblock.substring(begin, end));

    begin = v3KMolblock.indexOf('\n', begin) + 1;
  }

  return {atomIndices: atomIndices, atomTypes: atomTypes, x: x, y: y};


}
