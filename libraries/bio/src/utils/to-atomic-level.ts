/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '../../src/viewers/web-logo';
import {HELM_CORE_FIELDS, HELM_CORE_LIB_MONOMER_SYMBOL} from './monomer-utils';

const V2K_RGP_SHIFT = 8;
const V2K_RGP_LINE = 'M  RGP';

const V3K_ATOMS_DATA_SHIFT = 4;
const V3K_BOND_DATA_SHIFT = 4;
const V3K_COUNTS_SHIFT = 7;
const V3K_BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB\n';
const V3K_END_CTAB_BLOCK = 'M  V30 END CTAB\n';
const V3K_BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM\n';
const V3K_END_ATOM_BLOCK = 'M  V30 END ATOM\n';
const V3K_BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND\n';
const V3K_END_BOND_BLOCK = 'M  V30 END BOND\n';
const V3K_BEGIN_DATA_LINE = 'M  V30 ';
const V3K_END = 'M  END\n';

type AtomData = {
  atomIdx: number[],
  atomType: string[],
  x: number[], // Cartesian coordiantes
  y: number[],
}

type BondData = {
  // The node to which the radical is attached from the "left" (the order is specified by Macromolecule sequence):
  bondCount: number,
  leftNode: number,
  rightNode: number,
  leftRemovedNode: number, // "leftmost" r-group node of the monomer
  rightRemovedNode: number,
  leftRemovedBond: number, // idx of the bond between leftNode and leftRemovedNode
  rightRemovedBond: number,
}

type MonomerData = {
  molfileV3K: string,
  atomData: AtomData,
  bondData: BondData,
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
    reconstructed[row] = convertMacromolToMolfileV3K(monomerSeq, monomersDict);
    console.log(reconstructed[row]);
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
    const molfileV3K = convertMolfileToV3K(removeRGroupLine(molfileV2K), moduleRdkit);
    // todo: field name to constant
    const bondData = parseBondData(molfileV2K, molfileV3K);
    const atomData = parseAtomData(molfileV3K);
    adjustBackbone(atomData, bondData);
    // molfileV3K = getFromattedMolfileV3K(molfileV3K, atomData, bondData);
    return {molfileV3K: molfileV3K, atomData: atomData, bondData: bondData};
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
    // todo: verify that there are no multi-element rgroups, or consider how to
    // transform them
    rgroup = rgroup.replace(/(\[|\]|\*|:|\d)/g, '');
    rgroupsArray.push(rgroup);
  }
  return rgroupsArray;
}

function substituteRGroups(molfileV2K: string, rGroups: string[]) {
  let modifiedMolfile = molfileV2K;
  for (let value of rGroups) {
    // todo: handle hydrogens
    if (value = 'H')
      value = 'Li';
    modifiedMolfile = modifiedMolfile.replace('R#', value);
  }
  return modifiedMolfile;
}

// necessary to build a correct molfilev3k after substitution of r-groups in v2k
function removeRGroupLine(molfileV2K: string): string {
  // todo: consider the case of multiple rgp lines in a molfile
  const begin = molfileV2K.indexOf(V2K_RGP_LINE);
  const end = molfileV2K.indexOf('\n', begin) + 1;
  return molfileV2K.substring(0, begin) + molfileV2K.substring(end);
}

// todo: type of moduleRdkit
// todo: consider the use of unified converter (relies on creation of moduleRdkit
// on each iteration, though)
function convertMolfileToV3K(molfileV2K: string, moduleRdkit: any) {
  const molObj = moduleRdkit.get_mol(molfileV2K);
  const molfileV3K = molObj.get_v3Kmolblock();
  molObj.delete();
  return molfileV3K;
}

function parseBondData(molfileV2K: string, molfileV3K: string): BondData {
  // todo: consider the case when there is no simple leftmost/rightmost choice
  // todo: consider the case when there are multiple consequent M  RGP lines,
  // like in HELMCoreLibrary nucleotides

  const removedNodes = getRemovedNodes(molfileV2K);

  const counts = parseAtomAndBondCounts(molfileV3K); // todo: rename get to parse

  const indices = parseLRNodesAndBonds(molfileV3K, counts, removedNodes);

  return {
    bondCount: counts.bondCount,
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
  const begin = molfileV2K.indexOf(V2K_RGP_LINE, 0) + V2K_RGP_SHIFT;
  const end = molfileV2K.indexOf('\n', begin);

  const rgpStringParsed = molfileV2K.substring(begin, end).replaceAll('  ', ' ').replaceAll('  ', ' ').split(' ');
  const rgpIndices = rgpStringParsed.map((el) => parseInt(el));
  rgpIndices.shift(); // rgpIndices[0] is the number of R-groups
  return rgpIndices;
}

function parseLRNodesAndBonds(
  v3KMolblock: string,
  counts: {atomCount: number, bondCount: number},
  removedNodes: {left: number, right: number}
): {leftNode: number, rightNode: number, leftRemovedBond: number, rightRemovedBond: number} {
  let leftNode = 0;
  let rightNode = 0;
  let leftRemovedBond = 0;
  let rightRemovedBond = 0;
  // const counts = parseAtomAndBondCounts(v3KMolblock);

  let idxV3KBondBlock = v3KMolblock.indexOf(V3K_BEGIN_BOND_BLOCK);
  idxV3KBondBlock = v3KMolblock.indexOf('\n', idxV3KBondBlock);
  let begin = idxV3KBondBlock;
  let end = idxV3KBondBlock;

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
} // todo unify v3KMolblock and molfileV3K

function parseAtomAndBondCounts(v3KMolblock: string): {atomCount: number, bondCount: number} {
  v3KMolblock = v3KMolblock.replaceAll('\r', ''); // to handle old and new sdf standards

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

// function getFromattedMolfileV3K(
function adjustBackbone(
  // molfileV3K: string,
  atoms: AtomData,
  bondData: BondData
): void {
  const atomCount = atoms.atomIdx.length;
  const leftNode = bondData.leftNode;
  const rightNode = bondData.rightNode;

  // center tha backbone at origin
  const xCenter = (atoms.x[leftNode] + atoms.x[rightNode])/2;
  const yCenter = (atoms.y[leftNode] + atoms.y[rightNode])/2;
  for (let i = 0; i < atomCount; i++) {
    atoms.x[i] -= xCenter;
    atoms.y[i] -= yCenter;
  }

  rotateCenteredBackbone(atoms, leftNode, rightNode);

  // shift the backbone so that leftNode is at origin
  const xShift = atoms.x[rightNode];
  for (let i = 0; i < atomCount; ++i)
    atoms.x[i] += xShift;

  // return updateMolfileCoordinates(molfileV3K, atoms);
}

function parseAtomData(molfileV3K: string): AtomData {
  const counts = parseAtomAndBondCounts(molfileV3K);
  const atomIdx = new Array(counts.atomCount);
  const atomType = new Array(counts.atomCount);
  const x = new Array(counts.atomCount);
  const y = new Array(counts.atomCount);

  let begin = molfileV3K.indexOf(V3K_BEGIN_ATOM_BLOCK); // V3000 atoms block
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;

  for (let i = 0; i < counts.atomCount; i++) {
    begin = molfileV3K.indexOf('V30', begin) + V3K_ATOMS_DATA_SHIFT;
    end = molfileV3K.indexOf(' ', begin);
    atomIdx[i] = parseInt(molfileV3K.substring(begin, end));

    begin = end + 1; // todo: remove code copying
    end = molfileV3K.indexOf(' ', begin);
    atomType[i] = molfileV3K.substring(begin, end);

    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin);
    x[i] = parseFloat(molfileV3K.substring(begin, end));

    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin);
    y[i] = parseFloat(molfileV3K.substring(begin, end));

    begin = molfileV3K.indexOf('\n', begin) + 1;
  }

  return {atomIdx: atomIdx, atomType: atomType, x: x, y: y};
}

// Rotate the centered backbone so that "leftNode" is on the left and "rightNode" is on the right
function rotateCenteredBackbone(atoms: AtomData, leftNode: number, rightNode: number): void {
  let angle = 0;

  if (atoms.x[leftNode] === 0) { // both vertices are on OY
    angle = atoms.y[leftNode] > atoms.y[rightNode] ? Math.PI/2 : -Math.PI/2;
  } else if (atoms.y[leftNode] === 0) { // both vertices are on OX
    angle = atoms.x[leftNode] > atoms.x[rightNode] ? Math.PI : 0;
  } else {
    const tangent = atoms.y[leftNode]/atoms.x[leftNode];
    if (atoms.x[leftNode] < atoms.x[rightNode])
      angle = tangent > 0 ? -Math.atan(tangent) : Math.atan(tangent);
    else
      angle = tangent > 0 ? Math.PI - Math.atan(tangent) : Math.atan(tangent) - Math.PI;
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < atoms.x.length; ++i) {
    const xAdd = atoms.x[i];
    atoms.x[i] = xAdd*cos - atoms.y[i]*sin;
    atoms.y[i] = xAdd*sin + atoms.y[i]*cos;
  }
}

function convertMacromolToMolfileV3K(
  monomerSeq: string[],
  monomersDict: Map<string, MonomerData>
): string {
  if (monomerSeq.length === 0) {
    // todo: exception
    return '';
  } else if (monomerSeq.length === 1) {
    // todo: should it return rotated/aligned graph?
    return monomersDict.get(monomerSeq[0])!.molfileV3K; // todo: !
  } else {
    return concatenateMolfiles(monomerSeq, monomersDict);
  }
}

// there are two or more monomers in a sequence
function concatenateMolfiles(
  monomerSeq: string[],
  monomersDict: Map<string, MonomerData>
): string {
  let reconstructedAtomBlock = '';
  let reconstructedBondBlock = '';
  // todo: collection block?
  let totalNodeCount = 0;
  let totalBondCount = 0;
  let xShift = 0;

  // todo: handle first/last monomers separately

  for (let i = 0; i < monomerSeq.length; ++i) {
    const monomer = monomerSeq[i];
    const monomerData = monomersDict.get(monomer)!; // todo: exceptions
    const atomIdx = monomerData.atomData.atomIdx;
    const x = monomerData.atomData.x;
    const y = monomerData.atomData.y;
    const leftNode = monomerData.bondData.leftNode;
    const rightNode = monomerData.bondData.leftNode;
    const leftRemovedNode = monomerData.bondData.leftRemovedNode;
    const rightRemovedNode = monomerData.bondData.rightRemovedNode;
    const leftRemovedBond = monomerData.bondData.leftRemovedBond;
    const rightRemovedBond = monomerData.bondData.rightRemovedBond;
    const atomCount = atomIdx.length;
    const bondCount = monomerData.bondData.bondCount;
    const molfile = monomerData.molfileV3K;
    const totalShift = xShift - x[leftNode - 1]; // todo: BUG

    for (let j = 0; j < atomCount; ++j) {
      const node = atomIdx[j]; // todo: do we need this atomIdx at all? or should add bondIdx to BondData
      let updatedNode;
      if (node != leftRemovedNode || node != rightRemovedNode) {
        updatedNode = getIdxAfterLRItemsRemoval(node, leftRemovedNode, rightRemovedNode);
        updatedNode += totalNodeCount; // todo: hydrogen

        const updatedX = Math.round(10000*(x[j] + totalShift))/10000;
        const updatedY = Math.round(10000*(y[j]))/10000;

        const idxAtomBlock = molfile.indexOf(V3K_BEGIN_ATOM_BLOCK);
        const begin = molfile.indexOf(V3K_BEGIN_DATA_LINE + node.toString(), idxAtomBlock);
        const end = molfile.indexOf('\n', begin) + 1; // including \n
        let atomLine = molfile.slice(begin, end);
        atomLine = updateAtomLine(atomLine, updatedNode, updatedX, updatedY);

        reconstructedAtomBlock += atomLine;
      }
    }

    for (let j = 0; j < bondCount; ++j) {
      const bond = j + 1;
      let updatedBond;
      if (bond != leftRemovedBond || bond != rightRemovedBond) {
        updatedBond = getIdxAfterLRItemsRemoval(bond, leftRemovedBond, rightRemovedBond);
        updatedBond += totalBondCount;

        const idxBondBlock = molfile.indexOf(V3K_BEGIN_BOND_BLOCK);
        const begin = molfile.indexOf(V3K_BEGIN_DATA_LINE + bond.toString(), idxBondBlock);
        const end = molfile.indexOf('\n', begin) + 1; // including \n
        let bondLine = molfile.slice(begin, end);
        bondLine = updateBondLine(bondLine, updatedBond, leftRemovedNode, rightRemovedNode, totalNodeCount);

        reconstructedBondBlock += bondLine;
      }
    }

    totalNodeCount += atomCount - 2;
    totalBondCount += bondCount - 2;
    xShift += x[rightNode] - x[leftNode] + 1; // todo: ?


    // todo: refactor the following part
    if (i === monomerSeq.length - 1) {
      totalNodeCount++;
      const shift = xShift + 0.2;
      // todo: improve this part
      reconstructedAtomBlock += 'M  V30 ' + totalNodeCount + ' O ' + shift + ' 0 0.000000 0\n';
    }
    totalBondCount++;
    if (i === monomerSeq.length - 1) {
      const rightTerminal = (rightNode > leftRemovedNode && rightNode > rightRemovedNode) ?
        rightNode + totalNodeCount - (atomCount - 2) - 3:
        (rightNode > leftRemovedNode ||
          rightNode > rightRemovedNode) ?
          rightNode + totalNodeCount - (atomCount - 2) - 2 :
          rightNode + totalNodeCount - (atomCount - 2) - 1;
      reconstructedBondBlock += 'M  V30 ' + totalBondCount + ' 1 ' + rightTerminal + ' ' + totalNodeCount + '\n';
    } else {
      const rightTerminal = (rightNode > leftRemovedNode && rightNode > rightRemovedNode) ?
        rightNode + totalNodeCount - (atomCount - 2) - 2:
        (rightNode > leftRemovedNode ||
          rightNode > rightRemovedNode) ?
          rightNode + totalNodeCount - (atomCount - 2) - 1 :
          rightNode + totalNodeCount - (atomCount - 2);

      const nextMonomer = monomerSeq[i + 1];
      const nextMonomerData = monomersDict.get(nextMonomer)!; // todo: exceptions
      const nextLeftNode = nextMonomerData.bondData.leftNode;
      // const nextRightNode = nextMonomerData.bondData.leftNode;
      const nextLeftRemovedNode = monomerData.bondData.leftRemovedNode;
      const nextRightRemovedNode = monomerData.bondData.rightRemovedNode;

      const leftTerminal = (nextLeftNode > nextLeftRemovedNode &&
        nextLeftNode > nextRightRemovedNode) ?
        nextLeftNode + totalNodeCount - 2 :
        (nextLeftNode > nextLeftRemovedNode ||
          nextLeftNode > nextRightRemovedNode) ?
          nextLeftNode + totalNodeCount - 1 :
          nextLeftNode + totalNodeCount;

      reconstructedBondBlock += 'M  V30 ' + bondCount + ' 1 ' + rightTerminal + ' ' + leftTerminal + '\n';
    }
  }

  const reconstructedCountsLine = 'M  V30 COUNTS ' + totalNodeCount + ' ' + totalBondCount + ' 0 0 0\n';

  // produce molfile
  let reconstructedMolfile = '\nDatagrok macromolecule handler\n\n';
  reconstructedMolfile += '  0  0  0  0  0  0            999 V3000\n';
  reconstructedMolfile += V3K_BEGIN_CTAB_BLOCK;
  reconstructedMolfile += reconstructedCountsLine;
  reconstructedMolfile += V3K_BEGIN_ATOM_BLOCK;
  reconstructedMolfile += reconstructedAtomBlock;
  reconstructedMolfile += V3K_END_ATOM_BLOCK;
  reconstructedMolfile += V3K_BEGIN_BOND_BLOCK;
  reconstructedMolfile += reconstructedBondBlock;
  reconstructedMolfile += V3K_END_BOND_BLOCK;
  reconstructedMolfile += V3K_END_CTAB_BLOCK;
  reconstructedMolfile += V3K_END;

  return reconstructedMolfile;
}

function getIdxAfterLRItemsRemoval(idx: number, leftRemoved: number, rightRemoved: number) {
  if (idx > leftRemoved || idx > rightRemoved) {
    if (idx > leftRemoved && idx > rightRemoved)
      return idx - 2;
    else
      return idx - 1;
  } else {
    return idx;
  }
}

function updateAtomLine(
  atomLine: string,
  updatedNode: number,
  updatedX: number,
  updatedY: number
): string {
  // position of atom index in the line
  const beginNode = atomLine.indexOf('V30') + V3K_ATOMS_DATA_SHIFT;
  const endNode = atomLine.indexOf(' ', beginNode);

  // beginning of X coordinate in the line
  let beginX = endNode + 1;
  beginX = atomLine.indexOf(' ', beginX) + 1;

  // end of Y coordinate in the line
  let endY = atomLine.indexOf(' ', beginX) + 1;
  endY = atomLine.indexOf(' ', endY);

  const result = atomLine.slice(0, beginNode) + updatedNode +
    atomLine.slice(endNode, beginX) + updatedX + ' ' + updatedY +
    atomLine.slice(endY);

  return result;
}

function updateBondLine(
  bondLine: string,
  updatedBond: number,
  leftRemovedNode: number,
  rightRemovedNode: number,
  totalNodeCount: number
) {
  // position of bond index in the line
  const beginBond = bondLine.indexOf(V3K_BEGIN_DATA_LINE) + V3K_ATOMS_DATA_SHIFT;
  const endBond = bondLine.indexOf(' ', beginBond);

  // todo: remove the Indian code
  // begginning of the first atom of the bond
  let beginFirst = endBond + 1;
  beginFirst = bondLine.indexOf(' ', beginFirst) + 1;

  // end of the last atom of the bond
  let endSecond = bondLine.indexOf(' ', beginFirst) + 1;
  const tmp = bondLine.indexOf(' ', endSecond);
  // -1 in case there is no ' ' left
  endSecond = (tmp === -1) ? bondLine.indexOf('\n', endSecond) : tmp;

  // rewrite the value of the bond index
  let begin = bondLine.indexOf(V3K_BEGIN_DATA_LINE) + V3K_ATOMS_DATA_SHIFT;
  let end = bondLine.indexOf(' ', begin);
  let updatedLine = bondLine.slice(0, begin) + updatedBond + bondLine.slice(end);

  // rewrite the indices for the two atoms of the bond
  for (let k = 0; k < 2; ++k) {
    begin = end + 1; // todo: code copying?
    begin = bondLine.indexOf(' ', begin) + 1;
    // todo: exclude -1
    // we need min in case of the last atom idx, after which \n goes
    end = Math.min(bondLine.indexOf('\n', begin), bondLine.indexOf(' ', begin));
    let atomIdx = parseInt(bondLine.substring(begin, end));
    atomIdx = getIdxAfterLRItemsRemoval(atomIdx, leftRemovedNode, rightRemovedNode);
    atomIdx += totalNodeCount;

    updatedLine = bondLine.slice(0, begin) + atomIdx + bondLine.slice(end);
  }
  return updatedLine;
}
