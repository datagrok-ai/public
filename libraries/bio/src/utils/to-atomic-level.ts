/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '../../src/viewers/web-logo';
import {HELM_CORE_FIELDS, HELM_CORE_LIB_MONOMER_SYMBOL} from './monomer-utils';

const V2K_RGP_SHIFT = 8;
const V2K_RGP_LINE = 'M  RGP';

const V3K_COUNTS_SHIFT = 7;
const V3K_IDX_SHIFT = 7;
const V3K_BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB\n';
const V3K_END_CTAB_BLOCK = 'M  V30 END CTAB\n';
const V3K_BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM\n';
const V3K_END_ATOM_BLOCK = 'M  V30 END ATOM\n';
const V3K_BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND\n';
const V3K_END_BOND_BLOCK = 'M  V30 END BOND\n';
const V3K_BEGIN_DATA_LINE = 'M  V30 ';
const V3K_END = 'M  END\n';

// todo: docstrings

type AtomData = {
  atomType: string[],
  x: number[], // Cartesian coordiantes
  y: number[],
  leftNode: number, // node bound to leftRemovedNode, indexing starts from 1
  rightNode: number,
  leftRemovedNode: number, // "leftmost" r-group node of the monomer
  rightRemovedNode: number,
  kwargs: string[], // MDLV30 atom line may contain keyword args
}

type BondData = {
  // The node to which the radical is attached from the "left" (the order is specified by Macromolecule sequence):
  bondType: number[],
  atomPair: number[][],
  kwargs: string[], // MDLV30 atom line may contain keyword args
}

type MolGraph = {
  // molfileV3K: string,
  atoms: AtomData,
  bonds: BondData,
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
    const molGraph = concatenateMolGraphs(monomerSeq, monomersDict);
    reconstructed[row] = convertMolGraphToMolfileV3K(molGraph);
    console.log('Reconstructed' + reconstructed[row]);
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
): Promise<Map<string, MolGraph>> {
  const formattedMonomerLib = getFormattedMonomerLib(monomersLibList);
  const monomersDict = new Map<string, MolGraph>();

  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');

  for (let row = 0; row < monomerSequencesArray.length; ++row) {
    const monomerSeq: string[] = monomerSequencesArray[row];
    for (const sym of monomerSeq) {
      if (!monomersDict.has(sym)) {
        const monomerData: MolGraph | null = getMolGraph(sym, formattedMonomerLib, moduleRdkit);
        if (monomerData)
          monomersDict.set(sym, monomerData);
        // todo: handle exception when there is no monomer with symbol sym in
        // monomersLibList
      }
    }
  }
  return monomersDict;
}

function getMolGraph(
  monomerSymbol: string,
  formattedMonomerLib: Map<string, any>,
  moduleRdkit: any // todo: specify type
): MolGraph | null {
  if (!formattedMonomerLib.has(monomerSymbol)) {
    return null;
  } else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    // todo: create constants for HELM molfile fields
    const rgroups = parseRGroupTypes(libObject['rgroups']);
    const molfileV2K = substituteRGroups(libObject['molfile'], rgroups);
    // removeRGroupLine needed to reconstruct the correct V3K
    const molfileV3K = convertMolfileToV3K(removeRGroupLine(molfileV2K), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);
    const bondData = parseBondData(molfileV3K, counts.bondCount);
    const removedNodes = getRemovedNodes(molfileV2K);
    const atomData = parseAtomData(removedNodes, bondData, molfileV3K, counts.atomCount);
    removeIntermediateHydrogen(atomData, bondData);
    adjustGraph(atomData);
    return {atoms: atomData, bonds: bondData};
  }
}

function parseRGroupTypes(rgroupObjList: any[]): string[] {
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
  const rGroupIndices = parseRGroupIndices(molfileV2K);
  const idx = new Array((rGroupIndices.length - 1)/2);
  for (let i = 0; i < rGroups.length; ++i)
    // idx[i] is responsible for the correct order of substituted rgroups
    idx[i] = rGroupIndices[2*(i + 1)] - 1;
  for (let i = 0; i < rGroups.length - 1; ++i) {
    // an extra check of correctness of rgroups order
    if (rGroupIndices[2 * (i + 1) + 1] < rGroupIndices[2 * (i + 1) + 1])
      throw new Error('Wrong order of rgroups');
  }
  for (let i = 0; i < rGroups.length; ++i) {
    const value = rGroups[idx[i]];
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

function parseBondData(molfileV3K: string, bondCount: number): BondData {
  // todo: consider the case when there is no simple leftmost/rightmost choice
  // todo: consider the case when there are multiple consequent M  RGP lines,
  // like in HELMCoreLibrary nucleotides

  const bondType = new Array(bondCount);
  const atomPair = new Array(bondCount);
  const kwargs = new Array(bondCount);

  let begin = molfileV3K.indexOf(V3K_BEGIN_BOND_BLOCK);
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;
  for (let i = 0; i < bondCount; ++i) {
    const parsedValues = new Array(3);
    begin = molfileV3K.indexOf(V3K_BEGIN_DATA_LINE, end) + V3K_IDX_SHIFT;
    end = molfileV3K.indexOf(' ', begin);
    for (let k = 0; k < 3; ++k) {
      begin = end + 1;
      end = Math.min(molfileV3K.indexOf('\n', begin), molfileV3K.indexOf(' ', begin));
      parsedValues[k] = parseInt(molfileV3K.slice(begin, end));
    }
    bondType[i] = parsedValues[0];
    atomPair[i] = parsedValues.slice(1);
    const endOfLine = molfileV3K.indexOf('\n', begin) + 1;
    kwargs[i] = molfileV3K.slice(end, endOfLine);
  }

  return {
    bondType: bondType,
    atomPair: atomPair,
    kwargs: kwargs,
  };
}

function getRemovedNodes(molfileV2K: string): {left: number, right: number} {
  // todo: verify that in other monomer libraries the order of removed nodes in
  // RGP line is the same as in HELMCoreLibrary
  const rGroupIndices = parseRGroupIndices(molfileV2K);

  const leftRemovedNode = rGroupIndices[1]; // rgpIndices[0] is the number of R-groups
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
  return rgpIndices;
}

function parseRetainedLRNodes(
  atomPair: number[][],
  removedNodes: {left: number, right: number}
): {left: number, right: number} {
  let leftNode = 0;
  let rightNode = 0;

  let i = 0;
  while ((i < atomPair.length) && (leftNode === 0 || rightNode === 0)) {
    if (atomPair[i][0] === removedNodes.left)
      leftNode = atomPair[i][1];
    else if (atomPair[i][1] === removedNodes.left)
      leftNode = atomPair[i][0];
    else if (atomPair[i][0] === removedNodes.right)
      rightNode = atomPair[i][1];
    else if (atomPair[i][1] === removedNodes.right)
      rightNode = atomPair[i][0];
    ++i;
  }

  return {
    left: leftNode,
    right: rightNode,
  };
} // todo unify v3KMolblock and molfileV3K

function parseAtomAndBondCounts(v3KMolblock: string): {atomCount: number, bondCount: number} {
  // todo: rewrite using COUNTS const
  v3KMolblock = v3KMolblock.replaceAll('\r', ''); // to handle old and new sdf standards

  // parse atom count
  let begin = v3KMolblock.indexOf('COUNTS') + V3K_COUNTS_SHIFT;
  let end = v3KMolblock.indexOf(' ', begin);
  const numOfAtoms = parseInt(v3KMolblock.substring(begin, end));

  // parse bond count
  begin = end + 1;
  end = v3KMolblock.indexOf(' ', begin);
  const numOfBonds = parseInt(v3KMolblock.substring(begin, end));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

function parseAtomData(
  removedNodes: {left: number, right: number},
  bondData: BondData,
  molfileV3K: string,
  atomCount: number
): AtomData {
  const atomType = new Array(atomCount);
  const x = new Array(atomCount);
  const y = new Array(atomCount);
  const kwargs = new Array(atomCount);

  const atomPair = bondData.atomPair;
  const retainedNodes = parseRetainedLRNodes(atomPair, removedNodes);

  let begin = molfileV3K.indexOf(V3K_BEGIN_ATOM_BLOCK); // V3000 atoms block
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;

  for (let i = 0; i < atomCount; i++) {
    begin = molfileV3K.indexOf(V3K_BEGIN_DATA_LINE, begin) + V3K_IDX_SHIFT;
    end = molfileV3K.indexOf(' ', begin); // skip the idx row

    begin = end + 1; // proceed to atom type row
    end = molfileV3K.indexOf(' ', begin);
    atomType[i] = molfileV3K.substring(begin, end);

    begin = end + 1; // todo: remove code copying
    end = molfileV3K.indexOf(' ', begin);
    x[i] = parseFloat(molfileV3K.substring(begin, end));

    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin);
    y[i] = parseFloat(molfileV3K.substring(begin, end));

    begin = end;
    end = molfileV3K.indexOf('\n', begin) + 1;
    kwargs[i] = molfileV3K.slice(begin, end);

    begin = end;
  }

  return {
    atomType: atomType,
    x: x,
    y: y,
    leftNode: retainedNodes.left,
    rightNode: retainedNodes.right,
    leftRemovedNode: removedNodes.left,
    rightRemovedNode: removedNodes.right,
    kwargs: kwargs,
  };
}

// remove hydrogen nodes that are not L/R removed nodes
function removeIntermediateHydrogen(atomData: AtomData, bondData: BondData): void {
  let i = 0;
  while (i < atomData.atomType.length) {
    if ( atomData.atomType[i] == 'H' &&
      i != atomData.leftRemovedNode && // todo: update lr nodes!
      i != atomData.rightRemovedNode) {
      removeNodeAndBonds(atomData, bondData, i);
      --i;
    }
    ++i;
  }
}

// this must return a copy of the molGraph
function removeTerminalHydrogen(molGraph: MolGraph, removedNode: number): MolGraph {
  const molGraphCopy: MolGraph = structuredClone(molGraph);
  const removedIdx = removedNode - 1;
  if (molGraphCopy.atoms.atomType[removedIdx] === 'H')
    removeNodeAndBonds(molGraphCopy.atoms, molGraphCopy.bonds, removedIdx);
  return molGraphCopy;
}

// remove node 'idx' and the associated bonds
function removeNodeAndBonds(atomData: AtomData, bondData: BondData, removedIdx: number) {
  // update nodes
  const removedAtom = removedIdx + 1; // atom indices start from 1
  atomData.atomType.splice(removedIdx, 1);
  atomData.x.splice(removedIdx, 1);
  atomData.y.splice(removedIdx, 1);
  atomData.kwargs.splice(removedIdx, 1);

  atomData.leftNode = (atomData.leftNode > removedAtom) ? --atomData.leftNode : atomData.leftNode;
  atomData.rightNode = (atomData.rightNode > removedAtom) ? --atomData.rightNode : atomData.rightNode;
  atomData.leftRemovedNode = (atomData.leftRemovedNode > removedAtom) ?
    --atomData.leftRemovedNode : atomData.leftRemovedNode;
  atomData.rightRemovedNode = (atomData.rightRemovedNode > removedAtom) ?
    --atomData.rightRemovedNode : atomData.rightRemovedNode;

  // update bonds
  let i = 0;
  while (i < bondData.atomPair.length) {
    const firstAtom = bondData.atomPair[i][0];
    const secondAtom = bondData.atomPair[i][1];
    if (firstAtom === removedAtom ||
      secondAtom === removedAtom) {
      bondData.atomPair.splice(i, 1);
      bondData.bondType.splice(i, 1);
      bondData.kwargs.splice(i, 1);
      --i;
    } else {
      bondData.atomPair[i][0] = (firstAtom > removedAtom) ? firstAtom - 1 : firstAtom;
      bondData.atomPair[i][1] = (secondAtom > removedAtom) ? secondAtom - 1 : secondAtom;
    }
    ++i;
  }
}

// function getFromattedMolfileV3K(
function adjustGraph(atoms: AtomData): void {
  const atomCount = atoms.x.length;
  const leftNodeIdx = atoms.leftNode - 1;
  const rightNodeIdx = atoms.rightNode - 1;

  // center tha backbone at origin
  const xCenter = (atoms.x[leftNodeIdx] + atoms.x[rightNodeIdx])/2;
  const yCenter = (atoms.y[leftNodeIdx] + atoms.y[rightNodeIdx])/2;
  for (let i = 0; i < atomCount; i++) {
    atoms.x[i] -= xCenter;
    atoms.y[i] -= yCenter;
  }

  rotateCenteredGraph(atoms);

  // shift the backbone so that leftNode is at origin
  const xShift = atoms.x[rightNodeIdx];
  for (let i = 0; i < atomCount; ++i)
    atoms.x[i] += xShift;
}

// Rotate the centered backbone so that "leftNode" is on the left and "rightNode" is on the right
function rotateCenteredGraph(atoms: AtomData): void {
  let angle = 0;
  const leftNodeIdx = atoms.leftNode - 1;
  const rightNodeIdx = atoms.rightNode - 1;
  const x = atoms.x;
  const y = atoms.y;

  // console.log('Start transform');
  // console.log('atoms:' + atoms.atomType);
  // console.log('x:' + x);
  // console.log('y:' + y);

  if (x[leftNodeIdx] === 0) { // both vertices are on OY
    angle = y[leftNodeIdx] > y[rightNodeIdx] ? Math.PI/2 : -Math.PI/2;
  } else if (y[leftNodeIdx] === 0) { // both vertices are on OX
    angle = x[leftNodeIdx] > x[rightNodeIdx] ? Math.PI : 0;
  } else {
    const tangent = y[leftNodeIdx]/x[leftNodeIdx];
    if (x[leftNodeIdx] < x[rightNodeIdx])
      angle = tangent > 0 ? -Math.atan(tangent) : Math.atan(tangent);
    else
      angle = tangent > 0 ? Math.PI - Math.atan(tangent) : Math.atan(tangent) - Math.PI;
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  for (let i = 0; i < x.length; ++i) {
    const xAdd = x[i];
    x[i] = Math.round(10000 * (xAdd*cos - y[i]*sin))/10000;
    y[i] = Math.round(10000 * (xAdd*sin + y[i]*cos))/10000;
    // x[i] = x[i]*cos - y[i]*sin;
    // y[i] = x[i]*sin + y[i]*cos;
  }
  // console.log('End transform');
  // console.log('x:' + x);
  // console.log('y:' + y);
}

function concatenateMolGraphs(
  monomerSeq: string[],
  monomersDict: Map<string, MolGraph>
): MolGraph {
  if (monomerSeq.length === 0)
    throw new Error('The monomerSeq is empty');
    // todo: handle/modify exception

  // remove "left" hydrogens if any
  let result = removeTerminalHydrogen(
    monomersDict.get(monomerSeq[0])!,
    monomersDict.get(monomerSeq[0])!.atoms.leftRemovedNode
  );

  concatenateMonomers(result, monomerSeq, monomersDict);

  result = removeTerminalHydrogen(result, result.atoms.rightRemovedNode);
  return result;
}

function concatenateMonomers(
  result: MolGraph,
  monomerSeq: string[],
  monomersDict: Map<string, MolGraph>
): void {
  for (let i = 1; i < monomerSeq.length; i++) {
    const xShift = result.atoms.x[result.atoms.rightRemovedNode - 1];
    removeNodeAndBonds(result.atoms, result.bonds, result.atoms.rightRemovedNode);

    const nextMonomer: MolGraph = structuredClone(monomersDict.get(monomerSeq[i])!);
    removeNodeAndBonds(nextMonomer.atoms, nextMonomer.bonds, nextMonomer.atoms.leftRemovedNode);

    const nodeIdxShift = result.atoms.atomType.length;
    shiftNodes(nextMonomer, nodeIdxShift);
    shiftXCoordinate(nextMonomer, xShift);

    // bind the two monomers with a peptide bond
    // todo: eliminate hardcoding of the peptide bond
    result.bonds.bondType.push(1);
    result.bonds.atomPair.push([result.atoms.rightNode, nextMonomer.atoms.leftNode]);
    result.bonds.kwargs.push(result.bonds.kwargs[result.atoms.rightRemovedNode - 1]);

    // update result with the values of nextMonomer
    // todo: consider introduction of consts for the sake of readibilty
    result.atoms.rightNode = nextMonomer.atoms.rightNode;
    result.atoms.rightRemovedNode = nextMonomer.atoms.rightRemovedNode;
    result.atoms.atomType = result.atoms.atomType.concat(nextMonomer.atoms.atomType);
    result.atoms.kwargs = result.atoms.kwargs.concat(nextMonomer.atoms.kwargs);

    result.bonds.kwargs = result.bonds.kwargs.concat(nextMonomer.bonds.kwargs);
    result.bonds.atomPair = result.bonds.atomPair.concat(nextMonomer.bonds.atomPair);
    result.bonds.bondType = result.bonds.bondType.concat(nextMonomer.bonds.bondType);
  }
}

function shiftNodes(monomer: MolGraph, nodeIdxShift: number) {
  monomer.atoms.leftNode += nodeIdxShift;
  monomer.atoms.rightNode += nodeIdxShift;
  monomer.atoms.rightRemovedNode += nodeIdxShift;

  const bondsCount = monomer.bonds.atomPair.length;
  for (let i = 0; i < bondsCount; ++i) {
    const atomPair = monomer.bonds.atomPair[i];
    for (let j = 0; j < 2; ++j)
      atomPair[j] += nodeIdxShift;
  }
}

function shiftXCoordinate(monomer: MolGraph, xShift: number) {
  const x = monomer.atoms.x;
  for (let i = 0; i < x.length; ++i)
    x[i] = Math.round(10000*(x[i] + xShift))/10000;
}

function convertMolGraphToMolfileV3K(molGraph: MolGraph): string {
  // counts line
  const atomType = molGraph.atoms.atomType;
  const x = molGraph.atoms.x;
  const y = molGraph.atoms.y;
  const atomKwargs = molGraph.atoms.kwargs;
  const bondType = molGraph.bonds.bondType;
  const atomPair = molGraph.bonds.atomPair;
  const bondKwargs = molGraph.bonds.kwargs;
  const atomCount = atomType.length;
  const bondCount = molGraph.bonds.bondType.length;

  const reconstructedCountsLine = 'M  V30 COUNTS ' + atomCount + ' ' + bondCount + ' 0 0 0\n';

  // atom block
  let reconstructedAtomBlock = '';
  for (let i = 0; i < atomCount; ++i) {
    const atomIdx = i + 1;
    // todo: uncomment after debugging undefined
    // const coordinate = [x[i].toString(), y[i].toString()];
    // for (let k = 0; k < 2; ++k) {
    //   const formatted = coordinate[k].toString().split('.');
    //   if (formatted.length === 1)
    //     formatted.push('0');
    //   formatted[1] = formatted[1].padEnd(6, '0');
    //   coordinate[k] = formatted.join('.');
    // }
    // const atomLine = V3K_BEGIN_DATA_LINE + atomIdx + ' ' + atomType[i] + ' ' +
    //   coordinate[0] + ' ' + coordinate[1] + ' ' + atomKwargs[i];

    const atomLine = V3K_BEGIN_DATA_LINE + atomIdx + ' ' + atomType[i] + ' ' +
      x[i] + ' ' + y[i] + ' ' + atomKwargs[i];
    reconstructedAtomBlock += atomLine;
  }

  // bond block
  let reconstructedBondBlock = '';
  for (let i = 0; i < bondCount; ++i) {
    const bondIdx = i + 1;
    const firstAtom = atomPair[i][0];
    const secondAtom = atomPair[i][1];
    const bondLine = V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
      firstAtom + ' ' + secondAtom + ' ' + bondKwargs[i];
    reconstructedBondBlock += bondLine;
  }

  // todo: collection block?

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
