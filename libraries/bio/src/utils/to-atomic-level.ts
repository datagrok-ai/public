/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '../../src/viewers/web-logo';
import {HELM_CORE_FIELDS} from './monomer-utils';
import {RGROUP_FIELD, MOLFILE_FIELD, MONOMER_SYMBOL, CAP_GROUP_SMILES} from './const';

const V2K_RGP_SHIFT = 8;
const V2K_RGP_LINE = 'M  RGP';

const V3K_COUNTS_SHIFT = 14;
const V3K_IDX_SHIFT = 7;
const V3K_ATOM_COORDINATE_PRECISION = 6;
const V3K_HEADER_FIRST_LINE = '\nDatagrok macromolecule handler\n\n';
const V3K_HEADER_SECOND_LINE = '  0  0  0  0  0  0            999 V3000\n';
const V3K_BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB\n';
const V3K_END_CTAB_BLOCK = 'M  V30 END CTAB\n';
const V3K_BEGIN_COUNTS_LINE = 'M  V30 COUNTS ';
const V3K_COUNTS_LINE_ENDING = ' 0 0 0\n';
const V3K_BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM\n';
const V3K_END_ATOM_BLOCK = 'M  V30 END ATOM\n';
const V3K_BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND\n';
const V3K_END_BOND_BLOCK = 'M  V30 END BOND\n';
const V3K_BEGIN_DATA_LINE = 'M  V30 ';
const V3K_END = 'M  END\n';

// todo: test for separator/helm

type AtomData = {
  atomType: string[], // element symbol
  x: number[], // Cartesian coordiantes of the nodes
  y: number[],
  leftNode: number, // "leftmost" retained node (not removed upon chaining into polymer)
  rightNode: number,
  leftRemovedNode: number, // "leftmost" r-group node of the monomer (deleted upon chaining into polymer)
  rightRemovedNode: number,
  kwargs: string[], // MDLV30 atom line may contain keyword args
}

type BondData = {
  bondType: number[],
  atomPair: number[][], // indices of atoms, starting with 1
  kwargs: string[], // MDLV30 atom line may contain keyword args
}

/* Type storing the necessary information about atoms and bonds */
type MolGraph = {
  atoms: AtomData,
  bonds: BondData,
}

/* Convert Macromolecule column into Molecule column storing molfile V3000 with the help of a monomer library  */
export async function _toAtomicLevel(
  df: DG.DataFrame, macroMolCol: DG.Column<string>, monomersLibList: any[]
): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }

  // todo: remove?
  const pi = DG.TaskBarProgressIndicator.create('Restoring atomic structure...');

  const monomerSequencesArray: string[][] = getMonomerSequencesArray(macroMolCol);
  const monomersDict = await getMonomersDict(monomerSequencesArray, monomersLibList);
  const reconstructed: string[] = new Array(macroMolCol.length);
  const columnLength = monomerSequencesArray.length;
  for (let row = 0; row < columnLength; ++row) {
    const monomerSeq = monomerSequencesArray[row];
    const molGraph = monomerSeqToMolGraph(monomerSeq, monomersDict);
    reconstructed[row] = convertMolGraphToMolfileV3K(molGraph);
    pi.update(100 * (row + 1)/columnLength, '...complete');
  }

  // exclude name collisions
  const name = 'molfile(' + macroMolCol.name + ')';
  const newColName = df.columns.getUnusedName(name);
  const newCol = DG.Column.fromStrings(newColName, reconstructed);

  newCol.semType = DG.SEMTYPE.MOLECULE;
  newCol.tags[DG.TAGS.UNITS] = DG.UNITS.Molecule.MOLBLOCK;
  df.columns.add(newCol, true);
  await grok.data.detectSemanticTypes(df);

  pi.close();
}

/* Get a mapping of peptide symbols to HELM monomer library objects with
 * selectted fields  */
function getFormattedMonomerLib(monomersLibList: any[]): Map<string, any> {
  const map = new Map<string, any>();
  monomersLibList.forEach(
    (it) => {
      // todo: generalize for the case of nucleotides
      // todo: remove string literals
      if (it['polymerType'] === 'PEPTIDE') {
        const monomerObject: { [key: string]: any } = {};
        HELM_CORE_FIELDS.forEach((field) => {
          monomerObject[field] = it[field];
        });
        map.set(it[MONOMER_SYMBOL], monomerObject);
      }
    });
  return map;
}

/* Get jagged array of monomer symbols for the dataframe  */
function getMonomerSequencesArray(macroMolCol: DG.Column<string>): string[][] {
  const result: string[][] = new Array(macroMolCol.length);

  // todo: getTag, constants
  const colUnits = macroMolCol.tags[DG.TAGS.UNITS];
  // todo: separator to constant
  const separator = macroMolCol.getTag('separator');
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(colUnits, separator);
  for (let row = 0; row < macroMolCol.length; ++row) {
    const macroMolecule = macroMolCol.get(row);
    // todo: handle the exception case when macroMolecule can be null
    result[row] = macroMolecule ? splitterFunc(macroMolecule) : [];
  }
  return result;
}

/* Get a mapping of monomer symbols to MolGraph objects */
async function getMonomersDict(
  monomerSequencesArray: string[][], monomersLibList: any[]
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

/* Construct the MolGraph object for a specified monomerSymbol: the associated
 * graph is rotated and filled with default R-groups */
function getMolGraph(
  monomerSymbol: string, formattedMonomerLib: Map<string, any>, moduleRdkit: any // todo: specify type
): MolGraph | null {
  if (!formattedMonomerLib.has(monomerSymbol)) {
    return null;
  } else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    // todo: create constants for HELM molfile fields
    const rGroups = parseRGroupTypes(libObject[RGROUP_FIELD]);
    const rGroupIndices = parseRGroupIndices(libObject[MOLFILE_FIELD]);
    const molfileV2K = substituteRGroups(libObject[MOLFILE_FIELD], rGroups, rGroupIndices);
    const molfileV3K = convertMolfileToV3K(removeRGroupLine(molfileV2K), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);
    const bondData = parseBondData(molfileV3K, counts.bondCount);
    const removedNodes = getRemovedNodes(rGroupIndices);
    const atomData = parseAtomData(removedNodes, bondData, molfileV3K, counts.atomCount);

    removeIntermediateHydrogen(atomData, bondData);
    adjustGraph(atomData, atomData.leftNode, atomData.rightRemovedNode);
    // todo: delete after nucleotides debugging
    if (atomData.leftNode === atomData.leftRemovedNode || atomData.rightNode === atomData.rightRemovedNode)
      throw new Error('Non-intermediate hydrogens removed');
    return {atoms: atomData, bonds: bondData};
  }
}

/* Parse element symbols for R-groups from the HELM monomer library R-groups
 * field  */
function parseRGroupTypes(rgroupObjList: any[]): string[] {
  // specifically for HELMCoreLibrary
  // considered only monoatomic rgroups
  // supposing that elements in rgroupObjList are sorted w.r.t. the rgroups idx
  // todo: possible generalizations
  const rgroupsArray: string[] = [];
  for (const obj of rgroupObjList) {
    let rgroup: string = obj[CAP_GROUP_SMILES];
    // todo: verify that there are no multi-element rgroups, or consider how to
    // transform them
    rgroup = rgroup.replace(/(\[|\]|\*|:|\d)/g, '');
    if (rgroup.length > 1) // todo: check if such cases are possible, remove if not
      throw new Error('Default r-group has length more than one');
    rgroupsArray.push(rgroup);
  }
  return rgroupsArray;
}

/* Substitute the element symbols instead of R# in molfile  */
function substituteRGroups(molfileV2K: string, rGroups: string[], rGroupIndices: number[]) {
  let modifiedMolfile = molfileV2K;
  const idx: number[] = new Array((rGroupIndices.length - 1)/2);
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

/* Helper function necessary to build a correct V3000 molfile after substitution of r-groups in V2000 molfile  */
function removeRGroupLine(molfileV2K: string): string {
  // todo: consider the case of multiple rgp lines in a molfile
  const begin = molfileV2K.indexOf(V2K_RGP_LINE);
  const end = molfileV2K.indexOf('\n', begin) + 1;
  return molfileV2K.substring(0, begin) + molfileV2K.substring(end);
}

/* V2000 to V3000 converter  */
function convertMolfileToV3K(molfileV2K: string, moduleRdkit: any): string {
  // todo: type of moduleRdkit
  // todo: consider the use of unified converter (relies on creation of moduleRdkit
  // on each iteration, though)
  const molObj = moduleRdkit.get_mol(molfileV2K);
  const molfileV3K = molObj.get_v3Kmolblock();
  molObj.delete();
  return molfileV3K;
}

function parseBondData(molfileV3K: string, bondCount: number): BondData {
  // todo: consider the case when there is no simple leftmost/rightmost choice
  // todo: consider the case when there are multiple consequent M  RGP lines,
  // like in HELMCoreLibrary nucleotides

  const bondType: number[] = new Array(bondCount);
  const atomPair: number[][] = new Array(bondCount);
  const kwargs: string[] = new Array(bondCount);

  let begin = molfileV3K.indexOf(V3K_BEGIN_BOND_BLOCK);
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;
  for (let i = 0; i < bondCount; ++i) {
    const parsedValues: number[] = new Array(3);
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

/* Get the positions of nodes removed upon chaining the monomers into a polymer  */
function getRemovedNodes(rGroupIndices: number[]): {left: number, right: number} {
  // todo: verify that in other monomer libraries the order of removed nodes in
  // RGP line is the same as in HELMCoreLibrary
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
  atomPair: number[][], removedNodes: {left: number, right: number}
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

  return {left: leftNode, right: rightNode};
}

function parseAtomAndBondCounts(molfileV3K: string): {atomCount: number, bondCount: number} {
  molfileV3K = molfileV3K.replaceAll('\r', ''); // to handle old and new sdf standards

  // parse atom count
  let begin = molfileV3K.indexOf(V3K_BEGIN_COUNTS_LINE) + V3K_COUNTS_SHIFT;
  let end = molfileV3K.indexOf(' ', begin);
  const numOfAtoms = parseInt(molfileV3K.substring(begin, end));

  // parse bond count
  begin = end + 1;
  end = molfileV3K.indexOf(' ', begin);
  const numOfBonds = parseInt(molfileV3K.substring(begin, end));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

function parseAtomData(
  removedNodes: {left: number, right: number}, bondData: BondData, molfileV3K: string, atomCount: number
): AtomData {
  const atomType: string[] = new Array(atomCount);
  const x: number[] = new Array(atomCount);
  const y: number[] = new Array(atomCount);
  const kwargs: string[] = new Array(atomCount);

  const atomPair = bondData.atomPair;
  const retainedNodes = parseRetainedLRNodes(atomPair, removedNodes);

  let begin = molfileV3K.indexOf(V3K_BEGIN_ATOM_BLOCK); // V3000 atoms block
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;

  for (let i = 0; i < atomCount; i++) {
    begin = molfileV3K.indexOf(V3K_BEGIN_DATA_LINE, begin) + V3K_IDX_SHIFT;
    end = molfileV3K.indexOf(' ', begin); // skip the idx row

    // parse atom type
    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin);
    atomType[i] = molfileV3K.substring(begin, end);

    // parse X and Y coordinates of the atom
    const coordinate: number[] = new Array(2);
    for (let k = 0; k < 2; ++k) {
      begin = end + 1;
      end = molfileV3K.indexOf(' ', begin);
      coordinate[k] = parseFloat(molfileV3K.substring(begin, end));
    }
    x[i] = coordinate[0];
    y[i] = coordinate[1];

    // parse the remaining possible keyword arguments
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

/* Remove hydrogen nodes that are not left/right removed nodes */
function removeIntermediateHydrogen(atomData: AtomData, bondData: BondData): void {
  let i = 0;
  const leftRemovedIdx = atomData.leftRemovedNode - 1;
  const rightRemovedIdx = atomData.rightRemovedNode - 1;
  while (i < atomData.atomType.length) {
    if ( atomData.atomType[i] === 'H' &&
      i !== leftRemovedIdx &&
      i !== rightRemovedIdx) {
      removeNodeAndBonds(atomData, bondData, i + 1); // i + 1 because molfile node indexing starts from 1
      --i;
    }
    ++i;
  }
}

/* Remove H atoms that are L/R removed nodes */
function removeTerminalHydrogen(molGraph: MolGraph, removedNode: number): MolGraph {
  const molGraphCopy: MolGraph = structuredClone(molGraph);
  if (molGraphCopy.atoms.atomType[removedNode - 1] === 'H')
    removeNodeAndBonds(molGraphCopy.atoms, molGraphCopy.bonds, removedNode);
  return molGraphCopy;
}

/* Remove node 'removedNode' and the associated bonds. Notice, numeration of
 * nodes in molfiles starts from 1, not 0 */
function removeNodeAndBonds(atomData: AtomData, bondData: BondData, removedNode: number): void {
  const removedNodeIdx = removedNode - 1;

  // remove the node from atomData
  atomData.atomType.splice(removedNodeIdx, 1);
  atomData.x.splice(removedNodeIdx, 1);
  atomData.y.splice(removedNodeIdx, 1);
  atomData.kwargs.splice(removedNodeIdx, 1);

  // update the values of L/R retained and removed nodes if necessary
  atomData.leftNode = (atomData.leftNode > removedNode) ? --atomData.leftNode : atomData.leftNode;
  atomData.rightNode = (atomData.rightNode > removedNode) ? --atomData.rightNode : atomData.rightNode;
  atomData.leftRemovedNode = (atomData.leftRemovedNode > removedNode) ?
    --atomData.leftRemovedNode : atomData.leftRemovedNode;
  atomData.rightRemovedNode = (atomData.rightRemovedNode > removedNode) ?
    --atomData.rightRemovedNode : atomData.rightRemovedNode;

  // update indices of atoms in bonds
  let i = 0;
  while (i < bondData.atomPair.length) {
    const firstAtom = bondData.atomPair[i][0];
    const secondAtom = bondData.atomPair[i][1];
    if (firstAtom === removedNode || secondAtom === removedNode) {
      bondData.atomPair.splice(i, 1);
      bondData.bondType.splice(i, 1);
      bondData.kwargs.splice(i, 1);
      --i;
    } else {
      bondData.atomPair[i][0] = (firstAtom > removedNode) ? firstAtom - 1 : firstAtom;
      bondData.atomPair[i][1] = (secondAtom > removedNode) ? secondAtom - 1 : secondAtom;
    }
    ++i;
  }
}

/* Adjust the monomer's graph so that nodeOne is at origin, while nodeTwo
 * is at the positive ray of OX */
function adjustGraph(atoms: AtomData, nodeOne: number, nodeTwo: number): void {
  // todo: consider the case when d(l, r) < max(d)
  const nodeOneIdx = nodeOne - 1; // node indexing in molfiles starts from 1
  const nodeTwoIdx = nodeTwo - 1;

  // center the backbone
  const xCenter = (atoms.x[nodeOneIdx] + atoms.x[nodeTwoIdx])/2;
  const yCenter = (atoms.y[nodeOneIdx] + atoms.y[nodeTwoIdx])/2;
  shiftCoordinates(atoms, -xCenter, -yCenter);

  // rotate left node around the center
  rotateCenteredGraph(atoms, nodeOne);

  // shift the backbone so that nodeOne is at origin
  const xShift = atoms.x[nodeTwoIdx]; // equal to -atoms.x[nodeOneIdx]
  shiftCoordinates(atoms, xShift);
  // todo: remove after debugging
  if (atoms.y[nodeTwoIdx] !== 0)
    throw new Error('wrong rotation');
}

/*  Rotate 'node' around the origin, so that it ends up on OX's negative ray */
function rotateCenteredGraph(atoms: AtomData, node: number): void {
  let angle = 0;
  const nodeIdx = node - 1;
  const x = atoms.x;
  const y = atoms.y;
  const xRotated = x[nodeIdx];
  const yRotated = y[nodeIdx];

  if (xRotated === 0) { // the rotated node is on OY
    angle = yRotated > 0 ? Math.PI/2 : -Math.PI/2;
  } else if (yRotated === 0) {
    angle = xRotated > 0 ? Math.PI : 0;
  } else {
    const tangent = yRotated / xRotated;
    if (xRotated < 0) {
      angle = -Math.atan(tangent);
    } else {
      angle = tangent > 0 ?
        Math.PI - Math.atan(tangent) :
        -(Math.atan(tangent) + Math.PI);
    }
  }

  const cos = Math.cos(angle);
  const sin = Math.sin(angle);

  // rotation
  for (let i = 0; i < x.length; ++i) {
    const tmp = x[i];
    x[i] = Math.round(10_000 * (tmp*cos - y[i]*sin))/10_000;
    y[i] = Math.round(10_000 * (tmp*sin + y[i]*cos))/10_000;
  }
  // todo: remove after debugging
  if (atoms.x[nodeIdx] > 0)
    throw new Error(`Incorrect rotation: ${atoms.atomType[nodeIdx]}`);
}

/* Conctenate MolGraph objects according to the monomer symbol sequence, the
 * result is convertible to molfile V3000 */
function monomerSeqToMolGraph(monomerSeq: string[], monomersDict: Map<string, MolGraph>): MolGraph {
  if (monomerSeq.length === 0)
    throw new Error('The monomerSeq is empty');
    // todo: handle/modify exception

  // remove "left" hydrogens if any
  let result = removeTerminalHydrogen(
    monomersDict.get(monomerSeq[0])!,
    monomersDict.get(monomerSeq[0])!.atoms.leftRemovedNode
  );

  concatenateMonomerGraphs(result, monomerSeq, monomersDict);

  result = removeTerminalHydrogen(result, result.atoms.rightRemovedNode);
  return result;
}

/* Join two subsequent MolGraph objects */
function concatenateMonomerGraphs(result: MolGraph, monomerSeq: string[], monomersDict: Map<string, MolGraph>): void {
  // todo: improve naming
  for (let i = 1; i < monomerSeq.length; i++) {
    // todo: choose one of the shift types
    const xShift = result.atoms.x[result.atoms.rightRemovedNode - 1];
    // const yShift = result.atoms.y[result.atoms.rightRemovedNode - 1];
    removeNodeAndBonds(result.atoms, result.bonds, result.atoms.rightRemovedNode);

    const nextMonomer: MolGraph = structuredClone(monomersDict.get(monomerSeq[i])!);
    removeNodeAndBonds(nextMonomer.atoms, nextMonomer.bonds, nextMonomer.atoms.leftRemovedNode);

    const nodeIdxShift = result.atoms.atomType.length;
    shiftNodes(nextMonomer, nodeIdxShift);
    // shiftCoordinates(nextMonomer, xShift, yShift);
    shiftCoordinates(nextMonomer.atoms, xShift);

    // todo: improve chaining by considering max-area rectangles associated with
    // the monomers and making those non-overlapping, pairs of such rectangles
    // should be taken into account. If the rectangles intersect, flip the
    // current monomer (mirror w.r.t. OX). This is possible if we are interested
    // in preserving the metric relations

    // bind the two monomers with a peptide bond
    // todo: eliminate hardcoding of the peptide bond
    result.bonds.bondType.push(1);
    result.bonds.atomPair.push([result.atoms.rightNode, nextMonomer.atoms.leftNode]);
    result.bonds.kwargs.push('\n');

    // update result with the values of nextMonomer
    result.atoms.rightNode = nextMonomer.atoms.rightNode;
    result.atoms.rightRemovedNode = nextMonomer.atoms.rightRemovedNode;
    result.atoms.atomType = result.atoms.atomType.concat(nextMonomer.atoms.atomType);
    result.atoms.x = result.atoms.x.concat(nextMonomer.atoms.x);
    result.atoms.y = result.atoms.y.concat(nextMonomer.atoms.y);
    result.atoms.kwargs = result.atoms.kwargs.concat(nextMonomer.atoms.kwargs);

    result.bonds.kwargs = result.bonds.kwargs.concat(nextMonomer.bonds.kwargs);
    result.bonds.atomPair = result.bonds.atomPair.concat(nextMonomer.bonds.atomPair);
    result.bonds.bondType = result.bonds.bondType.concat(nextMonomer.bonds.bondType);
  }
}

/* Shift the molfile indices of the nodes  */
function shiftNodes(monomer: MolGraph, nodeIdxShift: number): void {
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

/* Shift the backbone's Cartesian coordinates, yShift is optional  */
function shiftCoordinates(atoms: AtomData, xShift: number, yShift?: number): void {
  const x = atoms.x;
  const y = atoms.y;
  for (let i = 0; i < x.length; ++i) {
    x[i] += xShift;
    if (typeof yShift !== 'undefined')
      y[i] += yShift;
  }
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

  // todo rewrite using constants
  const molfileCountsLine = V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + V3K_COUNTS_LINE_ENDING;

  // atom block
  let molfileAtomBlock = '';
  for (let i = 0; i < atomCount; ++i) {
    const atomIdx = i + 1;
    const coordinate = [x[i].toString(), y[i].toString()];
    for (let k = 0; k < 2; ++k) {
      const formatted = coordinate[k].toString().split('.');
      if (formatted.length === 1)
        formatted.push('0');
      formatted[1] = formatted[1].padEnd(V3K_ATOM_COORDINATE_PRECISION, '0');
      coordinate[k] = formatted.join('.');
    }
    const atomLine = V3K_BEGIN_DATA_LINE + atomIdx + ' ' + atomType[i] + ' ' +
      coordinate[0] + ' ' + coordinate[1] + ' ' + atomKwargs[i];
    molfileAtomBlock += atomLine;
  }

  // bond block
  let molfileBondBlock = '';
  for (let i = 0; i < bondCount; ++i) {
    const bondIdx = i + 1;
    const firstAtom = atomPair[i][0];
    const secondAtom = atomPair[i][1];
    const bondLine = V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
      firstAtom + ' ' + secondAtom + ' ' + bondKwargs[i];
    molfileBondBlock += bondLine;
  }

  const molfileParts = [
    V3K_HEADER_FIRST_LINE,
    V3K_HEADER_SECOND_LINE,
    V3K_BEGIN_CTAB_BLOCK,
    molfileCountsLine,
    V3K_BEGIN_ATOM_BLOCK,
    molfileAtomBlock,
    V3K_END_ATOM_BLOCK,
    V3K_BEGIN_BOND_BLOCK,
    molfileBondBlock,
    V3K_END_BOND_BLOCK,
    V3K_END_CTAB_BLOCK,
    V3K_END,
  ];

  return molfileParts.join('');
}
