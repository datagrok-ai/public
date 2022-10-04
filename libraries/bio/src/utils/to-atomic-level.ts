/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '../../src/viewers/web-logo';
import {HELM_FIELDS, HELM_CORE_FIELDS, HELM_POLYMER_TYPE, RGROUP_FIELDS} from './const';
import {ALPHABET, NOTATION, UnitsHandler} from './units-handler';
import {NotationConverter} from './notation-converter';

// constants for parsing of molfile V2000
const V2K_RGP_SHIFT = 8;
const V2K_RGP_LINE = 'M  RGP';

// constants for parsing/reconstruction of molfile V3000
const V3K_COUNTS_SHIFT = 14;
const V3K_IDX_SHIFT = 7;
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
const V3K_BOND_CONFIG = ' CFG=';
const V3K_BEGIN_DATA_LINE = 'M  V30 ';
const V3K_END = 'M  END\n';

const PRECISION_FACTOR = 10_000; // HELMCoreLibrary has 4 significant digits after decimal point in atom coordinates

const HELM_DEOXYRIBOSE = 'd';
const HELM_RIBOSE = 'r';
const HELM_PHOSPHATE = 'p';

type AtomData = {
  atomType: string[], // element symbol
  x: number[], // Cartesian coordiantes of the nodes
  y: number[],
  leftAttachmentNode: number, // "leftmost" attachment node, e.g. N-terminus in peptides
  rightAttachmentNode: number,
  leftRNode: number, // "leftmost" R-group node of the monomer (removed upon chaining into polymer)
  rightRNode: number,
  kwargs: string[], // MDLV30 atom line may contain keyword args
  shift: number[], // shift between leftAttachmentNode and rightRNode
}

type BondData = {
  bondType: number[],
  atomPair: number[][], // indices of atoms, starting with 1
  bondConfiguration: Map<number, number>,
  kwargs: Map<number, string>, // MDLV30 atom line may contain keyword args
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

  // todo: remove progress indicator?
  // const pi = DG.TaskBarProgressIndicator.create('Restoring atomic structure...');

  // verify the correctness of the semtype
  if (macroMolCol.semType !== DG.SEMTYPE.MACROMOLECULE) {
    grok.shell.warning(
      `Only the ${DG.SEMTYPE.MACROMOLECULE} columns can be converted to atomic
      level, the chosen column has semType ${macroMolCol.semType}`
    );
    return;
  }

  // convert 'helm' to 'separator' units
  if (macroMolCol.getTag(DG.TAGS.UNITS) === NOTATION.HELM) {
    const converter = new NotationConverter(macroMolCol);
    const separator = '/';
    macroMolCol = converter.convert(NOTATION.SEPARATOR, separator);
  }

  const alphabet = macroMolCol.getTag(UnitsHandler.TAGS.alphabet);
  // console.log(alphabet);

  // determine the polymer type according to HELM specifications
  let polymerType;
  // todo: an exception from dart comes before this check if the alphabet is UN
  if (alphabet === ALPHABET.PT) {
    polymerType = HELM_POLYMER_TYPE.PEPTIDE;
  } else if (alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA) {
    polymerType = HELM_POLYMER_TYPE.RNA;
  } else {
    // todo: improve the column selection mechanism
    grok.shell.warning(
      `Only PT, DNA and RNA alphabets are supported, the selected column has ${polymerType} alphabet`
    );
    return;
  }

  const monomerSequencesArray: string[][] = getMonomerSequencesArray(macroMolCol);
  // todo: consider separately backbone, terminal, branch monomer types
  const monomersDict = await getMonomersDict(monomerSequencesArray, monomersLibList, polymerType, alphabet);
  const columnLength = macroMolCol.length;
  const reconstructed: string[] = new Array(columnLength);
  for (let row = 0; row < columnLength; ++row) {
    const monomerSeq = monomerSequencesArray[row];
    reconstructed[row] = monomerSeqToMolfile(monomerSeq, monomersDict, polymerType);
    // console.log(reconstructed[row]);
    // if (row % 100 === 0)
    // pi.update(100 * (row + 1)/columnLength, '...complete');
  }

  // exclude name collisions
  const name = 'molfile(' + macroMolCol.name + ')';
  const newColName = df.columns.getUnusedName(name);
  const newCol = DG.Column.fromStrings(newColName, reconstructed);

  newCol.semType = DG.SEMTYPE.MOLECULE;
  newCol.setTag(DG.TAGS.UNITS, DG.UNITS.Molecule.MOLBLOCK);
  df.columns.add(newCol, true);
  await grok.data.detectSemanticTypes(df);

  // pi.close();
}

/* Get a mapping of peptide symbols to HELM monomer library objects with
 * selectted fields  */
function getFormattedMonomerLib(monomersLibList: any[], polymerType: HELM_POLYMER_TYPE): Map<string, any> {
  const map = new Map<string, any>();
  monomersLibList.forEach(
    (it) => {
      // filter out the monomers that are of the type we need
      if (it[HELM_FIELDS.POLYMER_TYPE] === polymerType) {
        const monomerObject: { [key: string]: any } = {};
        HELM_CORE_FIELDS.forEach((field) => {
          monomerObject[field] = it[field];
        });
        map.set(it[HELM_FIELDS.SYMBOL], monomerObject);
      }
    });
  return map;
}

/* Get jagged array of monomer symbols for the dataframe  */
function getMonomerSequencesArray(macroMolCol: DG.Column<string>): string[][] {
  const columnLength = macroMolCol.length;
  const result: string[][] = new Array(columnLength);

  // split the string into monomers
  const colUnits = macroMolCol.getTag(DG.TAGS.UNITS);
  const separator = macroMolCol.getTag(UnitsHandler.TAGS.separator);
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(colUnits, separator);

  for (let row = 0; row < columnLength; ++row) {
    const macroMolecule = macroMolCol.get(row);
    // todo: handle the exception case when macroMolecule is null
    result[row] = macroMolecule ? splitterFunc(macroMolecule) : [];
  }
  return result;
}

/* Get a mapping of monomer symbols to MolGraph objects */
async function getMonomersDict(
  monomerSequencesArray: string[][], monomersLibList: any[], polymerType: HELM_POLYMER_TYPE, alphabet: ALPHABET
): Promise<Map<string, MolGraph>> {
  const formattedMonomerLib = getFormattedMonomerLib(monomersLibList, polymerType);
  const monomersDict = new Map<string, MolGraph>();

  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');

  for (let row = 0; row < monomerSequencesArray.length; ++row) {
    const monomerSeq: string[] = monomerSequencesArray[row];
    for (const sym of monomerSeq)
      updateMonomersDict(monomersDict, sym, formattedMonomerLib, moduleRdkit, polymerType);
  }

  // add deoxyribose/ribose and phosphate for nucleotide sequences
  if (polymerType === HELM_POLYMER_TYPE.RNA) {
    const symbols = (alphabet === ALPHABET.RNA) ?
      [HELM_RIBOSE, HELM_PHOSPHATE] : [HELM_DEOXYRIBOSE, HELM_PHOSPHATE];
    for (const sym of symbols)
      updateMonomersDict(monomersDict, sym, formattedMonomerLib, moduleRdkit, polymerType);
  }

  return monomersDict;
}

function updateMonomersDict(
  monomersDict: Map<string, MolGraph>, sym: string, formattedMonomerLib: Map<string, any>,
  moduleRdkit: any, polymerType: HELM_POLYMER_TYPE
) {
  if (!monomersDict.has(sym)) {
    const monomerData: MolGraph | null = getMolGraph(sym, formattedMonomerLib, moduleRdkit, polymerType);
    if (monomerData)
      monomersDict.set(sym, monomerData);
    // todo: handle exception when there is no monomer with symbol sym in
    // monomersLibList
  }
}

/* Construct the MolGraph object for a specified monomerSymbol: the associated
 * graph is rotated and filled with default R-groups */
function getMolGraph(
  monomerSymbol: string, formattedMonomerLib: Map<string, any>,
  moduleRdkit: any, polymerType: HELM_POLYMER_TYPE // todo: specify type for moduleRdkit
): MolGraph | null {
  if (!formattedMonomerLib.has(monomerSymbol)) {
    return null;
  } else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    const rGroups = parseRGroupTypes(libObject[HELM_FIELDS.RGROUPS]);
    const rGroupIndices = parseRGroupIndices(libObject[HELM_FIELDS.MOLFILE]);
    const molfileV2K = substituteRGroups(libObject[HELM_FIELDS.MOLFILE], rGroups, rGroupIndices);
    const molfileV3K = convertMolfileToV3K(removeRGroupLine(molfileV2K), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);
    const bondData = parseBondData(molfileV3K, counts.bondCount);
    const removedNodes = getRemovedNodes(rGroupIndices);
    const atomData = parseAtomData(removedNodes, bondData, molfileV3K, counts.atomCount);

    const monomerGraph = {atoms: atomData, bonds: bondData};
    adjustPeptideMonomerGraph(monomerGraph);
    monomerGraph.atoms.shift = [
      keepPrecision(
        monomerGraph.atoms.x[monomerGraph.atoms.rightRNode - 1] -
        monomerGraph.atoms.x[monomerGraph.atoms.leftAttachmentNode - 1]
      ),
      keepPrecision(
        monomerGraph.atoms.y[monomerGraph.atoms.rightRNode - 1] -
        monomerGraph.atoms.y[monomerGraph.atoms.leftAttachmentNode - 1]
      ),
    ];
    removeHydrogen(monomerGraph);
    removeNodeAndBonds(monomerGraph, monomerGraph.atoms.rightRNode);

    return monomerGraph;
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
    let rgroup: string = obj[RGROUP_FIELDS.CAP_GROUP_SMILES];
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
    // todo: remove after debugging
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
  // todo: consider the use of standard Chem converter (relies on creation of moduleRdkit on each iteration, though)
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
  const bondConfiguration = new Map<number, number>();
  const kwargs = new Map<number, string>;

  let begin = molfileV3K.indexOf(V3K_BEGIN_BOND_BLOCK);
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;
  for (let i = 0; i < bondCount; ++i) {
    // parse bond type and atom pair
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

    // parse keyword arguments
    const endOfLine = molfileV3K.indexOf('\n', begin);
    let lineRemainder = molfileV3K.slice(end, endOfLine);
    let beginCfg = lineRemainder.indexOf(V3K_BOND_CONFIG);
    if (beginCfg !== -1) {
      beginCfg = lineRemainder.indexOf('=', beginCfg) + 1;
      let endCfg = lineRemainder.indexOf(' ', beginCfg);
      if (endCfg === -1)
        endCfg = lineRemainder.length;
      const bondConfig = parseInt(lineRemainder.slice(beginCfg, endCfg));
      bondConfiguration.set(i, bondConfig);
      const removedSubstring = V3K_BOND_CONFIG + bondConfig.toString();
      lineRemainder = lineRemainder.replace(removedSubstring, '');
    }
    if (!lineRemainder)
      kwargs.set(i, lineRemainder);
  }

  return {
    bondType: bondType,
    atomPair: atomPair,
    bondConfiguration: bondConfiguration,
    kwargs: kwargs,
  };
}

/* Get the positions of nodes removed upon chaining the monomers into a polymer  */
function getRemovedNodes(rGroupIndices: number[]): {left: number, right: number} {
  // todo: verify that in other monomer libraries the order of removed nodes in
  // RGP line is the same as in HELMCoreLibrary
  const leftRNode = rGroupIndices[1]; // rgpIndices[0] is the number of R-groups
  const rightRNode = rGroupIndices[rGroupIndices.length - 2];
  return {left: leftRNode, right: rightRNode};
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
  let leftAttachmentNode = 0;
  let rightAttachmentNode = 0;

  let i = 0;
  while ((i < atomPair.length) && (leftAttachmentNode === 0 || rightAttachmentNode === 0)) {
    if (atomPair[i][0] === removedNodes.left)
      leftAttachmentNode = atomPair[i][1];
    else if (atomPair[i][1] === removedNodes.left)
      leftAttachmentNode = atomPair[i][0];
    else if (atomPair[i][0] === removedNodes.right)
      rightAttachmentNode = atomPair[i][1];
    else if (atomPair[i][1] === removedNodes.right)
      rightAttachmentNode = atomPair[i][0];
    ++i;
  }

  return {left: leftAttachmentNode, right: rightAttachmentNode};
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
    leftAttachmentNode: retainedNodes.left,
    rightAttachmentNode: retainedNodes.right,
    leftRNode: removedNodes.left,
    rightRNode: removedNodes.right,
    kwargs: kwargs,
    shift: [],
  };
}

/* Remove hydrogen nodes */
function removeHydrogen(monomerGraph: MolGraph): void {
  let i = 0;
  while (i < monomerGraph.atoms.atomType.length) {
    if ( monomerGraph.atoms.atomType[i] === 'H') {
      removeNodeAndBonds(monomerGraph, i + 1); // i + 1 because molfile node indexing starts from 1
      --i;
    }
    ++i;
  }
}

/* Remove node 'removedNode' and the associated bonds. Notice, numeration of
 * nodes in molfiles starts from 1, not 0 */
function removeNodeAndBonds(monomerGraph: MolGraph, removedNode: number): void {
  const removedNodeIdx = removedNode - 1;
  const atomData = monomerGraph.atoms;
  const bondData = monomerGraph.bonds;

  // remove the node from atomData
  atomData.atomType.splice(removedNodeIdx, 1);
  atomData.x.splice(removedNodeIdx, 1);
  atomData.y.splice(removedNodeIdx, 1);
  atomData.kwargs.splice(removedNodeIdx, 1);

  // update the values of L/R retained and removed nodes if necessary
  atomData.leftAttachmentNode = (atomData.leftAttachmentNode > removedNode) ?
    --atomData.leftAttachmentNode : atomData.leftAttachmentNode;
  atomData.rightAttachmentNode = (atomData.rightAttachmentNode > removedNode) ?
    --atomData.rightAttachmentNode : atomData.rightAttachmentNode;
  atomData.leftRNode = (atomData.leftRNode > removedNode) ?
    --atomData.leftRNode : atomData.leftRNode;
  atomData.rightRNode = (atomData.rightRNode > removedNode) ?
    --atomData.rightRNode : atomData.rightRNode;

  // update indices of atoms in bonds
  let i = 0;
  while (i < bondData.atomPair.length) {
    const firstAtom = bondData.atomPair[i][0];
    const secondAtom = bondData.atomPair[i][1];
    if (firstAtom === removedNode || secondAtom === removedNode) {
      bondData.atomPair.splice(i, 1);
      bondData.bondType.splice(i, 1);
      if (bondData.bondConfiguration.has(i))
        bondData.bondConfiguration.delete(i);
      if (bondData.kwargs.has(i))
        bondData.kwargs.delete(i);
      --i;
    } else {
      bondData.atomPair[i][0] = (firstAtom > removedNode) ? firstAtom - 1 : firstAtom;
      bondData.atomPair[i][1] = (secondAtom > removedNode) ? secondAtom - 1 : secondAtom;
    }
    ++i;
  }

  // update bondConfiguration and kwargs keys
  let keys = Array.from(bondData.bondConfiguration.keys());
  keys.forEach((key) => {
    if (bondData.bondConfiguration.has(key) && key > removedNodeIdx) {
      const value = bondData.bondConfiguration.get(key)!;
      bondData.bondConfiguration.delete(key);
      bondData.bondConfiguration.set(key - 1, value);
    }
  });
  keys = Array.from(bondData.kwargs.keys());
  keys.forEach((key) => {
    if (bondData.kwargs.has(key) && key > removedNodeIdx) {
      const value = bondData.kwargs.get(key)!;
      bondData.kwargs.delete(key);
      bondData.kwargs.set(key - 1, value);
    }
  });
}

// todo: rewrite description
/* Adjust the (peptide) monomer graph so that it has standard form  */
function adjustPeptideMonomerGraph(monomer: MolGraph): void {
  const nodeOneIdx = monomer.atoms.leftAttachmentNode - 1; // node indexing in molfiles starts from 1
  const nodeTwoIdx = monomer.atoms.leftRNode - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // angle is measured between OY and the rotated node
  const angle = findAngleWithOY(x[nodeTwoIdx], y[nodeTwoIdx]);

  // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  rotateCenteredGraph(monomer.atoms, -angle);

  if (x[monomer.atoms.rightRNode - 1] < 0)
    flipMonomerAroundOY(monomer);

  const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);

  // flip carboxyl and R if necessary
  flipCarboxylAndRadical(monomer, doubleBondedOxygen);

  // flip hydroxyl group with double-bound O inside carboxyl group if necessary
  flipHydroxilGroup(monomer, doubleBondedOxygen);
}

function flipCarboxylAndRadical(monomer: MolGraph, doubleBondedOxygen: number): void {
  // verify that the carboxyl group is in the lower half-plane
  if (monomer.atoms.y[monomer.atoms.rightRNode - 1] < 0 &&
    monomer.atoms.y[doubleBondedOxygen - 1] < 0) {
    flipMonomerAroundOX(monomer);

    rotateCenteredGraph(monomer.atoms,
      -findAngleWithOX(
        monomer.atoms.x[monomer.atoms.rightAttachmentNode - 1],
        monomer.atoms.y[monomer.atoms.rightAttachmentNode - 1]
      )
    );
  }
}

/* Finds angle between OY and the ray joining origin with (x, y) */
function findAngleWithOY(x: number, y: number): number {
  let angle;
  if (x === 0) {
    angle = y > 0 ? 0 : Math.PI;
  } else if (y === 0) {
    angle = x > 0 ? -Math.PI/2 : Math.PI/2;
  } else {
    const tan = y / x;
    const atan = Math.atan(tan);
    angle = (x < 0) ? Math.PI/2 + atan : -Math.PI/2 + atan;
  }
  return angle;
}

/* Finds angle between OX and the ray joining origin with (x, y) */
function findAngleWithOX(x: number, y: number): number {
  return findAngleWithOY(x, y) + Math.PI/2;
}

/*  Rotate the graph around the origin by 'angle' */
function rotateCenteredGraph(atoms: AtomData, angle: number): void {
  if (angle !== 0) {
    const x = atoms.x;
    const y = atoms.y;

    const cos = Math.cos(angle);
    const sin = Math.sin(angle);

    for (let i = 0; i < x.length; ++i) {
      const tmp = x[i];
      x[i] = keepPrecision(tmp*cos - y[i]*sin);
      y[i] = keepPrecision(tmp*sin + y[i]*cos);
    }
  }
}

/* Flip monomer graph around OX axis preserving stereometry */
function flipMonomerAroundOX(monomer: MolGraph): void {
  flipMolGraph(monomer, true);
}

/* Flip monomer graph around OY axis preserving stereometry */
function flipMonomerAroundOY(monomer: MolGraph): void {
  flipMolGraph(monomer, false);
}

/* Flip graph around a specified axis: 'true' corresponds to OX, 'false' to OY */
function flipMolGraph(molGraph: MolGraph, axis: boolean): void {
  if (axis) { // flipping around OX
    const y = molGraph.atoms.y;
    for (let i = 0; i < y.length; i++)
      y[i] = -y[i];
  } else { // flipping around OY
    const x = molGraph.atoms.x;
    for (let i = 0; i < x.length; i++)
      x[i] = -x[i];
  }

  // preserve the stereometry
  const orientation = molGraph.bonds.bondConfiguration;
  for (const [key, value] of orientation) {
    const newValue = value === 1 ? 3 : 1;
    orientation.set(key, newValue);
  }
}

/* Flips double-bonded O in carbonyl group with OH in order for the monomers
 * to have standard representation simplifying their concatenation. The
 * monomer must already be adjusted with adjustPeptideMonomerGraph in order for this function to be implemented  */
function flipHydroxilGroup(monomer: MolGraph, doubleBondedOxygen: number): void {
  const x = monomer.atoms.x;
  // -1 below because indexing of nodes in molfiles starts from 1, unlike arrays
  if (x[monomer.atoms.rightRNode - 1] > x[doubleBondedOxygen - 1])
    swapNodes(monomer, doubleBondedOxygen, monomer.atoms.rightRNode);
}

/* Determine the number of node (starting from 1) corresponding to the
 * double-bonded oxygen of the carbonyl group  */
function findDoubleBondedCarbonylOxygen(monomer: MolGraph): number {
  const bondsMap = constructBondsMap(monomer);
  let doubleBondedOxygen = 0;
  let i = 0;
  // iterate over the nodes bonded to the carbon and find the double one
  while (doubleBondedOxygen === 0) {
    const node = bondsMap.get(monomer.atoms.rightAttachmentNode)![i];
    if (monomer.atoms.atomType[node - 1] === 'O' && node !== monomer.atoms.rightRNode)
      doubleBondedOxygen = node;
    i++;
  }
  return doubleBondedOxygen;
}

/* Swap the Cartesian coordinates of the two specified nodes in MolGraph  */
function swapNodes(monomer: MolGraph, nodeOne: number, nodeTwo: number): void {
  const nodeOneIdx = nodeOne - 1;
  const nodeTwoIdx = nodeTwo - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;
  const tmpX = x[nodeOneIdx];
  const tmpY = y[nodeOneIdx];
  x[nodeOneIdx] = x[nodeTwoIdx];
  y[nodeOneIdx] = y[nodeTwoIdx];
  x[nodeTwoIdx] = tmpX;
  y[nodeTwoIdx] = tmpY;
}

function constructBondsMap(monomer: MolGraph): Map<number, Array<number>> {
  const map = new Map<number, Array<number>>();
  for (const atomPair of monomer.bonds.atomPair) {
    for (let i = 0; i < 2; i++) {
      const key = atomPair[i];
      const value = atomPair[(i + 1)%2];
      if (map.has(key))
        map.get(key)?.push(value);
      else
        map.set(key, new Array<number>(1).fill(value));
    }
  }
  return map;
}

/* Shift the backbone's Cartesian coordinates, yShift is optional  */
function shiftCoordinates(molGraph: MolGraph, xShift: number, yShift?: number): void {
  const x = molGraph.atoms.x;
  const y = molGraph.atoms.y;
  for (let i = 0; i < x.length; ++i) {
    x[i] = keepPrecision(x[i] + xShift);
    if (typeof yShift !== 'undefined')
      y[i] = keepPrecision(y[i] + yShift);
  }
}

/* Translate a sequence of monomer symbols into Molfile V3000 */
function monomerSeqToMolfile(monomerSeq: string[], monomersDict: Map<string,
  MolGraph>, polymerType: HELM_POLYMER_TYPE): string {
  // todo: handle the case when the polymer is empty
  if (monomerSeq.length === 0)
    throw new Error('The monomerSeq is empty');

  // define atom and bond counts
  let atomCount = 0;
  let bondCount = 0;
  for (const monomerSymbol of monomerSeq) {
    const monomer = monomersDict.get(monomerSymbol)!;
    atomCount += monomer?.atoms.x.length;
    bondCount += (monomer.bonds.bondType.length + 1); // +1 because of C-NH
  }
  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE)
    atomCount += 1; // because of the terminal OH

  // create the arrays to store lines of the resulting molfile
  const molfileAtomBlock = new Array<string>(atomCount);
  const molfileBondBlock = new Array<string>(bondCount);

  const positionShift = new Array<number>(2).fill(0);
  let nodeShift = 0;
  let bondShift = 0;
  let attachNode = 0; // node to which the next monomer is attached
  let flipFactor = 0; // flip every odd monomer
  // todo: should we consider representations other than planar?

  for (let i = 0; i < monomerSeq.length; ++i) {
    const monomer = monomersDict.get(monomerSeq[i])!;
    flipFactor = (-1)**(i % 2); // to flip every even monomer over OX

    // construnct the lines of V3K molfile atom block
    for (let j = 0; j < monomer.atoms.atomType.length; ++j) {
      const atomIdx = nodeShift + j + 1;
      molfileAtomBlock[nodeShift + j] = V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
        monomer.atoms.atomType[j] + ' ' +
        keepPrecision(positionShift[0] + monomer.atoms.x[j]) + ' ' +
        keepPrecision(positionShift[1] + flipFactor * monomer.atoms.y[j]) +
        ' ' + monomer.atoms.kwargs[j];
    }

    // construct the lines of V3K molfile bond block
    for (let j = 0; j < monomer.bonds.atomPair.length; ++j) {
      const bondIdx = bondShift + j + 1;
      const firstAtom = monomer.bonds.atomPair[j][0] + nodeShift;
      const secondAtom = monomer.bonds.atomPair[j][1] + nodeShift;
      let bondCfg = '';
      if (monomer.bonds.bondConfiguration.has(j)) {
        // flip orientation when necessary
        let orientation = monomer.bonds.bondConfiguration.get(j);
        if (flipFactor < 0)
          orientation = (orientation === 1) ? 3 : 1;
        bondCfg = ' CFG=' + orientation;
      }
      const kwargs = monomer.bonds.kwargs.has(j) ?
        ' ' + monomer.bonds.kwargs.get(j) : '';
      molfileBondBlock[bondShift + j] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
        monomer.bonds.bondType[j] + ' ' +
      firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
      //console.log(`Reconstructed bond block at ${bondShift + j}:` + molfileBondBlock[bondShift + j]);
    }

    // todo: refactor, this predicate should only be checked once for the column
    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
      // hardcoded peptide bond
      if (attachNode !== 0) {
        const bondIdx = bondShift;
        const firstAtom = attachNode;
        const secondAtom = monomer.atoms.leftAttachmentNode + nodeShift;
        molfileBondBlock[bondShift - 1] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
          1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
        // console.log(`Reconstructed peptide bond at ${bondShift}:` + molfileBondBlock[bondShift]);
      }
    }

    attachNode = nodeShift + monomer.atoms.rightAttachmentNode;
    bondShift += monomer.bonds.atomPair.length + 1;

    nodeShift += monomer.atoms.atomType.length;
    positionShift[0] += monomer.atoms.shift[0];
    positionShift[1] += flipFactor * monomer.atoms.shift[1];
  }

  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
    // add terminal oxygen
    const atomIdx = nodeShift + 1;
    molfileAtomBlock[atomCount] = V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
      'O' + ' ' + keepPrecision(positionShift[0]) + ' ' +
      flipFactor * keepPrecision(positionShift[1]) + ' ' + '0.000000 0' + '\n';

    // add terminal bond
    const firstAtom = attachNode;
    const secondAtom = atomIdx;
    molfileBondBlock[bondCount] = V3K_BEGIN_DATA_LINE + bondShift + ' ' +
      1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
  }

  // todo rewrite using constants
  const molfileCountsLine = V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + V3K_COUNTS_LINE_ENDING;

  const molfileParts = [
    V3K_HEADER_FIRST_LINE,
    V3K_HEADER_SECOND_LINE,
    V3K_BEGIN_CTAB_BLOCK,
    molfileCountsLine,
    V3K_BEGIN_ATOM_BLOCK,
    molfileAtomBlock.join(''),
    V3K_END_ATOM_BLOCK,
    V3K_BEGIN_BOND_BLOCK,
    molfileBondBlock.join(''),
    V3K_END_BOND_BLOCK,
    V3K_END_CTAB_BLOCK,
    V3K_END,
  ];

  return molfileParts.join('');
}

/* Keep precision upon floating point operations over atom coordinates */
function keepPrecision(x: number) {
  return Math.round(PRECISION_FACTOR * x)/PRECISION_FACTOR;
}
