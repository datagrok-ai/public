/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HELM_FIELDS, HELM_CORE_FIELDS, HELM_POLYMER_TYPE, HELM_MONOMER_TYPE, RGROUP_FIELDS, MODE} from '../utils/const';
import {ALPHABET, getSplitter, NOTATION, SplitterFunc, TAGS} from '../utils/macromolecule';
// import {UnitsHandler} from '../utils/units-handler';
import {NotationConverter} from '../utils/notation-converter';
import {Monomer} from '../types';

// constants for parsing molfile V2000
const V2K_RGP_SHIFT = 8;
const V2K_RGP_LINE = 'M  RGP';
const V2K_A_LINE = 'A  ';

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

// symbols for the corresponding monomers in HELM library
const DEOXYRIBOSE = 'd';
const RIBOSE = 'r';
const PHOSPHATE = 'p';

const OXYGEN = 'O';
const HYDROGEN = 'H';

/* Stores necessary data about atoms of a monomer parsed from Molfile */
type Atoms = {
  /* element symbols for monomer's atoms */
  atomTypes: string[],
  /* Cartesian coordiantes of monomer's atoms */
  x: number[], // todo: convert to Float32
  y: number[],
  /* V3K atom line may contain keyword args */
  kwargs: string[],
}

/* Stores necessary data about bonds of a monomer parsed from Molfile */
type Bonds = {
  /* bond types for all lines of Molfile bond block */
  bondTypes: number[], // todo: convert to Ind32
  /* Indices of all atom pairs, indexing starting from 1  */
  atomPairs: number[][],
  /* If a bond has CFG=... keyword argument, it is parsed and sotred as a
   * value of the map, with the key being the bond's index  */
  bondConfiguration: Map<number, number>,
  /* V3K bond line may contain keyword args */
  kwargs: Map<number, string>,
}

/* Metadata associated with the monomer necessary to restore the resulting
 * molfile  */
type MonomerMetadata = {
  /* terminal nodes: 0-th corresponds to the "leftmost" one, 1st, to the "rightmost",
   * e.g. N-terminus and C-terminus in peptides */
  terminalNodes: number[],
  /* r-group nodes: 0-th corresponds to the "leftmost" one, 1st, to the "rightmost" */
  rNodes: number[],
  /* shift from the origin to the next backbone, null for branch monomers */
  backboneShift: number[] | null,
  /* shift from the origin to the next branch, null for branch monomers */
  branchShift: number[] | null
}

type MolGraph = {
  atoms: Atoms,
  bonds: Bonds,
  meta: MonomerMetadata,
}

/* Helper structure wrapping common arguments to several functions */
type LoopVariables = {
  i: number,
  nodeShift: number,
  bondShift: number,
  backbonePositionShift: number[],
  backboneAttachNode: number; // node to which the next backbone is attached
  branchPositionShift: number[],
  branchAttachNode: number,
  flipFactor: number,
  // todo: should we consider representations other than planar?
}

/* Helper structure wrapping common arguments to several functions */
type LoopConstants = {
  sugar: MolGraph | null,
  phosphate: MolGraph | null,
  seqLength: number,
  atomCount: number,
  bondCount: number,
}

// todo: verify that all functions have return types

/* Convert Macromolecule column into Molecule column storing molfile V3000 with the help of a monomer library  */
export async function _toAtomicLevel(
  df: DG.DataFrame, macroMolCol: DG.Column<string>, monomersLibList: any[]
): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }

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

  const alphabet = macroMolCol.getTag(TAGS.alphabet);

  // determine the polymer type according to HELM specifications
  let polymerType;
  // todo: an exception from dart comes before this check if the alphabet is UN
  if (alphabet === ALPHABET.PT) {
    polymerType = HELM_POLYMER_TYPE.PEPTIDE;
  } else if (alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA) {
    polymerType = HELM_POLYMER_TYPE.RNA;
  } else {
    grok.shell.warning(
      `Only PT, DNA and RNA alphabets are supported, while the selected column has ${polymerType} alphabet`
    );
    return;
  }

  // work in standard mode, where, as in HELMCoreLibrary:
  // - monomers with polymerType 'PEPTIDE' have monomer type 'backbone'
  // - monomers with polymerType 'RNA' have monomer type 'branch' and 'backbone'
  // - the library provides molfiles in format V2000
  const mode = MODE.STANDARD;

  const monomerSequencesArray: string[][] = getMonomerSequencesArray(macroMolCol);
  const monomersDict = await getMonomersDictFromLib(monomerSequencesArray, monomersLibList, polymerType, alphabet);
  const columnLength = macroMolCol.length;
  const reconstructed: string[] = new Array(columnLength);
  for (let row = 0; row < columnLength; ++row) {
    const monomerSeq = monomerSequencesArray[row];
    reconstructed[row] = monomerSeqToMolfile(monomerSeq, monomersDict, alphabet, polymerType, mode);
    // console.log(reconstructed[row]);
  }

  // exclude name collisions
  const name = 'molfile(' + macroMolCol.name + ')';
  const newColName = df.columns.getUnusedName(name);
  const newCol = DG.Column.fromStrings(newColName, reconstructed);

  newCol.semType = DG.SEMTYPE.MOLECULE;
  newCol.setTag(DG.TAGS.UNITS, DG.UNITS.Molecule.MOLBLOCK);
  df.columns.add(newCol, true);
  await grok.data.detectSemanticTypes(df);
}

/** Get a mapping of peptide symbols to HELM monomer library
 * objects with selected fields.
 */
function getFormattedMonomerLib(
  monomersLibList: any[], polymerType: HELM_POLYMER_TYPE, alphabet: ALPHABET
): Map<string, any> {
  const map = new Map<string, any>();
  monomersLibList.forEach(
    (it) => {
      if (it[HELM_FIELDS.POLYMER_TYPE] === polymerType) {
        if (
          polymerType === HELM_POLYMER_TYPE.RNA &&
          (it[HELM_FIELDS.MONOMER_TYPE] === HELM_MONOMER_TYPE.BRANCH ||
            alphabet === ALPHABET.DNA && it[HELM_FIELDS.SYMBOL] === DEOXYRIBOSE ||
            alphabet === ALPHABET.RNA && it[HELM_FIELDS.SYMBOL] === RIBOSE ||
            it[HELM_FIELDS.SYMBOL] === PHOSPHATE) ||
          polymerType === HELM_POLYMER_TYPE.PEPTIDE &&
          it[HELM_FIELDS.MONOMER_TYPE] !== HELM_MONOMER_TYPE.BRANCH
        ) {
          const monomerObject: { [key: string]: any } = {};
          HELM_CORE_FIELDS.forEach((field) => {
            monomerObject[field] = it[field];
          });
          map.set(it[HELM_FIELDS.SYMBOL], monomerObject);
        }
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
  const separator = macroMolCol.getTag(TAGS.separator);
  const splitterFunc: SplitterFunc = getSplitter(colUnits, separator);

  for (let row = 0; row < columnLength; ++row) {
    const macroMolecule = macroMolCol.get(row);
    // todo: handle the exception case when macroMolecule is null
    result[row] = macroMolecule ? splitterFunc(macroMolecule) : [];
  }
  return result;
}

/* Get a mapping of monomer symbols to MolGraph objects. Notice, the
 * transformation from molfile V2000 to V3000 takes place,
 * with the help of async function call from Chem (RdKit module) */
async function getMonomersDictFromLib(
  monomerSequencesArray: string[][], monomersLibList: any[], polymerType: HELM_POLYMER_TYPE, alphabet: ALPHABET
): Promise<Map<string, MolGraph>> {
  // todo: exception - no gaps, no empty string monomers
  const formattedMonomerLib = getFormattedMonomerLib(monomersLibList, polymerType, alphabet);
  const monomersDict = new Map<string, MolGraph>();

  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');

  // add deoxyribose/ribose and phosphate for nucleotide sequences
  if (polymerType === HELM_POLYMER_TYPE.RNA) {
    const symbols = (alphabet === ALPHABET.RNA) ?
      [RIBOSE, PHOSPHATE] : [DEOXYRIBOSE, PHOSPHATE];
    for (const sym of symbols)
      updateMonomersDict(monomersDict, sym, formattedMonomerLib, moduleRdkit, polymerType);
  }

  for (let row = 0; row < monomerSequencesArray.length; ++row) {
    const monomerSeq: string[] = monomerSequencesArray[row];
    for (const sym of monomerSeq)
      updateMonomersDict(monomersDict, sym, formattedMonomerLib, moduleRdkit, polymerType);
  }
  // console.log(monomersDict);

  return monomersDict;
}

/* Get a mapping of monomer symbols to MolGraph objects from a map whose keys
 * are symbols and values, V3000 molfiles */
function getMonomersDictFromMap(symbolToMolfileV3KMap: Map<string, string>): Map<string, MolGraph> {
  const monomersDict = new Map<string, MolGraph>();
  const mapKeyList = Array.from(symbolToMolfileV3KMap.keys());

  for (const sym of mapKeyList) {
    const molfileV3K = symbolToMolfileV3KMap.get(sym)!; // ! is guaranteed

    const counts = parseAtomAndBondCounts(molfileV3K);
    const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
    const bonds = parseBondBlock(molfileV3K, counts.bondCount);

    // the rNodes are set to 0th and the last atom!
    // this is used because as for now it is unclear how to include the R-groups
    // into the V3K molfile
    const meta = getMonomerMetadata(atoms, bonds);

    const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

    const leftNodeIdx = meta.rNodes[0] - 1;
    const rightNodeIdx = meta.terminalNodes[1] - 1;
    // todo: consider rotation?
    adjustBackboneMonomerGraph(monomerGraph, leftNodeIdx, rightNodeIdx);

    // set shifts
    // todo: wrap as a separate function?
    monomerGraph.meta.backboneShift = getShiftBetweenNodes(monomerGraph, rightNodeIdx, leftNodeIdx);

    removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);

    monomersDict.set(sym, monomerGraph);
  }
  // console.log(monomersDict);

  return monomersDict;
}

/* Adds MolGraph object for 'sym' to the monomers dict when necessary  */
function updateMonomersDict(
  monomersDict: Map<string, MolGraph>, sym: string,
  formattedMonomerLib: Map<string, any>, moduleRdkit: any, polymerType: HELM_POLYMER_TYPE
): void {
  if (!monomersDict.has(sym)) {
    const monomerData: MolGraph | null = getMolGraph(sym, formattedMonomerLib, moduleRdkit, polymerType);
    if (monomerData)
      monomersDict.set(sym, monomerData);
    else
      throw new Error(`Monomer with symbol '${sym}' is absent the monomer library`);
    // todo: handle exception
  }
}

/* Construct the MolGraph object for specified monomerSymbol: the associated
 * graph is adjusted in XY plane and filled with default R-groups */
function getMolGraph(
  monomerSymbol: string, formattedMonomerLib: Map<string, any>,
  moduleRdkit: any, polymerType: HELM_POLYMER_TYPE // todo: specify type for moduleRdkit
): MolGraph | null {
  if (!formattedMonomerLib.has(monomerSymbol)) {
    return null;
  } else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    const capGroups = parseCapGroups(libObject[HELM_FIELDS.RGROUPS]);
    const capGroupIdxMap = parseCapGroupIdxMap(libObject[HELM_FIELDS.MOLFILE]);
    const molfileV3K = convertMolfileToV3K(removeRGroupLines(libObject[HELM_FIELDS.MOLFILE]), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);

    const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
    const bonds = parseBondBlock(molfileV3K, counts.bondCount);
    const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);

    const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
      adjustPeptideMonomerGraph(monomerGraph);
    } else { // nucleotides
      if (monomerSymbol === RIBOSE || monomerSymbol === DEOXYRIBOSE)
        adjustSugarMonomerGraph(monomerGraph);
      else if (monomerSymbol === PHOSPHATE)
        adjustPhosphateMonomerGraph(monomerGraph);
      else
        adjustBaseMonomerGraph(monomerGraph);
    }

    // remove the 'rightmost' chain-extending r-group node in the backbone
    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
      setShifts(monomerGraph, polymerType);
      removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
    } else { // nucleotides
      if (monomerSymbol === RIBOSE || monomerSymbol === DEOXYRIBOSE) {
        // remove R2
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
        // set terminalNode2 (oxygen) as new R2
        monomerGraph.meta.rNodes[1] = monomerGraph.meta.terminalNodes[1];
        setTerminalNodes(monomerGraph.bonds, monomerGraph.meta); // set terminal nodes anew
        setShifts(monomerGraph, polymerType);
        // remove 'new' R2 (oxygen)
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
        // remove R1
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[0]);
        // remove the branching r-group
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[2]);
      } else if (monomerSymbol === PHOSPHATE) {
        monomerGraph.meta.terminalNodes[0] = monomerGraph.meta.rNodes[0];
        shiftCoordinates(
          monomerGraph,
          -monomerGraph.atoms.x[monomerGraph.meta.terminalNodes[0] - 1],
          -monomerGraph.atoms.y[monomerGraph.meta.terminalNodes[0] - 1]
        );
        setShifts(monomerGraph, polymerType);
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
      } else { // nucleobases
        removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[0]);
      }
    }
    removeHydrogen(monomerGraph);

    return monomerGraph;
  }
}

// todo: sdoc
function getMonomerMetadata(atoms: Atoms, bonds: Bonds, capGroups?: string[], capGroupIdxMap?: Map<number, number>
): MonomerMetadata {
  const meta: MonomerMetadata = {
    backboneShift: null,
    branchShift: null,
    terminalNodes: [],
    rNodes: [],
  };

  // corresponds to MODE.STANDARD
  const standardMode = typeof capGroups !== 'undefined' && typeof capGroupIdxMap !== 'undefined';

  if (standardMode) {
    substituteCapGroups(atoms, capGroups!, capGroupIdxMap!);
    setRNodes(capGroupIdxMap!, meta);
  } else { // the case used in SequenceTranslator
    // todo: verify that the monomers are prepared in such a way that this works
    meta.rNodes = [0, atoms.x.length];
  }

  setTerminalNodes(bonds, meta);
  return meta;
}

/* Parse element symbols for R-groups from the HELM monomer library R-groups
 * field  */
export function parseCapGroups(rGroupObjList: any[]): string[] {
  // specifically for HELMCoreLibrary
  // considered only monoatomic rgroups
  // supposing that elements in rGroupObjList are sorted w.r.t. the rgroups idx
  // todo: possible generalizations
  const capGroupsArray: string[] = [];
  for (const obj of rGroupObjList) {
    let capGroup: string = obj[RGROUP_FIELDS.CAP_GROUP_SMILES];

    // in some cases the smiles field is written with uppercase
    if (!capGroup)
      capGroup = obj[RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
    // todo: verify that there are no multi-element cap groups, or consider how to
    // transform them
    capGroup = capGroup.replace(/(\[|\]|\*|:|\d)/g, '');
    if (capGroup.length > 1) // todo: check if such cases are possible, remove if not
      throw new Error('Default cap group has length more than one');
    capGroupsArray.push(capGroup);
  }
  return capGroupsArray;
}

/* Substitute the cap group elements instead of R# */
function substituteCapGroups(
  atoms: Atoms, capGroups: string[], capGroupIdxMap: Map<number, number>
): void {
  for (const [node, capIdx] of capGroupIdxMap)
    atoms.atomTypes[node - 1] = capGroups[capIdx - 1]; // -1 because molfile indexing starts from 1
}

//todo: doc
function setRNodes(capGroupIdxMap: Map<number, number>, meta: MonomerMetadata): void {
  meta.rNodes = Array.from(capGroupIdxMap.keys());
  for (let i = 0; i < meta.rNodes.length; i++) {
    for (const j of [1, 2]) { // 1 and 2 by def. correspond to 'left/rightmost' r-nodes
      // swap the values if necessary, so that the "leftmost" r-node is at 0,
      // and the 'rightmost', at 1
      if (capGroupIdxMap.get(meta.rNodes[i]) === j) {
        const tmp = meta.rNodes[j - 1];
        meta.rNodes[j - 1] = meta.rNodes[i];
        meta.rNodes[i] = tmp;
      }
    }
  }
}

//todo: doc
function setTerminalNodes(bonds: Bonds, meta: MonomerMetadata): void {
  const rNodes = meta.rNodes;
  meta.terminalNodes = new Array<number>(rNodes.length).fill(0);
  const terminalNodes = meta.terminalNodes;
  const atomPairs = bonds.atomPairs;
  let i = 0;
  let j = 0;
  while ((i < atomPairs.length) && j < terminalNodes.length) {
    // rNodes array is sorted so that its 0th and 1st elements (if both
    // present) correspond to the chain extending (i.e. not branching) r-groups
    for (let k = 0; k < terminalNodes.length; ++k) {
      for (let l = 0; l < 2; ++l) {
        if (atomPairs[i][l] === rNodes[k]) {
          terminalNodes[k] = atomPairs[i][(l + 1) % 2];
          if (rNodes.length > 2) {
          }
          ++j;
        }
      }
    }
    ++i;
  }
}

//todo: doc
function setShifts(molGraph: MolGraph, polymerType: HELM_POLYMER_TYPE): void {
  if (molGraph.meta.rNodes.length > 1) {
    molGraph.meta.backboneShift = getShiftBetweenNodes(
      molGraph, molGraph.meta.rNodes[1] - 1,
      molGraph.meta.terminalNodes[0] - 1
    );
  }

  if (polymerType === HELM_POLYMER_TYPE.RNA && molGraph.meta.rNodes.length > 2) {
    molGraph.meta.branchShift = getShiftBetweenNodes(
      molGraph, molGraph.meta.rNodes[2] - 1,
      molGraph.meta.terminalNodes[0] - 1
    );
  }
}

/* Returns the pair [xShift, yShift] for specified node indices */
function getShiftBetweenNodes(
  molGraph: MolGraph, rightNodeIdx: number, leftNodeIdx: number
): number[] {
  return [
    keepPrecision(
      molGraph.atoms.x[rightNodeIdx] -
      molGraph.atoms.x[leftNodeIdx]
    ),
    keepPrecision(
      molGraph.atoms.y[rightNodeIdx] -
      molGraph.atoms.y[leftNodeIdx]
    ),
  ];
}

/* Helper function necessary to build a correct V3000 molfile out of V2000 with
 * specified r-groups*/
function removeRGroupLines(molfileV2K: string): string {
  let begin = molfileV2K.indexOf(V2K_A_LINE, 0);
  if (begin === -1)
    begin = molfileV2K.indexOf(V2K_RGP_LINE);
  const end = molfileV2K.indexOf(V3K_END, begin);
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

/* Parse V3000 bond block and construct the Bonds object */
function parseBondBlock(molfileV3K: string, bondCount: number): Bonds {
  // todo: consider the case when there is no simple leftmost/rightmost choice
  // todo: consider the case when there are multiple consequent M  RGP lines,
  // like in HELMCoreLibrary nucleotides

  const bondTypes: number[] = new Array(bondCount);
  const atomPairs: number[][] = new Array(bondCount);
  const bondConfiguration = new Map<number, number>();
  const kwargs = new Map<number, string>();

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
    bondTypes[i] = parsedValues[0];
    atomPairs[i] = parsedValues.slice(1);

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
    bondTypes: bondTypes,
    atomPairs: atomPairs,
    bondConfiguration: bondConfiguration,
    kwargs: kwargs,
  };
}

/* Constructs mapping of r-group nodes to default capGroups, all numeration starting from 1.
 * According to https://pubs.acs.org/doi/10.1021/ci3001925, R1 and R2 are the chain extending attachment points,
 * while R3 is the branching attachment point. */
function parseCapGroupIdxMap(molfileV2K: string): Map<number, number> {
  const capGroupIdxMap = new Map<number, number>();

  // parse A-lines (RNA)
  let begin = molfileV2K.indexOf(V2K_A_LINE, 0);
  let end = begin;
  while (begin !== -1) {
    // parse the rNode to which the cap group is attached
    end = molfileV2K.indexOf('\n', begin);
    const rNode = parseInt(molfileV2K.substring(begin, end).replace(/^A\s+/, ''));

    // parse the capGroup index
    begin = molfileV2K.indexOf('R', end);
    end = molfileV2K.indexOf('\n', begin);
    const capGroup = parseInt(molfileV2K.substring(begin, end).replace(/^R/, ''));
    capGroupIdxMap.set(rNode, capGroup);

    begin = molfileV2K.indexOf(V2K_A_LINE, end);
  }

  // parse RGP lines (may be more than one in RNA monomers)
  begin = molfileV2K.indexOf(V2K_RGP_LINE, 0);
  end = molfileV2K.indexOf('\n', begin);
  while (begin !== -1) {
    begin += V2K_RGP_SHIFT;
    end = molfileV2K.indexOf('\n', begin);
    const rgpStringParsed = molfileV2K.substring(begin, end)
      .replaceAll(/\s+/g, ' ')
      .split(' ');
    const rgpIndicesArray = rgpStringParsed.map((el) => parseInt(el))
      .slice(1); // slice from 1 because the 1st value is the number of pairs in the line
    for (let i = 0; i < rgpIndicesArray.length; i += 2) {
      // notice: there may be conflicting cap group definitions, like 3-O-Methylribose (2,5 connectivity)
      // (the last monomer in HELMCoreLibrary)
      // there the indices of cap groups are self-contradictory
      // todo: clarify why such situations occur in principle
      if (capGroupIdxMap.has(rgpIndicesArray[i]) && capGroupIdxMap.get(rgpIndicesArray[i]) !== rgpIndicesArray[i + 1])
        throw new Error(`r-group index ${rgpIndicesArray[i]} has already been added with a different value`);
      else
        capGroupIdxMap.set(rgpIndicesArray[i], rgpIndicesArray[i + 1]);
    }

    begin = molfileV2K.indexOf(V2K_RGP_LINE, end);
  }

  return capGroupIdxMap;
}

function parseAtomAndBondCounts(molfileV3K: string): { atomCount: number, bondCount: number } {
  molfileV3K = molfileV3K.replaceAll('\r', ''); // to handle old and new sdf standards

  // parse atom count
  let begin = molfileV3K.indexOf(V3K_BEGIN_COUNTS_LINE) + V3K_COUNTS_SHIFT;
  let end = molfileV3K.indexOf(' ', begin + 1);
  const numOfAtoms = parseInt(molfileV3K.substring(begin, end));

  // parse bond count
  begin = end + 1;
  end = molfileV3K.indexOf(' ', begin + 1);
  const numOfBonds = parseInt(molfileV3K.substring(begin, end));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

/* Parse V3000 atom block and return Atoms object. NOTICE: only atomTypes, x, y
 * and kwargs fields are set in the return value, with other fields dummy */
function parseAtomBlock(molfileV3K: string, atomCount: number): Atoms {
  const atomTypes: string[] = new Array(atomCount);
  const x: number[] = new Array(atomCount);
  const y: number[] = new Array(atomCount);
  const kwargs: string[] = new Array(atomCount);

  let begin = molfileV3K.indexOf(V3K_BEGIN_ATOM_BLOCK); // V3000 atoms block
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;

  for (let i = 0; i < atomCount; i++) {
    begin = molfileV3K.indexOf(V3K_BEGIN_DATA_LINE, begin) + V3K_IDX_SHIFT;
    end = molfileV3K.indexOf(' ', begin); // skip the idx row

    // parse atom type
    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin);
    atomTypes[i] = molfileV3K.substring(begin, end);

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
    atomTypes: atomTypes,
    x: x,
    y: y,
    kwargs: kwargs,
  };
}

/* Remove hydrogen nodes */
function removeHydrogen(monomerGraph: MolGraph): void {
  let i = 0;
  while (i < monomerGraph.atoms.atomTypes.length) {
    if (monomerGraph.atoms.atomTypes[i] === HYDROGEN) {
      removeNodeAndBonds(monomerGraph, i + 1); // i + 1 because molfile node indexing starts from 1
      --i;
      // monomerGraph.atoms.atomTypes[i] = 'Li';
    }
    ++i;
  }
}

/* Remove node 'removedNode' and the associated bonds. Notice, numeration of
 * nodes in molfiles starts from 1, not 0 */
function removeNodeAndBonds(monomerGraph: MolGraph, removedNode?: number): void {
  if (typeof removedNode !== 'undefined') {
    const removedNodeIdx = removedNode - 1;
    const atoms = monomerGraph.atoms;
    const bonds = monomerGraph.bonds;
    const meta = monomerGraph.meta;

    // remove the node from atoms
    atoms.atomTypes.splice(removedNodeIdx, 1);
    atoms.x.splice(removedNodeIdx, 1);
    atoms.y.splice(removedNodeIdx, 1);
    atoms.kwargs.splice(removedNodeIdx, 1);

    // update the values of terminal and r-group nodes if necessary
    for (let i = 0; i < meta.terminalNodes.length; ++i) {
      if (meta.terminalNodes[i] > removedNode)
        --meta.terminalNodes[i];
      else if (meta.terminalNodes[i] === removedNode)
        meta.terminalNodes[i] = -1; // sentinel to mark the value as removed
    }
    for (let i = 0; i < meta.rNodes.length; ++i) {
      if (meta.rNodes[i] > removedNode)
        --meta.rNodes[i];
      else if (meta.rNodes[i] === removedNode)
        meta.rNodes[i] = -1; // sentinel to mark the value as removed
    }

    // update indices of atoms in bonds
    let i = 0;
    while (i < bonds.atomPairs.length) {
      const firstAtom = bonds.atomPairs[i][0];
      const secondAtom = bonds.atomPairs[i][1];
      if (firstAtom === removedNode || secondAtom === removedNode) {
        bonds.atomPairs.splice(i, 1);
        bonds.bondTypes.splice(i, 1);
        if (bonds.bondConfiguration.has(i))
          bonds.bondConfiguration.delete(i);
        if (bonds.kwargs.has(i))
          bonds.kwargs.delete(i);
        --i;
      } else {
        bonds.atomPairs[i][0] = (firstAtom > removedNode) ? firstAtom - 1 : firstAtom;
        bonds.atomPairs[i][1] = (secondAtom > removedNode) ? secondAtom - 1 : secondAtom;
      }
      ++i;
    }

    // update bondConfiguration and kwargs keys
    let keys = Array.from(bonds.bondConfiguration.keys());
    keys.forEach((key) => {
      if (bonds.bondConfiguration.has(key) && key > removedNodeIdx) {
        const value = bonds.bondConfiguration.get(key)!;
        bonds.bondConfiguration.delete(key);
        bonds.bondConfiguration.set(key - 1, value);
      }
    });
    keys = Array.from(bonds.kwargs.keys());
    keys.forEach((key) => {
      if (bonds.kwargs.has(key) && key > removedNodeIdx) {
        const value = bonds.kwargs.get(key)!;
        bonds.kwargs.delete(key);
        bonds.kwargs.set(key - 1, value);
      }
    });
  }
}

// todo: rewrite description
/* Adjust the (peptide) monomer graph so that it has standard form  */
function adjustPeptideMonomerGraph(monomer: MolGraph): void {
  const nodeOneIdx = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  const nodeTwoIdx = monomer.meta.rNodes[0] - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // angle is measured between OY and the rotated node
  const angle = findAngleWithOY(x[nodeTwoIdx], y[nodeTwoIdx]);

  // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  rotateCenteredGraph(monomer.atoms, -angle);

  if (x[monomer.meta.rNodes[1] - 1] < 0)
    flipMonomerAroundOY(monomer);

  const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);

  // flip carboxyl and R if necessary
  flipCarboxylAndRadical(monomer, doubleBondedOxygen);

  // flip hydroxyl group with double-bound O inside carboxyl group if necessary
  flipHydroxilGroup(monomer, doubleBondedOxygen);
}

function adjustPhosphateMonomerGraph(monomer: MolGraph): void {
  const nodeOneIdx = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  // const nodeTwoIdx = monomer.meta.rNodes[0] - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // // angle is measured between OY and the rotated node
  // const angle = findAngleWithOY(x[nodeTwoIdx], y[nodeTwoIdx]);

  // // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  // rotateCenteredGraph(monomer.atoms, -angle);

  // if (x[monomer.meta.rNodes[1] - 1] < 0)
  //   flipMonomerAroundOY(monomer);

  // const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);

  // // flip carboxyl and R if necessary
  // flipCarboxylAndRadical(monomer, doubleBondedOxygen);

  // // flip hydroxyl group with double-bound O inside carboxyl group if necessary
  // flipHydroxilGroup(monomer, doubleBondedOxygen);
}

/* Adjust a backbone graph so that nodeOne is at origin and nodeTwo is at OX.
 * Notice: node indexing in molfiles starts from 1 */
function adjustBackboneMonomerGraph(
  monomer: MolGraph, nodeOneIdx: number, nodeTwoIdx: number
): void {
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // angle is measured between OX and the rotated node
  const angle = findAngleWithOX(x[nodeTwoIdx], y[nodeTwoIdx]);

  // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OX
  rotateCenteredGraph(monomer.atoms, -angle);
}

function adjustSugarMonomerGraph(monomer: MolGraph): void {
  const nodeOneIdx = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  // const nodeTwoIdx = monomer.meta.rNodes[0] - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // // angle is measured between OY and the rotated node
  // const angle = findAngleWithOY(x[nodeTwoIdx], y[nodeTwoIdx]);

  // // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  // rotateCenteredGraph(monomer.atoms, -angle);

  // if (x[monomer.meta.rNodes[1] - 1] < 0)
  //   flipMonomerAroundOY(monomer);

  // const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);

  // // flip carboxyl and R if necessary
  // flipCarboxylAndRadical(monomer, doubleBondedOxygen);

  // // flip hydroxyl group with double-bound O inside carboxyl group if necessary
  // flipHydroxilGroup(monomer, doubleBondedOxygen);
}

function adjustBaseMonomerGraph(monomer: MolGraph): void {
  const nodeOneIdx = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  // const nodeTwoIdx = monomer.meta.rNodes[0] - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[nodeOneIdx], -y[nodeOneIdx]);

  // // angle is measured between OY and the rotated node
  // const angle = findAngleWithOY(x[nodeTwoIdx], y[nodeTwoIdx]);

  // // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  // rotateCenteredGraph(monomer.atoms, -angle);

  // if (x[monomer.meta.rNodes[1] - 1] < 0)
  //   flipMonomerAroundOY(monomer);

  // const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);

  // // flip carboxyl and R if necessary
  // flipCarboxylAndRadical(monomer, doubleBondedOxygen);

  // // flip hydroxyl group with double-bound O inside carboxyl group if necessary
  // flipHydroxilGroup(monomer, doubleBondedOxygen);
}

/* Flip carboxyl group with the radical in a peptide monomer in case the
 * carboxyl group is in the lower half-plane */
function flipCarboxylAndRadical(monomer: MolGraph, doubleBondedOxygen: number): void {
  // verify that the carboxyl group is in the lower half-plane
  if (monomer.atoms.y[monomer.meta.rNodes[1] - 1] < 0 &&
    monomer.atoms.y[doubleBondedOxygen - 1] < 0) {
    flipMonomerAroundOX(monomer);

    rotateCenteredGraph(monomer.atoms,
      -findAngleWithOX(
        monomer.atoms.x[monomer.meta.terminalNodes[1] - 1],
        monomer.atoms.y[monomer.meta.terminalNodes[1] - 1]
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
    angle = x > 0 ? -Math.PI / 2 : Math.PI / 2;
  } else {
    const tan = y / x;
    const atan = Math.atan(tan);
    angle = (x < 0) ? Math.PI / 2 + atan : -Math.PI / 2 + atan;
  }
  return angle;
}

/* Finds angle between OX and the ray joining origin with (x, y) */
function findAngleWithOX(x: number, y: number): number {
  return findAngleWithOY(x, y) + Math.PI / 2;
}

/*  Rotate the graph around the origin by 'angle' */
function rotateCenteredGraph(atoms: Atoms, angle: number): void {
  if (angle !== 0) {
    const x = atoms.x;
    const y = atoms.y;

    const cos = Math.cos(angle);
    const sin = Math.sin(angle);

    for (let i = 0; i < x.length; ++i) {
      const tmp = x[i];
      x[i] = keepPrecision(tmp * cos - y[i] * sin);
      y[i] = keepPrecision(tmp * sin + y[i] * cos);
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

/* Flips double-bonded 'O' in carbonyl group with 'OH' in order for the monomers
 * to have standard representation simplifying their concatenation. The
 * monomer must already be adjusted with adjustPeptideMonomerGraph in order for this function to be implemented  */
function flipHydroxilGroup(monomer: MolGraph, doubleBondedOxygen: number): void {
  const x = monomer.atoms.x;
  // -1 below because indexing of nodes in molfiles starts from 1, unlike arrays
  if (x[monomer.meta.rNodes[1] - 1] > x[doubleBondedOxygen - 1])
    swapNodes(monomer, doubleBondedOxygen, monomer.meta.rNodes[1]);
}

/* Determine the number of node (starting from 1) corresponding to the
 * double-bonded oxygen of the carbonyl group  */
function findDoubleBondedCarbonylOxygen(monomer: MolGraph): number {
  const bondsMap = constructBondsMap(monomer);
  let doubleBondedOxygen = 0;
  let i = 0;
  // iterate over the nodes bonded to the carbon and find the double one
  while (doubleBondedOxygen === 0) {
    const node = bondsMap.get(monomer.meta.terminalNodes[1])![i];
    if (monomer.atoms.atomTypes[node - 1] === OXYGEN && node !== monomer.meta.rNodes[1])
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

// todo: doc
function constructBondsMap(monomer: MolGraph): Map<number, Array<number>> {
  const map = new Map<number, Array<number>>();
  for (const atomPairs of monomer.bonds.atomPairs) {
    for (let i = 0; i < 2; i++) {
      const key = atomPairs[i];
      const value = atomPairs[(i + 1) % 2];
      if (map.has(key))
        map.get(key)?.push(value);
      else
        map.set(key, new Array<number>(1).fill(value));
    }
  }
  return map;
}

/* Shift molGraph in the XOY plane  */
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
function monomerSeqToMolfile(
  monomerSeq: string[], monomersDict: Map<string, MolGraph>,
  alphabet: ALPHABET, polymerType: HELM_POLYMER_TYPE, mode: MODE
): string {
  // todo: handle the case when the polymer is empty
  if (monomerSeq.length === 0)
    throw new Error('monomerSeq is empty');

  // define atom and bond counts, taking into account the bond type
  const getAtomAndBondCounts = (mode === MODE.STANDARD) ?
    getResultingAtomBondCounts : getResultingAtomBondCountsST;
  const {atomCount, bondCount} = getAtomAndBondCounts(monomerSeq, monomersDict, alphabet, polymerType, mode);

  // create arrays to store lines of the resulting molfile
  const molfileAtomBlock = new Array<string>(atomCount);
  const molfileBondBlock = new Array<string>(bondCount);

  let addMonomerToMolblock; // todo: types?

  let sugar = null;
  let phosphate = null;

  if (mode === MODE.STANDARD) {
    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
      addMonomerToMolblock = addAminoAcidToMolblock;
    } else { // nucleotides
      addMonomerToMolblock = addNucleotideToMolblock;
      sugar = (alphabet === ALPHABET.DNA) ? monomersDict.get(DEOXYRIBOSE) : monomersDict.get(RIBOSE);
      phosphate = monomersDict.get(PHOSPHATE);
    }
  } else {
    addMonomerToMolblock = addMonomerToMolblockST;
  }

  const v: LoopVariables = {
    i: 0,
    nodeShift: 0,
    bondShift: 0,
    backbonePositionShift: new Array<number>(2).fill(0),
    branchPositionShift: new Array<number>(2).fill(0),
    backboneAttachNode: 0,
    branchAttachNode: 0,
    flipFactor: 1,
  };

  const C: LoopConstants = {
    sugar: sugar!,
    phosphate: phosphate!,
    seqLength: monomerSeq.length,
    atomCount: atomCount,
    bondCount: bondCount,
  };

  for (v.i = 0; v.i < C.seqLength; ++v.i) {
    const monomer = monomersDict.get(monomerSeq[v.i])!;
    addMonomerToMolblock(monomer, molfileAtomBlock, molfileBondBlock, v, C);
  }

  capResultingMolblock(molfileAtomBlock, molfileBondBlock, v, C, mode);

  const molfileCountsLine = V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + V3K_COUNTS_LINE_ENDING;

  // todo: optimize concatenation using Alexander's hint
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

/* Cap the resulting (after sewing up all the monomers) molfile with 'O' */
function capResultingMolblock(
  molfileAtomBlock: string[], molfileBondBlock: string[],
  v: LoopVariables, C: LoopConstants, mode: MODE
): void {
  if (mode === MODE.STANDARD) {
    // add terminal oxygen
    const atomIdx = v.nodeShift + 1;
    molfileAtomBlock[C.atomCount] = V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
      OXYGEN + ' ' + keepPrecision(v.backbonePositionShift[0]) + ' ' +
      v.flipFactor * keepPrecision(v.backbonePositionShift[1]) + ' ' + '0.000000 0' + '\n';

    // add terminal bond
    const firstAtom = v.backboneAttachNode;
    const secondAtom = atomIdx;
    molfileBondBlock[C.bondCount] = V3K_BEGIN_DATA_LINE + v.bondShift + ' ' +
      1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
  }
}

function addAminoAcidToMolblock(monomer: MolGraph, molfileAtomBlock: string[],
  molfileBondBlock: string[], v: LoopVariables, C: LoopConstants
): void {
  v.flipFactor = (-1) ** (v.i % 2); // to flip every even monomer over OX
  addBackboneMonomerToMolblock(monomer, molfileAtomBlock, molfileBondBlock, v, C);
}

function addMonomerToMolblockST(monomer: MolGraph, molfileAtomBlock: string[],
  molfileBondBlock: string[], v: LoopVariables, C: LoopConstants
): void {
  addBranchMonomerToMolblock(monomer, molfileAtomBlock, molfileBondBlock, v, C);
}

function addBackboneMonomerToMolblock(
  monomer: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables, C: LoopConstants
): void {
  // todo: remove these comments to the docstrings of the corr. functions
  // construnct the lines of V3K molfile atom block
  fillAtomLines(monomer, molfileAtomBlock, v);

  // construct the lines of V3K molfile bond block
  fillBondLines(monomer, molfileBondBlock, v);

  // peptide bond
  fillChainExtendingBond(monomer, molfileBondBlock, v);

  // update branch variables if necessary
  if (monomer.meta.branchShift !== null && monomer.meta.terminalNodes.length > 2)
    updateBranchVariables(monomer, v);

  // update loop variables
  updateChainExtendingVariables(monomer, v, C);
}

/* Defined for 'standard mode' RNA monomers, i.e. of 'branch' type, as in
 * HELMCoreLibrary. Will not work in SEQ_TRAN mode */
function addNucleotideToMolblock(
  nucleobase: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables, C: LoopConstants
): void {
  // construnct the lines of V3K molfile atom block corresponding to phosphate
  // and sugar
  for (const monomer of [C.phosphate, C.sugar])
    addBackboneMonomerToMolblock(monomer!, molfileAtomBlock, molfileBondBlock, v, C);

  addBranchMonomerToMolblock(nucleobase, molfileAtomBlock, molfileBondBlock, v, C);
}

function addBranchMonomerToMolblock(
  monomer: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables, C: LoopConstants
): void {
  fillBranchAtomLines(monomer, molfileAtomBlock, v);
  fillBondLines(monomer, molfileBondBlock, v);
  fillBackboneToBranchBond(monomer, molfileBondBlock, v);

  // C-N bond
  const bondIdx = v.bondShift;
  const firstAtom = v.branchAttachNode;
  const secondAtom = monomer.meta.terminalNodes[0] + v.nodeShift;
  molfileBondBlock[bondIdx - 1] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
    1 + ' ' + firstAtom + ' ' + secondAtom + '\n';

  // update loop variables
  v.bondShift += monomer.bonds.atomPairs.length + 1;
  v.nodeShift += monomer.atoms.atomTypes.length;
}

function updateChainExtendingVariables(monomer: MolGraph, v: LoopVariables, C: LoopConstants): void {
  v.backboneAttachNode = v.nodeShift + monomer.meta.terminalNodes[1];
  v.bondShift += monomer.bonds.atomPairs.length + 1;

  v.nodeShift += monomer.atoms.atomTypes.length;
  v.backbonePositionShift[0] += monomer.meta.backboneShift![0]; // todo: non-null check
  v.backbonePositionShift[1] += v.flipFactor * monomer.meta.backboneShift![1];
}

function updateBranchVariables(monomer: MolGraph, v: LoopVariables) {
  v.branchAttachNode = v.nodeShift + monomer.meta.terminalNodes[2];
  for (let i = 0; i < 2; ++i)
    v.branchPositionShift[i] = v.backbonePositionShift[i] + monomer.meta.branchShift![i];
}

function fillAtomLines(monomer: MolGraph, molfileAtomBlock: string[], v: LoopVariables): void {
  for (let j = 0; j < monomer.atoms.atomTypes.length; ++j) {
    const atomIdx = v.nodeShift + j + 1;
    molfileAtomBlock[v.nodeShift + j] = V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
      monomer.atoms.atomTypes[j] + ' ' +
      keepPrecision(v.backbonePositionShift[0] + monomer.atoms.x[j]) + ' ' +
      keepPrecision(v.backbonePositionShift[1] + v.flipFactor * monomer.atoms.y[j]) +
      ' ' + monomer.atoms.kwargs[j];
  }
}

// todo: remove as quickfix
function fillBranchAtomLines(monomer: MolGraph, molfileAtomBlock: string[], v: LoopVariables): void {
  for (let j = 0; j < monomer.atoms.atomTypes.length; ++j) {
    const atomIdx = v.nodeShift + j + 1;
    molfileAtomBlock[v.nodeShift + j] = V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
      monomer.atoms.atomTypes[j] + ' ' +
      keepPrecision(v.branchPositionShift[0] + monomer.atoms.x[j]) + ' ' +
      keepPrecision(v.branchPositionShift[1] + v.flipFactor * monomer.atoms.y[j]) +
      ' ' + monomer.atoms.kwargs[j];
  }
}

function fillBondLines(monomer: MolGraph, molfileBondBlock: string[], v: LoopVariables): void {
  // construct the lines of V3K molfile bond block
  for (let j = 0; j < monomer.bonds.atomPairs.length; ++j) {
    const bondIdx = v.bondShift + j + 1;
    const firstAtom = monomer.bonds.atomPairs[j][0] + v.nodeShift;
    const secondAtom = monomer.bonds.atomPairs[j][1] + v.nodeShift;
    let bondCfg = '';
    if (monomer.bonds.bondConfiguration.has(j)) {
      // flip orientation when necessary
      let orientation = monomer.bonds.bondConfiguration.get(j);
      if (v.flipFactor < 0)
        orientation = (orientation === 1) ? 3 : 1;
      bondCfg = ' CFG=' + orientation;
    }
    const kwargs = monomer.bonds.kwargs.has(j) ?
      ' ' + monomer.bonds.kwargs.get(j) : '';
    molfileBondBlock[v.bondShift + j] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
      monomer.bonds.bondTypes[j] + ' ' +
      firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
  }
}

function fillChainExtendingBond(monomer: MolGraph, molfileBondBlock: string[], v: LoopVariables): void {
  if (v.backboneAttachNode !== 0) {
    const bondIdx = v.bondShift;
    const firstAtom = v.backboneAttachNode;
    const secondAtom = monomer.meta.terminalNodes[0] + v.nodeShift;
    molfileBondBlock[v.bondShift - 1] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
      1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
  }
}

// todo: remove
function fillBackboneToBranchBond(branchMonomer: MolGraph, molfileBondBlock: string[], v: LoopVariables): void {
  const bondIdx = v.bondShift;
  const firstAtom = v.branchAttachNode;
  const secondAtom = branchMonomer.meta.terminalNodes[0] + v.nodeShift;
  molfileBondBlock[bondIdx - 1] = V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
    1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
}

/* Compute the atom/bond counts for the resulting molfile, depending on the
 * type of polymer (peptide/nucleotide) */
function getResultingAtomBondCounts(
  monomerSeq: string[], monomersDict: Map<string, MolGraph>,
  alphabet: ALPHABET, polymerType: HELM_POLYMER_TYPE,
  mode: MODE
): { atomCount: number, bondCount: number } {
  let atomCount = 0;
  let bondCount = 0;

  // sum up all the atoms/nodes provided by the sequence
  for (const monomerSymbol of monomerSeq) {
    const monomer = monomersDict.get(monomerSymbol)!;
    atomCount += monomer.atoms.x.length;
    bondCount += monomer.bonds.bondTypes.length;
  }

  // add extra values depending on the polymer type
  if (mode === MODE.STANDARD) {
    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
      // add the rightmost/terminating cap group 'OH' (i.e. 'O')
      atomCount += 1;
      // add chain-extending bonds (C-NH per each monomer pair and terminal C-OH)
      bondCount += monomerSeq.length;
    } else { // nucleotides
      const sugar = (alphabet === ALPHABET.DNA) ?
        monomersDict.get(DEOXYRIBOSE)! : monomersDict.get(RIBOSE)!;
      const phosphate = monomersDict.get(PHOSPHATE)!;

      // add phosphate and sugar per each nucleobase symbol
      atomCount += monomerSeq.length * (phosphate.atoms.x.length + sugar.atoms.x.length);
      // add the leftmost cap group 'OH' (i.e. 'O') to the first phosphate
      atomCount += 1;

      // add bonds from phosphate and sugar
      bondCount += monomerSeq.length * (phosphate.bonds.bondTypes.length + sugar.bonds.bondTypes.length);

      // add chain-extending and branch bonds (O-P, C-O and C-N per each nucleotide)
      bondCount += monomerSeq.length * 3;
    }
  } else {
    // todo: fill for SequenceTranslator
  }

  return {atomCount, bondCount};
}

/* Keep precision upon floating point operations over atom coordinates */
function keepPrecision(x: number) {
  return Math.round(PRECISION_FACTOR * x) / PRECISION_FACTOR;
}

function convertMolGraphToMolfileV3K(molGraph: MolGraph): string {
  // counts line
  const atomType = molGraph.atoms.atomTypes;
  const x = molGraph.atoms.x;
  const y = molGraph.atoms.y;
  const atomKwargs = molGraph.atoms.kwargs;
  const bondType = molGraph.bonds.bondTypes;
  const atomPair = molGraph.bonds.atomPairs;
  const bondKwargs = molGraph.bonds.kwargs;
  const bondConfig = molGraph.bonds.bondConfiguration;
  const atomCount = atomType.length;
  const bondCount = molGraph.bonds.bondTypes.length;

  // todo rewrite using constants
  const molfileCountsLine = V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + V3K_COUNTS_LINE_ENDING;

  // atom block
  let molfileAtomBlock = '';
  for (let i = 0; i < atomCount; ++i) {
    const atomIdx = i + 1;
    const coordinate = [x[i].toString(), y[i].toString()];

    // format coordinates so that they have 6 digits after decimal point
    // for (let k = 0; k < 2; ++k) {
    //   const formatted = coordinate[k].toString().split('.');
    //   if (formatted.length === 1)
    //     formatted.push('0');
    //   formatted[1] = formatted[1].padEnd(V3K_ATOM_COORDINATE_PRECISION, '0');
    //   coordinate[k] = formatted.join('.');
    // }

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
    const kwargs = bondKwargs.has(i) ? ' ' + bondKwargs.get(i) : '';
    const bondCfg = bondConfig.has(i) ? ' CFG=' + bondConfig.get(i) : '';
    const bondLine = V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
      firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
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
  const resultingMolfile = molfileParts.join('');
  // console.log(resultingMolfile);

  return resultingMolfile;
}

export async function getSymbolToCappedMolfileMap(monomersLibList: any[]): Promise<Map<string, string> | undefined> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }

  const symbolToCappedMolfileMap = new Map<string, string>();
  const moduleRdkit = await grok.functions.call('Chem:getRdKitModule');

  for (const monomerLibObject of monomersLibList) {
    const monomerSymbol = monomerLibObject[HELM_FIELDS.SYMBOL];
    const capGroups = parseCapGroups(monomerLibObject[HELM_FIELDS.RGROUPS]);
    const capGroupIdxMap = parseCapGroupIdxMap(monomerLibObject[HELM_FIELDS.MOLFILE]);

    const molfileV3K = convertMolfileToV3K(removeRGroupLines(monomerLibObject[HELM_FIELDS.MOLFILE]), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);

    const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
    const bonds = parseBondBlock(molfileV3K, counts.bondCount);
    const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);

    const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

    removeHydrogen(monomerGraph);

    const molfile = convertMolGraphToMolfileV3K(monomerGraph);
    symbolToCappedMolfileMap.set(monomerSymbol, molfile);
  }
  return symbolToCappedMolfileMap;
}

/* Get the V3K molfile corresponding to the capped Monomer (default cap groups)  */
export function capPeptideMonomer(monomer: Monomer): string {
  const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'getRdKitModule'});
  const moduleRdkit = funcList[0].apply();

  const capGroups = parseCapGroups(monomer[HELM_FIELDS.RGROUPS]);
  const capGroupIdxMap = parseCapGroupIdxMap(monomer[HELM_FIELDS.MOLFILE]);
  const molfileV3K = convertMolfileToV3K(removeRGroupLines(monomer[HELM_FIELDS.MOLFILE]), moduleRdkit);
  const counts = parseAtomAndBondCounts(molfileV3K);

  const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
  const bonds = parseBondBlock(molfileV3K, counts.bondCount);
  const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);

  const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

  adjustPeptideMonomerGraph(monomerGraph);

  const molfile = convertMolGraphToMolfileV3K(monomerGraph);
  return molfile;
}

///////////////  Sequence translator /////////////////

/** Currently the ST has peculiar types of monomers, which do not fit the HELM
 * Core library types (in particular, the RNA monomers are backbones only,
 * and presented in Molfile V3K format).
 * TODO: integrate this part with the above functionality
 * Custom _toAtomicLevel version for SequenceTranslator
 */
export function sequenceToMolFileST(
  monomerSeq: string[], // sequence of values of 'symbol' field for monomers
  symbolToMolfileV3KObj: { [symbol: string]: string } // mapping of symbol to molfile V3000
): string | null {
  // work in SEQ_TRAN mode, where:
  // - monomers with polymerType 'RNA' have monomer type 'backbone'
  // - the library provides molfiles in format V3000
  const mode = MODE.SEQ_TRAN;
  const alphabet = ALPHABET.PT; // dummy value! todo: make the argument optional
  const polymerType = HELM_POLYMER_TYPE.RNA; // dummy value! todo: make the argument optional

  // todo: consider refactoring from obj to map in monomer-worls
  const symbolToMolfileV3KMap = new Map<string, string>();
  for (const sym in symbolToMolfileV3KObj)
    symbolToMolfileV3KMap.set(sym, symbolToMolfileV3KObj[sym]);

  const monomersDict = getMonomersDictFromMap(symbolToMolfileV3KMap);
  const result = monomerSeqToMolfile(monomerSeq, monomersDict, alphabet, polymerType, mode);
  // console.log(reconstructed[row]);

  return result;
}

/* Compute the atom/bond counts for the resulting molfile, depending on the
 * type of polymer (peptide/nucleotide) */
function getResultingAtomBondCountsST(
  monomerSeq: string[], monomersDict: Map<string, MolGraph>,
  alphabet: ALPHABET, polymerType: HELM_POLYMER_TYPE
): { atomCount: number, bondCount: number } {
  let atomCount = 0;
  let bondCount = 0;

  // sum up all the atoms/nodes provided by the sequence
  for (const monomerSymbol of monomerSeq) {
    const monomer = monomersDict.get(monomerSymbol)!;
    atomCount += monomer.atoms.x.length;
    bondCount += monomer.bonds.bondTypes.length;
  }

  // add extra values depending on the polymer type
  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
    // add the rightmost/terminating cap group 'OH' (i.e. 'O')
    atomCount += 1;
    // add chain-extending bonds (C-NH per each monomer pair and terminal C-OH)
    bondCount += monomerSeq.length;
  } else { // nucleotides
    const sugar = (alphabet === ALPHABET.DNA) ?
      monomersDict.get(DEOXYRIBOSE)! : monomersDict.get(RIBOSE)!;
    const phosphate = monomersDict.get(PHOSPHATE)!;

    // add phosphate and sugar per each nucleobase symbol
    atomCount += monomerSeq.length * (phosphate.atoms.x.length + sugar.atoms.x.length);
    // add the leftmost cap group 'OH' (i.e. 'O') to the first phosphate
    atomCount += 1;

    // add bonds from phosphate and sugar
    bondCount += monomerSeq.length * (phosphate.bonds.bondTypes.length + sugar.bonds.bondTypes.length);

    // add chain-extending and branch bonds (O-P, C-O and C-N per each nucleotide)
    bondCount += monomerSeq.length * 3;
  }

  return {atomCount, bondCount};
}
