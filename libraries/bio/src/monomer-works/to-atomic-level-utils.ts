/* eslint-disable max-len */
import {monomerWorksConsts as C} from './consts';
import {getMolGraph, LoopConstants, LoopVariables, MolfileWithMap, MolGraph, MonomerMap, MonomerMolGraphMap,
  NucleotideRole} from './types';
import {HELM_FIELDS, HELM_CORE_FIELDS, HELM_POLYMER_TYPE, HELM_MONOMER_TYPE,} from '../utils/const';
import {ALPHABET, GAP_SYMBOL} from '../utils/macromolecule/consts';
import {IMonomerLibBase, Monomer} from '../types/monomer-library';
import {ISeqMonomer, PolymerType} from '../helm/types';
import {helmTypeToPolymerType} from './monomer-works';


/** Get a mapping of peptide symbols to HELM monomer library objects with selected fields.
 * @param {IMonomerLib} monomerLib - Monomer library
 * @param {HELM_POLYMER_TYPE} polymerType - Polymer type
 * @param {ALPHABET} alphabet - Alphabet of the column
 * @return {Map<string, any>} - Mapping of peptide symbols to HELM monomer library objects with selected fields*/
export function getFormattedMonomerLib(
  monomerLib: IMonomerLibBase, polymerType: PolymerType, alphabet: ALPHABET
): Map<string, any> {
  const map = new Map<string, any>();
  for (const monomerSymbol of monomerLib.getMonomerSymbolsByType(polymerType)) {
    const it: Monomer = monomerLib.getMonomer(polymerType, monomerSymbol)!;
    if (
      // RNA: all branch monomers (bases) and all backbone monomers (sugars +
      // phosphates, including modified ones). Modifications are looked up by
      // symbol at assembly time, so they MUST be present in the formatted lib.
      polymerType === HELM_POLYMER_TYPE.RNA &&
      (it[HELM_FIELDS.MONOMER_TYPE] === HELM_MONOMER_TYPE.BRANCH ||
        it[HELM_FIELDS.MONOMER_TYPE] === HELM_MONOMER_TYPE.BACKBONE) ||
      polymerType === HELM_POLYMER_TYPE.PEPTIDE &&
      it[HELM_FIELDS.MONOMER_TYPE] !== HELM_MONOMER_TYPE.BRANCH
    ) {
      const monomerObject: { [key: string]: any } = {};
      HELM_CORE_FIELDS.forEach((field) => {
        //@ts-ignore
        monomerObject[field] = it[field];
      });
      map.set(monomerSymbol, monomerObject);
    }
  }
  return map;
}

/** Translate a sequence of monomer symbols into Molfile V3000
 * @param {ISeqMonomer[]} monomerSeq - Sequence of monomer symbols (canonical)
 * @param {MonomerMolGraphMap} monomersDict - Mapping of monomer symbols to MolGraph objects
 * @param {ALPHABET} alphabet - Alphabet of the column
 * @param {PolymerType} polymerType - Polymer type
 * @param {Array} roles - Optional per-position NucleotideRole tags. When set, RNA assembly
 *   uses per-position sugars/phosphates from monomerSeq directly (HELM triples mode).
 * @return {MolfileWithMap} - Molfile V3000 + per-position monomer index map */
export function monomerSeqToMolfile(
  monomerSeq: ISeqMonomer[], monomersDict: MonomerMolGraphMap,
  alphabet: ALPHABET, polymerType: PolymerType,
  roles?: NucleotideRole[]
): MolfileWithMap {
  if (monomerSeq.length === 0) {
    // throw new Error('monomerSeq is empty');
    return MolfileWithMap.createEmpty();
  }

  // Triples mode is on only when (a) caller flagged the row with roles, and
  // (b) the row's length is a positive multiple of 3 (sugar/base/phosphate triples).
  // in edge cases, terminal phosphate may be missing, so we also allow length % 3 === 2, but only if roles are present to disambiguate the monomers.
  const triplesMode = polymerType === HELM_POLYMER_TYPE.RNA && !!roles &&
    roles.length === monomerSeq.length && (monomerSeq.length % 3 === 0 || monomerSeq.length % 3 === 2);

  // define atom and bond counts, taking into account the bond type
  const {atomCount, bondCount, needsCapping} =
    getResultingAtomBondCounts(monomerSeq, monomersDict, alphabet, polymerType, triplesMode);

  // create arrays to store lines of the resulting molfile
  const molfileAtomBlock = new Array<string>(atomCount);
  const molfileBondBlock = new Array<string>(bondCount);

  let addMonomerToMolblock: (monomer: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables, LC: LoopConstants) => void;

  let sugar = null;
  let phosphate = null;

  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE)
    addMonomerToMolblock = addAminoAcidToMolblock;
  else { // nucleotides
    addMonomerToMolblock = addNucleotideToMolblock;
    // Default sugar/phosphate are only consulted in bases-only mode. In
    // triples mode, every nucleotide carries its own.
    if (!triplesMode) {
      sugar = (alphabet === ALPHABET.DNA) ? getMolGraph(monomersDict, C.DEOXYRIBOSE) : getMolGraph(monomersDict, C.RIBOSE);
      phosphate = getMolGraph(monomersDict, C.PHOSPHATE);
    }
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

  const LC: LoopConstants = {
    sugar: sugar!,
    phosphate: phosphate!,
    // In triples mode, the "logical" sequence length is the nucleotide count.
    seqLength: triplesMode ? Math.ceil(monomerSeq.length / 3) : monomerSeq.length,
    atomCount: atomCount,
    bondCount: bondCount,
  };

  const monomers: MonomerMap = new MonomerMap();
  const steabsCollection: number [] = [];
  let nAtoms = 0;
  let lastMonomerCappingAtom: string | undefined = undefined;

  if (triplesMode) {
    runTriplesAssembly(
      monomerSeq, monomersDict, molfileAtomBlock, molfileBondBlock,
      v, LC, monomers, steabsCollection,
      (a) => { nAtoms += a; }, () => nAtoms);
  } else {
    for (v.i = 0; v.i < LC.seqLength; ++v.i) {
      const seqMonomer = monomerSeq[v.i];
      if (seqMonomer.symbol === GAP_SYMBOL) continue;
      const monomer = getMolGraph(monomersDict, {symbol: seqMonomer.symbol, polymerType: helmTypeToPolymerType(seqMonomer.biotype)})!;
      lastMonomerCappingAtom = monomer.terminalR2Atom;
      const mAtomFirst = v.nodeShift;
      const mBondFirst = v.bondShift;
      addMonomerToMolblock(monomer, molfileAtomBlock, molfileBondBlock, v, LC);
      //adding stereo atoms to array for further STEABS block generation
      monomer.stereoAtoms?.forEach((i) => steabsCollection.push(i + nAtoms));
      nAtoms += monomer.atoms.x.length;

      const mAtomCount = v.nodeShift - mAtomFirst;
      const mAtomList: number[] = new Array<number>(mAtomCount);
      for (let maI = 0; maI < mAtomCount; ++maI) mAtomList[maI] = mAtomFirst + maI;

      const mBondCount = v.bondShift - mBondFirst;
      const mBondList: number[] = new Array<number>(mBondCount);
      for (let mbI = 0; mbI < mBondCount; ++mbI) mBondList[mbI] = mBondFirst + mbI;

      monomers.set(v.i, {
        biotype: seqMonomer.biotype,
        symbol: seqMonomer.symbol,
        atoms: mAtomList, bonds: mBondList
      });
    }
  }

  // if the last monomer needs to be capped, add the terminal OH to the resulting molfile
  if (needsCapping)
    capResultingMolblock(molfileAtomBlock, molfileBondBlock, v, LC, lastMonomerCappingAtom ?? C.OXYGEN);

  const molfileCountsLine = C.V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + C.V3K_COUNTS_LINE_ENDING;

  // todo: possible optimization may be achieved by replacing .join('') with +=
  // since counterintuitively joining an array into a new string is reportedly
  // slower than using += as below

  let result = '';
  result += C.V3K_HEADER_FIRST_LINE;
  result += C.V3K_HEADER_SECOND_LINE;
  result += C.V3K_BEGIN_CTAB_BLOCK;
  result += molfileCountsLine;
  result += C.V3K_BEGIN_ATOM_BLOCK;
  result += molfileAtomBlock.join('');
  result += C.V3K_END_ATOM_BLOCK;
  result += C.V3K_BEGIN_BOND_BLOCK;
  result += molfileBondBlock.join('');
  result += C.V3K_END_BOND_BLOCK;
  if (steabsCollection.length > 0)
    result += getCollectionBlock(steabsCollection);
  result += C.V3K_END_CTAB_BLOCK;
  result += C.V3K_END;

  // return molfileParts.join('');
  return {molfile: result, monomers: monomers};
}


function getCollectionBlock(collection: number[]): string {
  //one row in STEABS block can be no longer than 80 symbols
  //maxSymbols = 80 symbols minus ' -\n' (4 symbols)
  const maxSymbols = 76;
  const rowsArray = [];

  let newCollectionRow = `M  V30 MDLV30/STEABS ATOMS=(${collection.length}`;
  for (let i = 0; i < collection.length; i++) {
    const updatedRow = `${newCollectionRow} ${collection[i]}`;
    if (updatedRow.length > maxSymbols) {
      rowsArray.push(`${newCollectionRow} -\n`);
      newCollectionRow = `M  V30 ${collection[i]}`;
    } else
      newCollectionRow = updatedRow;
    //in case last atom was added - close the block
    if (i === collection.length - 1)
      rowsArray.push(`${newCollectionRow})\n`);
  }
  return `M  V30 BEGIN COLLECTION\n${rowsArray.join('')}M  V30 END COLLECTION\n`;
}

/** Cap the resulting (after sewing up all the monomers) molfile with 'O'
 * @param {string[]} molfileAtomBlock - Array of lines of the resulting molfile atom block
 * @param {string[]} molfileBondBlock - Array of lines of the resulting molfile bond block
 * @param {LoopVariables} v - Loop variables
 * @param {LoopConstants} LC - Loop constants*/
function capResultingMolblock(
  molfileAtomBlock: string[], molfileBondBlock: string[],
  v: LoopVariables, LC: LoopConstants, cappingAtomType: string = C.OXYGEN
): void {
  // add terminal oxygen
  const atomIdx = v.nodeShift + 1;
  molfileAtomBlock[LC.atomCount] = C.V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
    (cappingAtomType ?? C.OXYGEN) + ' ' + keepPrecision(v.backbonePositionShift[0]) + ' ' +
    v.flipFactor * keepPrecision(v.backbonePositionShift[1]) + ' ' + '0.000000 0' + '\n';

  // add terminal bond
  const firstAtom = v.backboneAttachNode;
  const secondAtom = atomIdx;
  molfileBondBlock[LC.bondCount] = C.V3K_BEGIN_DATA_LINE + v.bondShift + ' ' +
    1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
}

function addAminoAcidToMolblock(monomer: MolGraph, molfileAtomBlock: string[],
  molfileBondBlock: string[], v: LoopVariables
): void {
  v.flipFactor = (-1) ** (v.i % 2); // to flip every even monomer over OX
  addBackboneMonomerToMolblock(monomer, molfileAtomBlock, molfileBondBlock, v);
}

function addBackboneMonomerToMolblock(
  monomer: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables
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
  updateChainExtendingVariables(monomer, v);
}

function addNucleotideToMolblock(
  nucleobase: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables, LC: LoopConstants
): void {
  // construnct the lines of V3K molfile atom block corresponding to phosphate
  // and sugar
  if (v.i === 0)
    addBackboneMonomerToMolblock(LC.sugar!, molfileAtomBlock, molfileBondBlock, v);
  else {
    for (const monomer of [LC.phosphate, LC.sugar])
      addBackboneMonomerToMolblock(monomer!, molfileAtomBlock, molfileBondBlock, v);
  }

  addBranchMonomerToMolblock(nucleobase, molfileAtomBlock, molfileBondBlock, v);
}

function addBranchMonomerToMolblock(
  monomer: MolGraph, molfileAtomBlock: string[], molfileBondBlock: string[], v: LoopVariables
): void {
  fillBranchAtomLines(monomer, molfileAtomBlock, v);
  fillBondLines(monomer, molfileBondBlock, v);
  fillBackboneToBranchBond(monomer, molfileBondBlock, v);

  // C-N bond
  const bondIdx = v.bondShift;
  const firstAtom = v.branchAttachNode;
  const secondAtom = monomer.meta.terminalNodes[0] + v.nodeShift;
  molfileBondBlock[bondIdx - 1] = C.V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
    1 + ' ' + firstAtom + ' ' + secondAtom + '\n';

  // update loop variables
  v.bondShift += monomer.bonds.atomPairs.length + 1;
  v.nodeShift += monomer.atoms.atomTypes.length;
}

function updateChainExtendingVariables(monomer: MolGraph, v: LoopVariables): void {
  v.backboneAttachNode = v.nodeShift + monomer.meta.terminalNodes[1];
  v.bondShift += monomer.bonds.atomPairs.length + 1;

  v.nodeShift += monomer.atoms.atomTypes.length;
  v.backbonePositionShift[0] += monomer.meta.backboneShift?.[0] ?? 0; // todo: non-null check
  v.backbonePositionShift[1] += v.flipFactor * (monomer.meta.backboneShift?.[1] ?? 0);
}

function updateBranchVariables(monomer: MolGraph, v: LoopVariables): void {
  v.branchAttachNode = v.nodeShift + monomer.meta.terminalNodes[2];
  for (let i = 0; i < 2; ++i)
    v.branchPositionShift[i] = v.backbonePositionShift[i] + monomer.meta.branchShift![i];
}

function fillAtomLines(monomer: MolGraph, molfileAtomBlock: string[], v: LoopVariables): void {
  for (let j = 0; j < monomer.atoms.atomTypes.length; ++j) {
    const atomIdx = v.nodeShift + j + 1;
    molfileAtomBlock[v.nodeShift + j] = C.V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
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
    molfileAtomBlock[v.nodeShift + j] = C.V3K_BEGIN_DATA_LINE + atomIdx + ' ' +
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
    molfileBondBlock[v.bondShift + j] = C.V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
      monomer.bonds.bondTypes[j] + ' ' +
      firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
  }
}

function fillChainExtendingBond(monomer: MolGraph, molfileBondBlock: string[], v: LoopVariables): void {
  if (v.backboneAttachNode !== 0) {
    const bondIdx = v.bondShift;
    const firstAtom = v.backboneAttachNode;
    const secondAtom = monomer.meta.terminalNodes[0] + v.nodeShift;
    molfileBondBlock[v.bondShift - 1] = C.V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
      1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
  }
}

// todo: remove
function fillBackboneToBranchBond(branchMonomer: MolGraph, molfileBondBlock: string[], v: LoopVariables): void {
  const bondIdx = v.bondShift;
  const firstAtom = v.branchAttachNode;
  const secondAtom = branchMonomer.meta.terminalNodes[0] + v.nodeShift;
  molfileBondBlock[bondIdx - 1] = C.V3K_BEGIN_DATA_LINE + bondIdx + ' ' +
    1 + ' ' + firstAtom + ' ' + secondAtom + '\n';
}

/** Compute the atom/bond counts for the resulting molfile, depending on the
 * type of polymer (peptide/nucleotide)
 * @param {ISeqMonomer[]} monomerSeq - the sequence of monomers
 * @param {MonomerMolGraphMap} monomersDict - the dictionary of monomers
 * @param {ALPHABET} alphabet - the alphabet of the monomers
 * @param {PolymerType} polymerType - the type of polymer
 * @param {boolean} triplesMode - true when monomerSeq is a flat list of HELM RNA triples
 * @return {Object} the atom/bond counts plus needsCapping flag */
function getResultingAtomBondCounts(
  monomerSeq: ISeqMonomer[], monomersDict: MonomerMolGraphMap,
  alphabet: ALPHABET, polymerType: PolymerType,
  triplesMode: boolean
): { atomCount: number, bondCount: number, needsCapping: boolean } {
  let atomCount = 0;
  let bondCount = 0;

  let monomerCount: number = 0;
  let needsCapping = true;
  let lastMonomerGraph: MolGraph | null = null;
  let lastPhosphateGraph: MolGraph | null = null;
  // sum up all the atoms/nodes provided by the sequence
  for (let i = 0; i < monomerSeq.length; ++i) {
    const seqMonomer = monomerSeq[i];
    if (seqMonomer.symbol === GAP_SYMBOL) continue; // Skip for gap/empty monomer in MSA
    if (seqMonomer.symbol == '*')
      throw new Error(`Gap canonical symbol is '', not '*`);
    lastMonomerGraph = getMolGraph(monomersDict, {symbol: seqMonomer.symbol, polymerType: helmTypeToPolymerType(seqMonomer.biotype)})!;
    atomCount += lastMonomerGraph.atoms.x.length;
    bondCount += lastMonomerGraph.bonds.bondTypes.length;
    monomerCount++;
    // In triples mode, every 3rd entry (index 2 mod 3) is a phosphate. Track
    // the LAST one — its atoms/bonds are dropped at the 3'-terminus.
    if (triplesMode && i % 3 === 2) lastPhosphateGraph = lastMonomerGraph;
  }

  // add extra values depending on the polymer type
  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
    // add the rightmost/terminating cap group 'OH' (i.e. 'O')
    atomCount += 1;
    // add chain-extending bonds (C-NH per each monomer pair and terminal C-OH)
    bondCount += monomerCount;
    // if the last monomer is something like NH2, which only has R1, there is no need to cap it
    // although, this should never happen, but hey... in other bits of code, there is a chunk that adds pseudo-R2 as hydrogen
    // we should also check, if the R2 of the last monomer is not hydrogen, that case should also be omitted
    if (monomerCount > 0) {
      if ((lastMonomerGraph?.meta?.rNodes?.length ?? 0) < 2 || lastMonomerGraph?.terminalR2Atom?.toLowerCase() === C.HYDROGEN.toLowerCase()) {
        needsCapping = false;
        atomCount -= 1; // remove the last atom (the terminal 'O')
        bondCount -= 1; // remove the last bond (the terminal C-OH)
      }
    }
  } else if (triplesMode) { // nucleotides — HELM triples (per-position sugar/base/phosphate)
    void lastPhosphateGraph; // kept for symmetry; trailing P is now retained when HELM wrote it

    // Per-monomer-loop already summed sugars + bases + phosphates for every
    // entry in monomerSeq. Reservation: each backbone (sugar/phosphate)
    // reserves +1 bond slot for the chain-extending bond, each branch (base)
    // reserves +1 for its branch bond. Total reservations = monomerCount.
    //   length === 3N      → trailing P emitted, monomerCount = 3N.
    //   length === 3N - 1  → trailing P absent, monomerCount = 3N - 1.
    // Either way, bondCount = sum_bonds + monomerCount.
    atomCount += 1; // OH cap atom (rides on trailing P or on last sugar's R2)
    bondCount += monomerCount;
  } else { // nucleotides — bases-only legacy path with default sugar/phosphate
    const sugar = (alphabet === ALPHABET.DNA) ?
      getMolGraph(monomersDict, C.DEOXYRIBOSE)! : getMolGraph(monomersDict, C.RIBOSE)!;
    const phosphate = getMolGraph(monomersDict, C.PHOSPHATE)!;

    // add phosphate per each pair of nucleobase symbols
    atomCount += (monomerSeq.length - 1) * phosphate.atoms.x.length;

    // add sugar per each nucleobase symbol
    atomCount += monomerSeq.length * sugar.atoms.x.length;

    // add the leftmost cap group 'OH' (i.e. 'O')
    atomCount += 1;

    // add bonds from phosphate monomers
    bondCount += (monomerSeq.length - 1) * phosphate.bonds.bondTypes.length;

    // add bonds from sugar monomers
    bondCount += monomerSeq.length * sugar.bonds.bondTypes.length;

    // exclude the first chain-extending bond O-P (absent, no 'leftmost' phosphate)
    bondCount -= 1;

    // add chain-extending and branch bonds (O-P, C-O and C-N per each nucleotide)
    bondCount += monomerSeq.length * 3;
  }

  return {atomCount, bondCount, needsCapping};
}

/** Triples-mode RNA assembly. Iterate by NUCLEOTIDE; each nucleotide normally
 * contributes 3 entries from monomerSeq in order [sugar, base, phosphate].
 * The LAST nucleotide may omit its trailing phosphate (length === 3N - 1) —
 * in that case HELM specified an unphosphorylated 3'-end and we cap on the
 * last sugar's R2 instead of on a final phosphate.
 *
 * @param {ISeqMonomer[]} monomerSeq - flat sugar/base/[phosphate] entries
 * @param {MonomerMolGraphMap} monomersDict - monomer dictionary
 * @param {string[]} molfileAtomBlock - output atom block buffer
 * @param {string[]} molfileBondBlock - output bond block buffer
 * @param {LoopVariables} v - mutable loop variables
 * @param {LoopConstants} LC - loop constants (seqLength = nucleotide count)
 * @param {MonomerMap} monomers - per-position monomer index map (filled in)
 * @param {number[]} steabsCollection - stereo atoms collected for STEABS block
 * @param {Function} addAtoms - callback to bump per-row atom counter for stereo offsets
 * @param {Function} getAtoms - callback to read current atom counter */
function runTriplesAssembly(
  monomerSeq: ISeqMonomer[], monomersDict: MonomerMolGraphMap,
  molfileAtomBlock: string[], molfileBondBlock: string[],
  v: LoopVariables, LC: LoopConstants,
  monomers: MonomerMap, steabsCollection: number[],
  addAtoms: (n: number) => void, getAtoms: () => number
): void {
  const N = LC.seqLength;
  // length === 3N → trailing P present. length === 3N-1 → trailing P absent.
  const hasTrailingP = monomerSeq.length === 3 * N;
  for (let n = 0; n < N; ++n) {
    v.i = n;
    const sugarSm = monomerSeq[3 * n];
    const baseSm = monomerSeq[3 * n + 1];

    const sugarG = getMolGraph(monomersDict, {symbol: sugarSm.symbol, polymerType: helmTypeToPolymerType(sugarSm.biotype)})!;
    const baseG = getMolGraph(monomersDict, {symbol: baseSm.symbol, polymerType: helmTypeToPolymerType(baseSm.biotype)})!;
    // The previous nucleotide's trailing phosphate is the linker into THIS
    // nucleotide's sugar. It always exists for n >= 1 because only the
    // very LAST nucleotide is permitted to omit its phosphate.
    const prevPhosG = (n === 0) ? null :
      getMolGraph(monomersDict, {symbol: monomerSeq[3 * (n - 1) + 2].symbol,
        polymerType: helmTypeToPolymerType(monomerSeq[3 * (n - 1) + 2].biotype)})!;

    // emit prevPhosphate (chain linker), then sugar, then base
    const seqStartIdx = 3 * n; // for keying the per-row monomer map by HELM splitter index
    if (prevPhosG) {
      const prevAtomFirst = v.nodeShift;
      const prevBondFirst = v.bondShift;
      addBackboneMonomerToMolblock(prevPhosG, molfileAtomBlock, molfileBondBlock, v);
      prevPhosG.stereoAtoms?.forEach((i) => steabsCollection.push(i + getAtoms()));
      addAtoms(prevPhosG.atoms.x.length);
      // The phosphate "belongs to" nucleotide n-1's HELM index (3*(n-1)+2).
      const prevPhosKey = 3 * (n - 1) + 2;
      const aList: number[] = [];
      for (let a = prevAtomFirst; a < v.nodeShift; ++a) aList.push(a);
      const bList: number[] = [];
      for (let b = prevBondFirst; b < v.bondShift; ++b) bList.push(b);
      monomers.set(prevPhosKey, {
        biotype: monomerSeq[prevPhosKey].biotype,
        symbol: monomerSeq[prevPhosKey].symbol,
        atoms: aList, bonds: bList,
      });
    }

    const sugarAtomFirst = v.nodeShift;
    const sugarBondFirst = v.bondShift;
    addBackboneMonomerToMolblock(sugarG, molfileAtomBlock, molfileBondBlock, v);
    sugarG.stereoAtoms?.forEach((i) => steabsCollection.push(i + getAtoms()));
    addAtoms(sugarG.atoms.x.length);
    const sAList: number[] = [];
    for (let a = sugarAtomFirst; a < v.nodeShift; ++a) sAList.push(a);
    const sBList: number[] = [];
    for (let b = sugarBondFirst; b < v.bondShift; ++b) sBList.push(b);
    monomers.set(seqStartIdx, {
      biotype: sugarSm.biotype, symbol: sugarSm.symbol,
      atoms: sAList, bonds: sBList,
    });

    const baseAtomFirst = v.nodeShift;
    const baseBondFirst = v.bondShift;
    addBranchMonomerToMolblock(baseG, molfileAtomBlock, molfileBondBlock, v);
    baseG.stereoAtoms?.forEach((i) => steabsCollection.push(i + getAtoms()));
    addAtoms(baseG.atoms.x.length);
    const bAList: number[] = [];
    for (let a = baseAtomFirst; a < v.nodeShift; ++a) bAList.push(a);
    const bBList: number[] = [];
    for (let b = baseBondFirst; b < v.bondShift; ++b) bBList.push(b);
    monomers.set(seqStartIdx + 1, {
      biotype: baseSm.biotype, symbol: baseSm.symbol,
      atoms: bAList, bonds: bBList,
    });
  }

  // Emit the trailing phosphate (the LAST nucleotide's P) as the final
  // backbone monomer ONLY if HELM explicitly wrote it. The OH cap then
  // rides on its free oxygen (3'-phosphate-OH terminus). If HELM omitted
  // the trailing phosphate, the cap rides on the last sugar's R2 — that
  // requires no extra emission here, since v.backboneAttachNode already
  // points at the last sugar's R2 attach.
  if (hasTrailingP) {
    const lastPhosKey = 3 * (N - 1) + 2;
    const lastPhosSm = monomerSeq[lastPhosKey];
    const lastPhosG = getMolGraph(monomersDict,
      {symbol: lastPhosSm.symbol, polymerType: helmTypeToPolymerType(lastPhosSm.biotype)})!;
    const trailAtomFirst = v.nodeShift;
    const trailBondFirst = v.bondShift;
    addBackboneMonomerToMolblock(lastPhosG, molfileAtomBlock, molfileBondBlock, v);
    lastPhosG.stereoAtoms?.forEach((i) => steabsCollection.push(i + getAtoms()));
    addAtoms(lastPhosG.atoms.x.length);
    const tAList: number[] = [];
    for (let a = trailAtomFirst; a < v.nodeShift; ++a) tAList.push(a);
    const tBList: number[] = [];
    for (let b = trailBondFirst; b < v.bondShift; ++b) tBList.push(b);
    monomers.set(lastPhosKey, {
      biotype: lastPhosSm.biotype, symbol: lastPhosSm.symbol,
      atoms: tAList, bonds: tBList,
    });
  }
}

/** Keep precision upon floating point operations over atom coordinates
 * @param {number}x - the floating point number
 * @return {number} - the floating point number with the same precision
 */
export function keepPrecision(x: number): number {
  return Math.round(C.PRECISION_FACTOR * x) / C.PRECISION_FACTOR;
}

