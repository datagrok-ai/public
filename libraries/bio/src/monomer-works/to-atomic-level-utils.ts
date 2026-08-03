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
      polymerType === HELM_POLYMER_TYPE.RNA || (
        polymerType === HELM_POLYMER_TYPE.PEPTIDE &&
      it[HELM_FIELDS.MONOMER_TYPE] !== HELM_MONOMER_TYPE.BRANCH
      )) {
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
 * @param {number[]} chainStarts - Optional 0-based start positions of disjoint HELM chains.
 *   When more than one chain is present each chain is assembled independently and the
 *   chains are stacked vertically into a single molfile (multi-strand HELM support).
 * @return {MolfileWithMap} - Molfile V3000 + per-position monomer index map */
export function monomerSeqToMolfile(
  monomerSeq: ISeqMonomer[], monomersDict: MonomerMolGraphMap,
  alphabet: ALPHABET, polymerType: PolymerType,
  roles?: NucleotideRole[], chainStarts?: number[]
): MolfileWithMap {
  if (monomerSeq.length === 0) {
    // throw new Error('monomerSeq is empty');
    return MolfileWithMap.createEmpty();
  }

  // Multi-chain HELM (e.g. RNA1{...}|RNA2{...}): assemble each disjoint chain
  // on its own, then stack them vertically into a single molfile. Each chain
  // reuses the exact single-chain assembly below, so its geometry, capping and
  // per-monomer atom/bond map are unchanged; combineMolfilesVertically shifts
  // the coordinates and re-bases the atom/bond indices and monomer-map keys.
  if (chainStarts && chainStarts.length > 1) {
    const parts: MolfileWithMap[] = new Array(chainStarts.length);
    for (let c = 0; c < chainStarts.length; ++c) {
      const s = chainStarts[c];
      const e = c + 1 < chainStarts.length ? chainStarts[c + 1] : monomerSeq.length;
      const subSeq = monomerSeq.slice(s, e);
      const subRoles = roles ? roles.slice(s, e) : undefined;
      parts[c] = monomerSeqToMolfile(subSeq, monomersDict, alphabet, polymerType, subRoles);
    }
    return combineMolfilesVertically(parts, chainStarts);
  }

  // Triples mode is on only when the caller flagged the row with roles
  // (built and validated by `buildRolesForHelmRna` in to-atomic-level.ts).
  // The roles array carries the per-position semantics — including
  // TERMINAL_5P / TERMINAL_3P for non-canonical chain ends — so we don't
  // re-validate the length here.
  const triplesMode = polymerType === HELM_POLYMER_TYPE.RNA && !!roles &&
    roles.length === monomerSeq.length;

  // define atom and bond counts, taking into account the bond type
  const {atomCount, bondCount, needsCapping} =
    getResultingAtomBondCounts(monomerSeq, monomersDict, alphabet, polymerType, triplesMode, roles);

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
      monomerSeq, roles!, monomersDict, molfileAtomBlock, molfileBondBlock,
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
 * @param {Array} roles - per-position role tags (only meaningful when triplesMode is true)
 * @return {Object} the atom/bond counts plus needsCapping flag */
function getResultingAtomBondCounts(
  monomerSeq: ISeqMonomer[], monomersDict: MonomerMolGraphMap,
  alphabet: ALPHABET, polymerType: PolymerType,
  triplesMode: boolean, roles?: NucleotideRole[]
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
    // entry in monomerSeq (terminals included). Reservation: each backbone
    // emit (sugar / phosphate / terminal) reserves +1 bond slot for the
    // chain-extending bond, each branch (base) reserves +1 for its branch
    // bond. Total reservations = monomerCount.
    bondCount += monomerCount;

    // OH cap atom rides on the 3'-end. Skip it when HELM specified a
    // 3'-terminal modifier (e.g. GalNAc) that IS the chain end.
    const has3pTerm = !!roles && roles.length > 0 &&
      roles[roles.length - 1] === NucleotideRole.TERMINAL_3P;
    if (has3pTerm) {
      needsCapping = false;
      // Mirror the peptide branch (which does both atomCount-- and bondCount--
      // when it skips the terminal cap). `bondCount += monomerCount` above
      // reserves one chain-extending slot per monomer; the LAST monomer's slot
      // is normally filled by the terminal OH cap bond. With a 3'-terminal
      // modifier there is no cap, so that last reserved slot stays empty —
      // leaving the declared bond count one higher than the emitted bond
      // lines. Drop it so the V3000 COUNTS line matches the bond block exactly
      // (otherwise the pre-OCL molfile is malformed and only survives because
      // the OCL chirality pass re-derives the counts).
      bondCount -= 1;
    } else
      atomCount += 1; // OH cap atom (rides on trailing P or on last sugar's R2)
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

// Role-driven RNA assembly. Walks the monomers in chain order and emits
// each one according to its NucleotideRole (assigned by `buildRolesForHelmRna`
// from the library, not from a fixed triple index):
//   - BASE                                  → branch monomer, attached to the
//                                             branch point of the sugar most
//                                             recently emitted.
//   - SUGAR / PHOSPHATE / TERMINAL_5P / 3P  → backbone monomer, chained to the
//                                             previous backbone unit.
//
// Because the role comes from chemistry, this transparently handles every
// backbone layout — standard [sugar, base, phosphate] triples, a 5'-leading
// phosphate, several phosphates / linkers in a row, a linker dropped in the
// middle of the chain, missing trailing phosphate, and 5'/3' terminal
// modifiers — without any index arithmetic. The first backbone monomer
// naturally has no incoming chain bond (v.backboneAttachNode starts at 0),
// and the trailing OH cap is added by the caller via needsCapping (skipped
// when a TERMINAL_3P ends the chain).
function runTriplesAssembly(
  monomerSeq: ISeqMonomer[], roles: NucleotideRole[],
  monomersDict: MonomerMolGraphMap,
  molfileAtomBlock: string[], molfileBondBlock: string[],
  v: LoopVariables, LC: LoopConstants,
  monomers: MonomerMap, steabsCollection: number[],
  addAtoms: (n: number) => void, getAtoms: () => number
): void {
  void LC; // assembly is driven by per-monomer roles, not by LC.seqLength
  for (let i = 0; i < monomerSeq.length; ++i) {
    const sm = monomerSeq[i];
    if (sm.symbol === GAP_SYMBOL) continue;
    const role = roles[i];
    const mG = getMolGraph(monomersDict,
      {symbol: sm.symbol, polymerType: helmTypeToPolymerType(sm.biotype)})!;

    const aFirst = v.nodeShift;
    const bFirst = v.bondShift;
    v.i = i; // keep loop counter monotone (unused by nucleotide geometry)
    if (role === NucleotideRole.BASE)
      // Branch: attaches to v.branchAttachNode, set when the preceding sugar
      // was emitted. Does not advance the backbone, so the next backbone unit
      // still chains from that sugar's 3' side.
      addBranchMonomerToMolblock(mG, molfileAtomBlock, molfileBondBlock, v);
    else
      // Backbone: sugar / phosphate / linker / terminal modifier. Chains from
      // the previous backbone unit (none for the very first one) and, for a
      // sugar, sets up the branch attach point for the next base.
      addBackboneMonomerToMolblock(mG, molfileAtomBlock, molfileBondBlock, v);

    mG.stereoAtoms?.forEach((s) => steabsCollection.push(s + getAtoms()));
    addAtoms(mG.atoms.x.length);

    const aList: number[] = [];
    for (let a = aFirst; a < v.nodeShift; ++a) aList.push(a);
    const bList: number[] = [];
    for (let b = bFirst; b < v.bondShift; ++b) bList.push(b);
    monomers.set(i, {biotype: sm.biotype, symbol: sm.symbol, atoms: aList, bonds: bList});
  }
}

/** Keep precision upon floating point operations over atom coordinates
 * @param {number}x - the floating point number
 * @return {number} - the floating point number with the same precision
 */
export function keepPrecision(x: number): number {
  return Math.round(C.PRECISION_FACTOR * x) / C.PRECISION_FACTOR;
}

// Vertical gap (in molfile coordinate units) inserted between stacked chains
// of a multi-strand HELM. Roughly two bond lengths — enough to read the strands
// as separate without pushing them far apart.
const CHAIN_VERTICAL_GAP = 3.0;

const STEABS_MARKER = 'M  V30 MDLV30/STEABS ATOMS=(';

type ParsedAtom = { type: string, xStr: string, y: number, tail: string };
type ParsedBond = { btype: string, a1: number, a2: number, tail: string };
type ParsedMol = { atoms: ParsedAtom[], bonds: ParsedBond[], steabs: number[] };

/** Split `body` into its first `n` space-separated tokens plus a trailing
 * remainder token (so the result has n+1 entries). The remainder keeps any
 * embedded spaces of the tail verbatim.
 * @param {string} body - the line body to split
 * @param {number} n - number of leading tokens to peel off
 * @return {string[]} - the n tokens followed by the remainder */
function splitLeadingTokens(body: string, n: number): string[] {
  const out: string[] = [];
  let s = body;
  for (let i = 0; i < n; ++i) {
    const sp = s.indexOf(' ');
    if (sp === -1) { out.push(s); s = ''; continue; }
    out.push(s.substring(0, sp));
    s = s.substring(sp + 1);
  }
  out.push(s);
  return out;
}

/** Extract the data lines (`M  V30 ...`) between a BEGIN/END marker pair.
 * @param {string} molfile - the V3000 molfile
 * @param {string} beginMarker - the block BEGIN marker line
 * @param {string} endMarker - the block END marker line
 * @return {string[]} - the data lines of the block, newline stripped */
function extractBlockDataLines(molfile: string, beginMarker: string, endMarker: string): string[] {
  const begin = molfile.indexOf(beginMarker);
  if (begin === -1) return [];
  const start = molfile.indexOf('\n', begin) + 1;
  const end = molfile.indexOf(endMarker, start);
  return molfile.substring(start, end === -1 ? molfile.length : end)
    .split('\n')
    .filter((l) => l.startsWith(C.V3K_BEGIN_DATA_LINE))
    .map((l) => l.substring(C.V3K_BEGIN_DATA_LINE.length));
}

/** Parse a molfile produced by monomerSeqToMolfile back into atoms / bonds /
 * stereo-atom indices. Robust because we control the exact serialization
 * format; only the fields we need to renumber (indices) and shift (y) are
 * parsed, everything past them is preserved verbatim as an opaque tail.
 * @param {string} molfile - the V3000 molfile
 * @return {ParsedMol} - parsed atoms, bonds and STEABS atom indices */
function parseGeneratedMolfile(molfile: string): ParsedMol {
  const atoms: ParsedAtom[] = [];
  for (const body of extractBlockDataLines(molfile, C.V3K_BEGIN_ATOM_BLOCK, C.V3K_END_ATOM_BLOCK)) {
    // body: "<idx> <type> <x> <y> <z> <...kwargs>"
    const t = splitLeadingTokens(body, 4); // [idx, type, x, y, tail]
    atoms.push({type: t[1], xStr: t[2], y: parseFloat(t[3]), tail: t[4]});
  }

  const bonds: ParsedBond[] = [];
  for (const body of extractBlockDataLines(molfile, C.V3K_BEGIN_BOND_BLOCK, C.V3K_END_BOND_BLOCK)) {
    // body: "<idx> <btype> <a1> <a2> <...cfg/kwargs>"
    const t = splitLeadingTokens(body, 4); // [idx, btype, a1, a2, tail]
    bonds.push({btype: t[1], a1: parseInt(t[2]), a2: parseInt(t[3]), tail: t[4]});
  }

  const steabs: number[] = [];
  let idx = molfile.indexOf(STEABS_MARKER);
  while (idx !== -1) {
    idx += STEABS_MARKER.length;
    const closeIdx = molfile.indexOf(')', idx);
    // first entry inside the parens is the atom count, not an atom index
    steabs.push(...molfile.substring(idx, closeIdx).split(' ').slice(1)
      .map((s) => parseInt(s)).filter((n) => !Number.isNaN(n)));
    idx = molfile.indexOf(STEABS_MARKER, closeIdx);
  }

  return {atoms, bonds, steabs};
}

/** Combine several already-assembled single-chain molfiles into one molfile,
 * stacking them vertically: the first chain keeps its coordinates, every
 * subsequent chain is shifted straight down so its top-most atom sits just
 * below the running bottom of everything placed so far. Atom / bond indices
 * and per-monomer maps are re-based, STEABS collections merged, and the
 * monomer-map keys offset by each chain's start position so the combined map
 * still keys by original flat HELM position.
 * @param {MolfileWithMap[]} parts - per-chain assembled molfiles (chain order)
 * @param {number[]} chainStarts - flat start position of each chain
 * @return {MolfileWithMap} - single combined molfile + merged monomer map */
export function combineMolfilesVertically(parts: MolfileWithMap[], chainStarts: number[]): MolfileWithMap {
  const combinedAtomLines: string[] = [];
  const combinedBondLines: string[] = [];
  const combinedSteabs: number[] = [];
  const combinedMap: MonomerMap = new MonomerMap();

  let atomOffset = 0;
  let bondOffset = 0;
  let runningBottom = Number.POSITIVE_INFINITY; // lowest y placed so far
  let anyPlaced = false;

  for (let c = 0; c < parts.length; ++c) {
    const part = parts[c];
    if (!part || !part.molfile) continue;
    const pm = parseGeneratedMolfile(part.molfile);
    if (pm.atoms.length === 0) continue;

    let minY = Number.POSITIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;
    for (const a of pm.atoms) {
      if (a.y < minY) minY = a.y;
      if (a.y > maxY) maxY = a.y;
    }

    // First placed chain keeps its coordinates; later chains drop below the
    // current bottom, leaving CHAIN_VERTICAL_GAP between the strands.
    const yShift = anyPlaced ? keepPrecision(runningBottom - CHAIN_VERTICAL_GAP - maxY) : 0;

    for (let ai = 0; ai < pm.atoms.length; ++ai) {
      const a = pm.atoms[ai];
      const newIdx = atomOffset + ai + 1;
      const newY = keepPrecision(a.y + yShift);
      combinedAtomLines.push(
        C.V3K_BEGIN_DATA_LINE + newIdx + ' ' + a.type + ' ' + a.xStr + ' ' + newY + ' ' + a.tail + '\n');
    }

    for (let bi = 0; bi < pm.bonds.length; ++bi) {
      const b = pm.bonds[bi];
      const newIdx = bondOffset + bi + 1;
      const tail = b.tail ? ' ' + b.tail : '';
      combinedBondLines.push(
        C.V3K_BEGIN_DATA_LINE + newIdx + ' ' + b.btype + ' ' + (b.a1 + atomOffset) + ' ' + (b.a2 + atomOffset) + tail + '\n');
    }

    for (const s of pm.steabs) combinedSteabs.push(s + atomOffset);

    const posOffset = chainStarts[c] ?? 0;
    for (const [pos, val] of part.monomers) {
      combinedMap.set(posOffset + pos, {
        biotype: val.biotype,
        symbol: val.symbol,
        atoms: val.atoms.map((x) => x + atomOffset),
        bonds: val.bonds.map((x) => x + bondOffset),
      });
    }

    const newBottom = keepPrecision(minY + yShift);
    runningBottom = anyPlaced ? Math.min(runningBottom, newBottom) : newBottom;
    anyPlaced = true;
    atomOffset += pm.atoms.length;
    bondOffset += pm.bonds.length;
  }

  if (!anyPlaced) return MolfileWithMap.createEmpty();

  let result = '';
  result += C.V3K_HEADER_FIRST_LINE;
  result += C.V3K_HEADER_SECOND_LINE;
  result += C.V3K_BEGIN_CTAB_BLOCK;
  result += C.V3K_BEGIN_COUNTS_LINE + atomOffset + ' ' + bondOffset + C.V3K_COUNTS_LINE_ENDING;
  result += C.V3K_BEGIN_ATOM_BLOCK;
  result += combinedAtomLines.join('');
  result += C.V3K_END_ATOM_BLOCK;
  result += C.V3K_BEGIN_BOND_BLOCK;
  result += combinedBondLines.join('');
  result += C.V3K_END_BOND_BLOCK;
  if (combinedSteabs.length > 0)
    result += getCollectionBlock(combinedSteabs);
  result += C.V3K_END_CTAB_BLOCK;
  result += C.V3K_END;

  return {molfile: result, monomers: combinedMap};
}

