/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {HELM_FIELDS, HELM_MONOMER_TYPE, HELM_POLYMER_TYPE, HELM_RGROUP_FIELDS} from '../utils/const';
import {ALPHABET, NOTATION} from '../utils/macromolecule/consts';
import {IMonomerLib, IMonomerLibBase, Monomer} from '../types/monomer-library';
import {getFormattedMonomerLib, keepPrecision} from './to-atomic-level-utils';
import {seqToMolFileWorker} from './seq-to-molfile';
import {Atoms, Bonds, hasMolGraph, ITypedArray, LibMonomerKey, MolGraph,
  MonomerMetadata, MonomerMolGraphMap, NucleotideRole, NumberWrapper, Point, setMolGraph} from './types';
import {ISeqHelper, ToAtomicLevelRes} from '../utils/seq-helper';
import {errInfo} from '../utils/err-info';
import {alphabetToPolymerType} from './utils';
import {HelmType, ISeqMonomer, PolymerType} from '../helm/types';
import {HelmTypes, PolymerTypes} from '../helm/consts';
import {monomerWorksConsts as C} from './consts';

// todo: verify that all functions have return types

/** Convert Macromolecule column into Molecule column storing molfile V3000 with the help of a monomer library
 * @param {DG.DataFrame} df - DataFrame containing the column to be converted
 * @param {DG.Column} seqCol - Column containing the macromolecule sequence
 * @param {IMonomerLib} monomerLib - Monomer library
 */
export async function _toAtomicLevel(
  df: DG.DataFrame, seqCol: DG.Column<string>,
  monomerLib: IMonomerLib, seqHelper: ISeqHelper, rdKitModule: RDModule
): Promise<ToAtomicLevelRes> {
  if (seqCol.semType !== DG.SEMTYPE.MACROMOLECULE) {
    const msg: string = `Only the ${DG.SEMTYPE.MACROMOLECULE} columns can be converted to atomic level, ` +
      `the chosen column has semType '${seqCol.semType}'`;
    return {molCol: null, warnings: [msg]};
  }

  let srcCol: DG.Column<string> = seqCol;
  const seqUh = seqHelper.getSeqHandler(seqCol);

  // Decide whether we must keep HELM as-is for nucleotide assembly. The
  // alphabet is unreliable here: a HELM RNA column with non-canonical bases
  // can detect as ALPHABET.UN. Inspect the splitter's per-monomer polymer
  // types instead — every monomer in every (non-empty) row must be RNA.
  const keepHelmTriples = seqUh.notation === NOTATION.HELM &&
    isAllHelmRna(seqCol, seqHelper);

  // convert 'helm' to 'separator' units (peptide HELM, anything non-RNA/DNA)
  if (!keepHelmTriples &&
    seqUh.notation !== NOTATION.SEPARATOR && seqUh.notation !== NOTATION.BILN) {
    srcCol = seqUh.convert(NOTATION.SEPARATOR, '.');
    srcCol.name = seqCol.name; // Replace converted col name 'separator(<original>)' to '<original>';
  }

  let polymerType: PolymerType;
  let alphabet: ALPHABET;
  try {
    if (keepHelmTriples) {
      // We already proved the column is HELM RNA. Pin the polymer type to
      // RNA without leaning on the alphabet detector (which can yield UN
      // when only modified monomers appear).
      polymerType = HELM_POLYMER_TYPE.RNA;
      alphabet = pickRnaAlphabetFromHelm(seqCol, seqHelper);
    } else {
      const srcSh = seqHelper.getSeqHandler(srcCol);
      alphabet = srcSh.alphabet as ALPHABET;
      polymerType = alphabetToPolymerType(alphabet);
    }
  } catch (err: any) {
    const [errMsg, _errStack] = errInfo(err);
    return {molCol: null, warnings: [errMsg]};
  }

  const monomerSequencesArray: ISeqMonomer[][] = getMonomerSequencesArray(srcCol, seqHelper);
  // Per-row disjoint-chain start positions. For HELM with several chains
  // (e.g. RNA1{...}|RNA2{...}) each row carries multiple starts; every other
  // notation (and single-chain HELM) carries [0]. Computed from srcCol — the
  // exact column that seqToMolFileWorker splits — so the indices line up.
  const chainStartsArray: number[][] = getChainStartsArray(srcCol, seqHelper);
  // Per-row roles: only set in HELM RNA/DNA mode where the splitter emits triples.
  const rolesArray: (NucleotideRole[] | undefined)[] = keepHelmTriples ?
    buildRolesForHelmRna(monomerSequencesArray, monomerLib, polymerType, chainStartsArray) :
    new Array(monomerSequencesArray.length).fill(undefined);
  // In triples mode, override the per-row biotype so it carries NUCLEOTIDE
  // (rather than whatever the seq-handler defaults to for an UN-alphabet column).
  if (keepHelmTriples) {
    for (const row of monomerSequencesArray)
      for (const m of row) m.biotype = HelmTypes.NUCLEOTIDE as HelmType;
  }
  const monomersDict = getMonomersDictFromLib(
    monomerSequencesArray, rolesArray, polymerType, alphabet, monomerLib, rdKitModule);
  const srcColLength = srcCol.length;

  const res = await seqToMolFileWorker(
    srcCol, monomersDict, alphabet, polymerType, monomerLib, seqHelper, rdKitModule, rolesArray, chainStartsArray);
  if (res.warnings.length > 0.05 * srcColLength)
    grok.shell.warning(`Molfile conversion resulted in ${res.warnings.length} errors`);

  // use chirality engine to fix chirality of linear sequences, i.e. add STEABS
  const chiralityFunc = DG.Func.find({name: 'convertToV3KViaOCL'})[0];
  if (chiralityFunc && res.molCol) {
    try {
      const mols = res.molCol.toList();
      const molsV3K = await chiralityFunc.apply({mols: mols});
      res.molCol.init((i) => {
        return molsV3K[i] ? molsV3K[i] : mols[i];
      });
    } catch (err: any) {
      console.error(err);
    }
  }

  return res;
}

// True iff every row in the column is HELM with a single disjoint
// sequence whose every monomer carries polymerType=RNA. The check goes
// through the splitter, not the alphabet, because non-canonical HELM RNA
// sequences can be mis-detected as ALPHABET.UN.
function isAllHelmRna(seqCol: DG.Column<string>, seqHelper: ISeqHelper): boolean {
  const sh = seqHelper.getSeqHandler(seqCol);
  if (sh.notation !== NOTATION.HELM) return false;

  let sawAny = false;
  for (let rowI = 0; rowI < seqCol.length; ++rowI) {
    const seqSS = sh.getSplitted(rowI);
    if (seqSS.length === 0) continue;
    const gi = seqSS.graphInfo;
    // graphInfo.polymerTypes is per-monomer. If the row has any non-RNA
    // monomer (peptide chain, CHEM, BLOB, ...), bail out.
    const pts = gi?.polymerTypes;
    if (!pts || pts.length !== seqSS.length) return false;
    for (const pt of pts)
      if (pt !== PolymerTypes.RNA) return false;
    sawAny = true;
  }
  return sawAny;
}

// When we've already confirmed all rows are HELM RNA via `isAllHelmRna`,
// pick a column-level alphabet for the (rare) bases-only fallback rows.
// Heuristic: if any sugar position references the deoxyribose symbol `d`,
// use DNA, otherwise RNA. Modified-only sequences default to RNA.
function pickRnaAlphabetFromHelm(seqCol: DG.Column<string>, seqHelper: ISeqHelper): ALPHABET {
  const sh = seqHelper.getSeqHandler(seqCol);
  for (let rowI = 0; rowI < seqCol.length; ++rowI) {
    const seqSS = sh.getSplitted(rowI);
    for (let i = 0; i < seqSS.length; ++i) {
      // Any deoxyribose sugar anywhere in the column → DNA. Don't rely on a
      // fixed triple index (i % 3 === 0): a leading phosphate, repeated
      // phosphates, or mid-chain linkers shift every monomer position, so the
      // sugar is no longer guaranteed to sit at index 0 modulo 3.
      if (seqSS.getCanonical(i) === C.DEOXYRIBOSE.symbol) return ALPHABET.DNA;
    }
  }
  return ALPHABET.RNA;
}

// Structural role of a single HELM RNA monomer, derived from its library
// definition rather than its position. The caller resolves 'terminal' to
// TERMINAL_5P / TERMINAL_3P from the monomer's position in the row.
type RnaRoleClass =
  NucleotideRole.SUGAR | NucleotideRole.BASE | NucleotideRole.PHOSPHATE | 'terminal' | null;

// Classify a HELM RNA monomer by its library properties (monomer type +
// R-groups), NOT by its position in the sequence:
//   - Branch monomers (A/C/G/U/T and modified bases) → BASE. They carry a
//     single R1 that attaches to the preceding sugar's R3 branch point.
//   - Backbone monomers carrying a branch R-group (R3) → SUGAR (ribose,
//     deoxyribose, and modified sugars such as fl2r / lna / m all expose
//     R1/R2/R3).
//   - Backbone monomers with two chain-extending R-groups and no branch →
//     PHOSPHATE / linker (p, sp, Rsp, …). This is what makes a leading
//     phosphate, several phosphates in a row, or a linker in the middle of
//     the chain just work: each is recognised as a chain-extending unit,
//     never mistaken for a sugar.
//   - Single-R-group non-Branch monomers (Chol / GalNAc / Bio, which the
//     library marks monomerType='Undefined') → 'terminal' chain-end
//     modifiers. The caller decides 5' vs 3' from the boundary.
//   - Unknown symbol (absent from lib) → null.
//
// Returning roles from chemistry instead of index removes the rigid
// "[sugar, base, phosphate] repeated" assumption that produced NaN
// coordinates whenever a phosphate did not sit where the triple index
// expected it.
function classifyRnaRole(
  monomerLib: IMonomerLibBase, polymerType: PolymerType, symbol: string
): RnaRoleClass {
  const m = monomerLib.getMonomer(polymerType, symbol);
  if (!m) return null;
  if (m[HELM_FIELDS.MONOMER_TYPE] === HELM_MONOMER_TYPE.BRANCH)
    return NucleotideRole.BASE;
  const rgroups = m.rgroups ?? [];
  // A chain-end modifier exposes exactly one attachment point.
  if (rgroups.length <= 1) return 'terminal';
  // A sugar exposes a branch attachment (R3) for the nucleobase; a
  // phosphate / linker only has chain-extending R-groups.
  const hasBranchR = rgroups.some((rg) => rg.label === 'R3') || rgroups.length >= 3;
  return hasBranchR ? NucleotideRole.SUGAR : NucleotideRole.PHOSPHATE;
}

// Assign a NucleotideRole to every monomer of every HELM RNA row by
// classifying each monomer from the library (see `classifyRnaRole`), then
// resolving single-R-group terminal modifiers to TERMINAL_5P / TERMINAL_3P
// based on the boundary they sit at.
//
// This is intentionally position-agnostic for the backbone: sugars, bases
// and phosphates are tagged by what they ARE, so an arbitrary number of
// phosphates / linkers may appear at the 5' end, between nucleotides, or
// anywhere else, and still be assembled as a connected chain.
//
// A row falls back to bases-only mode (undefined) only when it cannot be
// laid out as a linear chain: an unknown monomer, or a single-R-group
// terminal modifier appearing mid-chain (it has nothing to extend to).
export function buildRolesForHelmRna(
  monomerSequencesArray: ISeqMonomer[][],
  monomerLib: IMonomerLibBase,
  polymerType: PolymerType,
  chainStartsArray?: number[][],
): (NucleotideRole[] | undefined)[] {
  const out: (NucleotideRole[] | undefined)[] = new Array(monomerSequencesArray.length);
  for (let rowI = 0; rowI < monomerSequencesArray.length; ++rowI) {
    const row = monomerSequencesArray[rowI];
    if (row.length === 0) { out[rowI] = undefined; continue; }

    // Build roles one chain at a time. Terminal 5'/3' detection is scoped to
    // each chain's own boundaries — otherwise a chain end that sits in the
    // middle of the flat multi-chain row would be misread as a mid-chain
    // modifier and drop the whole row to bases-only mode.
    const starts = chainStartsArray?.[rowI];
    const chainBounds: number[] = (starts && starts.length > 0) ? starts : [0];
    const roles: NucleotideRole[] = new Array(row.length);
    let fellBack = false;
    for (let c = 0; c < chainBounds.length && !fellBack; ++c) {
      const s = chainBounds[c];
      const e = c + 1 < chainBounds.length ? chainBounds[c + 1] : row.length;
      const sliceRoles = buildRolesForChainSlice(row, s, e, monomerLib, polymerType);
      if (!sliceRoles) { fellBack = true; break; }
      for (let i = s; i < e; ++i) roles[i] = sliceRoles[i - s];
    }
    out[rowI] = fellBack ? undefined : roles;
  }
  return out;
}

// Classify the monomers of a single chain slice [start, end) into
// NucleotideRoles, resolving single-R-group terminal modifiers to
// TERMINAL_5P / TERMINAL_3P against this chain's own ends. Returns null when
// the slice cannot be laid out as a linear chain (unknown monomer, or a
// terminal modifier sitting mid-chain).
function buildRolesForChainSlice(
  row: ISeqMonomer[], start: number, end: number,
  monomerLib: IMonomerLibBase, polymerType: PolymerType,
): NucleotideRole[] | null {
  const len = end - start;
  if (len <= 0) return [];

  const classes: RnaRoleClass[] = new Array(len);
  for (let i = 0; i < len; ++i) {
    classes[i] = classifyRnaRole(monomerLib, polymerType, row[start + i].symbol);
    if (classes[i] === null) return null;
  }

  const roles: NucleotideRole[] = new Array(len);
  for (let i = 0; i < len; ++i) {
    if (classes[i] === 'terminal') {
      // A single-R-group modifier can only cap a chain end. HELM authors
      // sometimes put the "wrong" one at either boundary (e.g. GalNAc with
      // only R1 at the 5' end); getMolGraph swaps the rNodes so the real
      // R lands where the role expects it.
      if (i === 0) roles[i] = NucleotideRole.TERMINAL_5P;
      else if (i === len - 1) roles[i] = NucleotideRole.TERMINAL_3P;
      else return null;
    } else
      roles[i] = classes[i] as NucleotideRole;
  }
  return roles;
}

// Per-row disjoint-chain start positions (0-based indices into the flat
// monomer list). Reads the splitter's graphInfo; falls back to [0] for
// notations or rows without chain info.
function getChainStartsArray(seqCol: DG.Column<string>, seqHelper: ISeqHelper): number[][] {
  const sh = seqHelper.getSeqHandler(seqCol);
  const out: number[][] = new Array(seqCol.length);
  for (let rowI = 0; rowI < seqCol.length; ++rowI) {
    const seqSS = sh.getSplitted(rowI);
    const starts = seqSS.graphInfo?.disjointSeqStarts;
    out[rowI] = (starts && starts.length > 0) ? starts.slice() : [0];
  }
  return out;
}

/** Get jagged array of monomer symbols for the dataframe
 * @param {DG.Column} macroMolCol - Column with macro-molecules
 * @return {string[]} - Jagged array of monomer symbols for the dataframe */
export function getMonomerSequencesArray(macroMolCol: DG.Column<string>, seqHelper: ISeqHelper): ISeqMonomer[][] {
  const rowCount = macroMolCol.length;
  const result: ISeqMonomer[][] = new Array(rowCount);

  // split the string into monomers
  const sh = seqHelper.getSeqHandler(macroMolCol);

  let containsEmptyValues = false;
  const biotype: HelmType = sh.defaultBiotype;

  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const seqSS = sh.getSplitted(rowIdx);
    containsEmptyValues ||= seqSS.length === 0;
    result[rowIdx] = wu.count(0).take(seqSS.length).filter((posIdx) => !seqSS.isGap(posIdx))
      .map((posIdx) => {
        return {position: posIdx, biotype: biotype, symbol: seqSS.getCanonical(posIdx)} as ISeqMonomer;
      }).toArray();
  }

  if (containsEmptyValues)
    grok.shell.warning(`Some values in the "${macroMolCol.name}" column are empty`);

  return result;
}


/** Get a mapping of monomer symbols to MolGraph objects. Notice, the
 * transformation from molfile V2000 to V3000 takes place,
 * with the help of async function call from Chem (RdKit module)
 * @param {ISeqMonomer[][]} monomerSequencesArray - Jagged array of monomer symbols for the dataframe
 * @param {Array} rolesArray - Per-row NucleotideRole tags (undefined for legacy bases-only rows)
 * @param {PolymerType} polymerType - Polymer type
 * @param {ALPHABET} alphabet - Alphabet
 * @param {IMonomerLibBase} monomerLib - Monomer library
 * @param {RDModule} rdKitModule - RDKit module
 * @return {MonomerMolGraphMap} - Mapping of monomer symbols to MolGraph objects */
export function getMonomersDictFromLib(
  monomerSequencesArray: ISeqMonomer[][],
  rolesArray: (NucleotideRole[] | undefined)[],
  polymerType: PolymerType, alphabet: ALPHABET,
  monomerLib: IMonomerLibBase, rdKitModule: RDModule
): MonomerMolGraphMap {
  // todo: exception - no gaps, no empty string monomers
  const formattedMonomerLib = getFormattedMonomerLib(monomerLib, polymerType, alphabet);
  const monomersDict: MonomerMolGraphMap = {};

  const pointerToBranchAngle: NumberWrapper = {
    value: null
  };

  // Pre-load default sugar/phosphate when ANY row falls back to bases-only
  // mode (FASTA, separator, or HELM RNA rows the splitter could not split).
  // The default monomer adjusters depend on these being present in the dict.
  const anyBasesOnly = polymerType === HELM_POLYMER_TYPE.RNA &&
    rolesArray.some((r) => r === undefined);
  if (anyBasesOnly) {
    const symbols = (alphabet === ALPHABET.RNA) ?
      [C.RIBOSE, C.PHOSPHATE] : [C.DEOXYRIBOSE, C.PHOSPHATE];
    for (const sym of symbols) {
      addMonomerToDict(monomersDict, sym.symbol, formattedMonomerLib, rdKitModule, polymerType,
        pointerToBranchAngle, undefined);
    }
  }

  for (let rowI = 0; rowI < monomerSequencesArray.length; ++rowI) {
    const monomerSeq: ISeqMonomer[] = monomerSequencesArray[rowI];
    const roles = rolesArray[rowI];
    for (let posI = 0; posI < monomerSeq.length; ++posI) {
      const seqMonomer = monomerSeq[posI];
      const sym = seqMonomer.symbol;
      if (sym === '') continue; // Skip gap/empty monomer for MSA
      const role = roles ? roles[posI] : undefined;
      try {
        if (polymerType === HELM_POLYMER_TYPE.RNA && role === undefined &&
          sym.split(/\(|\)/).filter((e) => !!e).length === 3) {
          // legacy: a single token of the form sugar(base)phosphate slipped
          // through unsplit; index just the base.
          const nsym = sym.split(/\(|\)/)[1];
          addMonomerToDict(monomersDict, nsym, formattedMonomerLib,
            rdKitModule, polymerType, pointerToBranchAngle, NucleotideRole.BASE);
          if (monomersDict[polymerType]?.[nsym])
            monomersDict[polymerType][sym] = monomersDict[polymerType][nsym];
        } else {
          addMonomerToDict(monomersDict, sym, formattedMonomerLib,
            rdKitModule, polymerType, pointerToBranchAngle, role);
        }
      } catch (err: any) {
        const errTxt = err instanceof Error ? err.message : err.toString();
        const errStack = err instanceof Error ? err.stack : undefined;
        console.error(`bio lib: getMonomersDictFromLib() sym='${sym}', error:\n${errTxt}\n${errStack}`);
        const errMsg = `Can't get monomer '${sym}' from library: ${errTxt}`; // Text for Datagrok error baloon
        throw new Error(errMsg);
      }
    }
  }

  return monomersDict;
}

/** Adds MolGraph object for 'sym' to the monomers dict when necessary
 * @param {MonomerMolGraphMap} monomersDict - Monomers dictionary
 * @param {string} sym - Monomer symbol
 * @param {Map} formattedMonomerLib - Formatted monomer library
 * @param {any} moduleRdkit - RDKit module
 * @param {PolymerType} polymerType - Polymer type
 * @param {NumberWrapper} pointerToBranchAngle - Pointer to branch angle
 * @param {NucleotideRole} role - Optional role hint for RNA (sugar/base/phosphate) */
function addMonomerToDict(
  monomersDict: MonomerMolGraphMap, sym: string,
  formattedMonomerLib: Map<string, any>, moduleRdkit: any,
  polymerType: PolymerType, pointerToBranchAngle: NumberWrapper,
  role: NucleotideRole | undefined
): void {
  const symKey: LibMonomerKey = {polymerType, symbol: sym};
  if (!hasMolGraph(monomersDict, symKey)) {
    const monomerData: MolGraph | null =
      getMolGraph(sym, formattedMonomerLib, moduleRdkit, polymerType, pointerToBranchAngle, role);
    if (monomerData)
      setMolGraph(monomersDict, symKey, monomerData);
    else {
      // todo: handle exception
      throw new Error(`Monomer with symbol '${sym}' is absent the monomer library`);
    }
  }
}

/** Construct the MolGraph object for specified monomerSymbol: the associated
 * graph is adjusted in XY plane and filled with default R-groups
 * @param {string} monomerSymbol - Monomer symbol
 * @param {Map} formattedMonomerLib - Formatted monomer library
 * @param {any} moduleRdkit - RDKit module
 * @param {PolymerType} polymerType - Polymer type
 * @param {NumberWrapper} pointerToBranchAngle - Pointer to branch angle
 * @param {NucleotideRole} role - Optional explicit role (sugar/base/phosphate); when set
 *   it takes priority over the legacy symbol-based dispatch.
 * @return {MolGraph} - MolGraph object or null if monomerSymbol is absent in the library */
function getMolGraph(
  monomerSymbol: string, formattedMonomerLib: Map<string, any>,
  moduleRdkit: any, polymerType: PolymerType,
  pointerToBranchAngle: NumberWrapper,
  role: NucleotideRole | undefined
): MolGraph | null {
  if (!formattedMonomerLib.has(monomerSymbol))
    return null;
  else {
    const libObject = formattedMonomerLib.get(monomerSymbol);
    const capGroups = parseCapGroups(libObject[HELM_FIELDS.RGROUPS]);
    const capGroupIdxMap = parseCapGroupIdxMap(libObject[HELM_FIELDS.MOLFILE]);
    //in case molfile is molV2000 need to convert it to molV3000
    const molfileV3K = libObject[HELM_FIELDS.MOLFILE].includes('V3000') ? libObject[HELM_FIELDS.MOLFILE] :
      convertMolfileToV3K(removeRGroupLines(libObject[HELM_FIELDS.MOLFILE]), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);

    const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
    const bonds = parseBondBlock(molfileV3K, counts.bondCount);
    const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);
    const steabs = getAbsStereocenters(molfileV3K);

    const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta, stereoAtoms: steabs};

    if (polymerType === HELM_POLYMER_TYPE.PEPTIDE)
      adjustPeptideMonomerGraph(monomerGraph);
    else { // nucleotides
      // Prefer explicit role (HELM-RNA triples mode); fall back to symbol
      // matching against the canonical r/d/p for the legacy bases-only path.
      const isSugar = role === NucleotideRole.SUGAR ||
        (role === undefined && (monomerSymbol === C.RIBOSE.symbol || monomerSymbol === C.DEOXYRIBOSE.symbol));
      const isPhosphate = role === NucleotideRole.PHOSPHATE ||
        (role === undefined && monomerSymbol === C.PHOSPHATE.symbol);
      const isTerminal = role === NucleotideRole.TERMINAL_5P || role === NucleotideRole.TERMINAL_3P;
      if (isTerminal) {
        // Terminal monomers have ONE real R-group from the lib; the other
        // is a pseudo synthesized by getMonomerMetadata and points at a
        // real atom of the monomer body. The TERMINAL_5P/TERMINAL_3P
        // logic in setShiftsAndTerminalNodes assumes:
        //   TERMINAL_5P → real R is at rNodes[1] (chain-extending right).
        //   TERMINAL_3P → real R is at rNodes[0] (chain-incoming left).
        // When HELM places a single-R monomer at the "wrong" boundary
        // (e.g. Bio/GalNAc with only R1 at row[0], or Chol with only R2
        // at row[last]), the real R sits at the OTHER index. Swap the
        // rNodes / terminalNodes pair so the existing role logic works.
        const realRLabel = libObject[HELM_FIELDS.RGROUPS]?.[0]?.label;
        const realRIdx = realRLabel === 'R1' ? 0 : (realRLabel === 'R2' ? 1 : -1);
        const expectedRealRIdx = role === NucleotideRole.TERMINAL_5P ? 1 : 0;
        if (realRIdx !== -1 && realRIdx !== expectedRealRIdx) {
          [monomerGraph.meta.rNodes[0], monomerGraph.meta.rNodes[1]] =
            [monomerGraph.meta.rNodes[1], monomerGraph.meta.rNodes[0]];
          [monomerGraph.meta.terminalNodes[0], monomerGraph.meta.terminalNodes[1]] =
            [monomerGraph.meta.terminalNodes[1], monomerGraph.meta.terminalNodes[0]];
        }
        // Peptide-style adjustment correctly orients a 2-R monomer.
        adjustPeptideMonomerGraph(monomerGraph);
      } else if (isSugar)
        adjustSugarMonomerGraph(monomerGraph, pointerToBranchAngle);
      else if (isPhosphate)
        adjustPhosphateMonomerGraph(monomerGraph);
      else
        adjustBaseMonomerGraph(monomerGraph, pointerToBranchAngle);
    }

    setShiftsAndTerminalNodes(polymerType, monomerGraph, monomerSymbol, role);
    // todo: restore after debugging
    removeHydrogen(monomerGraph);
    replaceWrongfulRGroups(monomerGraph);

    removeRgroupKwargs(monomerGraph);

    return monomerGraph;
  }
}

function getAbsStereocenters(molfileV3K: string): number[] {
  let collection: number[] = [];
  let indexCollection = molfileV3K.indexOf('M  V30 MDLV30/STEABS ATOMS=('); // V3000 index for collections

  while (indexCollection !== -1) {
    indexCollection += 28;
    const collectionEnd = molfileV3K.indexOf(')', indexCollection);
    //using slice(1) because first digit in brackets is a number of stereo atoms
    collection = collection.concat(molfileV3K.substring(indexCollection, collectionEnd).split(' ')
      .slice(1).map((it) => parseInt(it)));
    indexCollection = collectionEnd;
    indexCollection = molfileV3K.indexOf('M  V30 MDLV30/STEABS ATOMS=(', indexCollection);
  }
  return collection;
}

function setShiftsAndTerminalNodes(
  polymerType: PolymerType, monomerGraph: MolGraph, monomerSymbol: string,
  role: NucleotideRole | undefined
): void {
  // remove the 'rightmost' chain-extending r-group node in the backbone
  if (polymerType === HELM_POLYMER_TYPE.PEPTIDE) {
    setShifts(monomerGraph, polymerType);
    const removedR2Atom = removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
    if (removedR2Atom?.removedAtom)
      monomerGraph.terminalR2Atom = removedR2Atom.removedAtom; // store the removed R2 Atom type, if any
  } else if (role === NucleotideRole.TERMINAL_5P) {
    // 5'-end terminal modifier (e.g. Chol): R2 is the REAL R-group (chain
    // attach point); R1 is a pseudo synthesized by getMonomerMetadata that
    // points at a real atom of the monomer — DO NOT remove it. Remove the
    // R2 placeholder so the chain bond can attach at terminalNodes[1].
    setShifts(monomerGraph, polymerType);
    removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
  } else if (role === NucleotideRole.TERMINAL_3P) {
    // 3'-end terminal modifier (e.g. GalNAc): R1 is the REAL R-group
    // placeholder. After substituteCapGroups it carries the cap atom type
    // (e.g. 'O' from "OH" or 'H' from "H") — that atom represents the
    // FREE-form residue at this position and MUST be removed when the
    // monomer is in a chain (otherwise we leave a stray OH/H bonded to the
    // chain-attach carbon). After removal, terminalNodes[0] still points
    // at the correct neighbor (the chain-attach atom).
    // R2 is the PSEUDO synthesized by getMonomerMetadata; it points at a
    // real atom of the monomer — DO NOT remove it.
    removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[0]);
  } else { // nucleotides
    // Prefer explicit role (HELM RNA triples mode); fall back to symbol
    // matching against canonical r/d/p for the legacy bases-only path.
    const isSugar = role === NucleotideRole.SUGAR ||
      (role === undefined && (monomerSymbol === C.RIBOSE.symbol || monomerSymbol === C.DEOXYRIBOSE.symbol));
    const isPhosphate = role === NucleotideRole.PHOSPHATE ||
      (role === undefined && monomerSymbol === C.PHOSPHATE.symbol);
    if (isSugar) {
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
    } else if (isPhosphate) {
      // Chain bond comes IN at the R1 cap atom. For canonical phosphates
      // (R1 cap "OH" → cap atom is O) this is the bridging oxygen between
      // the previous sugar's C3' and the central atom — exactly the
      // chemistry we want.
      //
      // Variants with R1 cap "H" (sp, en, ...) substitute to H instead;
      // the previous sugar (in the standard flow) has already given up its
      // own 3'-O on the assumption that the linker will bring a bridging
      // atom, and removeHydrogen would later strip the H, leaving the
      // chain bond attached directly to the central atom (e.g., C3'-P
      // instead of C3'-O-P). Promote the cap H to an O so the bridging
      // atom exists in the assembled chain regardless of how the library
      // happened to spell the cap. After the promotion both branches go
      // through the same terminalNodes[0]=cap-atom code path.
      const r1Idx = monomerGraph.meta.rNodes[0] - 1;
      const r1AtomType = monomerGraph.atoms.atomTypes[r1Idx];
      if (r1AtomType?.toUpperCase() === C.HYDROGEN)
        monomerGraph.atoms.atomTypes[r1Idx] = C.OXYGEN;
      monomerGraph.meta.terminalNodes[0] = monomerGraph.meta.rNodes[0];
      shiftCoordinates(
        monomerGraph,
        -monomerGraph.atoms.x[monomerGraph.meta.terminalNodes[0] - 1],
        -monomerGraph.atoms.y[monomerGraph.meta.terminalNodes[0] - 1]
      );
      setShifts(monomerGraph, polymerType);
      removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[1]);
    } else { // nucleobases
      // removeNodeAndBonds(monomerGraph, monomerGraph.meta.rNodes[0]);
    }
  }
}


/**
 * Get monomer metadata object
 * @param {Atoms} atoms - Atoms object
 * @param {Bonds} bonds - Bonds object
 * @param {string[]} capGroups - Cap groups
 * @param {Map<number, number>} capGroupIdxMap - Cap group index map
 * @return {MonomerMetadata}*/
function getMonomerMetadata(atoms: Atoms, bonds: Bonds, capGroups?: string[], capGroupIdxMap?: Map<number, number>
): MonomerMetadata {
  const meta: MonomerMetadata = {
    backboneShift: null,
    branchShift: null,
    terminalNodes: [],
    rNodes: [],
  };

  // R1 and R2 are required, but we might have cases where one of those is absent,
  // e.g. in starting monomers or ending monomers. in those cases,
  // mark the missing r-node as the furthest atom from non-missing one
  if (capGroups && capGroupIdxMap && [1, 2].some((i) => !Array.from(capGroupIdxMap.values()).find((j) => j == i))) {
    const missingRGroup = [1, 2].find((i) => !Array.from(capGroupIdxMap.values()).find((j) => j == i))!;
    const existingRgroup = [1, 2].find((i) => Array.from(capGroupIdxMap.values()).find((j) => j == i))!;
    const existingRgroupIdx = Array.from(capGroupIdxMap.keys()).find((i) => capGroupIdxMap.get(i) === existingRgroup)! - 1; // -1 because molfile indexing starts from 1
    const existingRgroupX = atoms.x[existingRgroupIdx];
    const existingRgroupY = atoms.y[existingRgroupIdx];
    const atomBondedToExistingRgroupIdx = bonds.atomPairs.find((pair) => pair.includes(existingRgroupIdx + 1))!.
      find((i) => i !== existingRgroupIdx + 1)! - 1; // -1 because molfile indexing starts from 1

    let furthestAtomIdx = atoms.x.reduce((furthestIdx, x, idx) => {
      if (idx === existingRgroupIdx)
        return furthestIdx; // skip existing r-group
      // also skip the atom bonded to the existing r-group
      if (idx === atomBondedToExistingRgroupIdx)
        return furthestIdx;
      const y = atoms.y[idx];
      const dist = Math.sqrt((x - existingRgroupX) ** 2 + (y - existingRgroupY) ** 2);
      return dist > Math.sqrt((atoms.x[furthestIdx] - existingRgroupX) ** 2 +
        (atoms.y[furthestIdx] - existingRgroupY) ** 2) ? idx : furthestIdx;
    }, -1);
    if (furthestAtomIdx === -1) {
      furthestAtomIdx = atoms.x.length;
      // add pseudo-hydrogen atom, very relevant for stuff like NH2-R
      atoms.x = new Float32Array([...atoms.x, -existingRgroupX]);
      atoms.y = new Float32Array([...atoms.y, -existingRgroupY]);
      atoms.atomTypes = [...atoms.atomTypes, 'H'];
      atoms.kwargs = [...atoms.kwargs, ''];
      bonds.atomPairs.push([atomBondedToExistingRgroupIdx + 1, furthestAtomIdx + 1]); // +1 because molfile indexing starts from 1
      bonds.bondTypes = new Uint32Array([...bonds.bondTypes, 1]); // single bond
      bonds.kwargs.set(bonds.atomPairs.length - 1, ''); // empty kwargs for the new bond
    }


    capGroupIdxMap.set(furthestAtomIdx + 1, missingRGroup); // +1 because molfile indexing starts from 1
    // finaly splice the capGroups array in correct place
    if (missingRGroup === 1)
      capGroups.unshift(atoms.atomTypes[furthestAtomIdx]); // add to the beginning
    else if (missingRGroup === 2)
      capGroups.splice(1, 0, atoms.atomTypes[furthestAtomIdx]); // add to the second position
    else
      throw new Error(`Unexpected missing R-group: ${missingRGroup}`);
  }

  substituteCapGroups(atoms, capGroups!, capGroupIdxMap!);
  setRNodes(capGroupIdxMap!, meta);

  setTerminalNodes(bonds, meta);
  return meta;
}

/** Parse element symbols for R-groups from the HELM monomer library R-group field
 * @param {any[]} rGroupObjList - R-group object list
 * @return {Map<number, number>} - Cap group index map*/
export function parseCapGroups(rGroupObjList: any[]): string[] {
  // specifically for HELMCoreLibrary
  // considered only monoatomic rgroups
  // supposing that elements in rGroupObjList are sorted w.r.t. the rgroups idx
  const capGroupsArray: string[] = [];
  for (const obj of rGroupObjList) {
    let capGroup: string = obj[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES];

    // in some cases the smiles field is written with uppercase
    if (!capGroup)
      capGroup = obj[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
    // some nice person can set the smiles field in r groups to have OH, which is incorrect smiles
    capGroup = capGroup.replace(/(\[|\]|\*|:|\d)/g, '').replace('OH', 'O').replace('Oh', 'O');
    capGroupsArray.push(capGroup);
  }
  return capGroupsArray;
}


/** Substitute the cap group elements instead of R#
 * @param {Atoms} atoms - Atoms object
 * @param {string[]} capGroups - Cap groups
 * @param {Map<number, number>} capGroupIdxMap - Cap group index map*/
function substituteCapGroups(
  atoms: Atoms, capGroups: string[], capGroupIdxMap: Map<number, number>
): void {
  for (const [node, capIdx] of capGroupIdxMap)
    atoms.atomTypes[node - 1] = capGroups[capIdx - 1]; // -1 because molfile indexing starts from 1
}

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

/** Sets shifts in 'meta' attribute of MolGraph
 * @param {MolGraph} molGraph - MolGraph object
 * @param {PolymerType} polymerType - Polymer type*/
function setShifts(molGraph: MolGraph, polymerType: PolymerType): void {
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

/** Returns the pair [xShift, yShift] for specified node indices
 * @param {MolGraph} molGraph - MolGraph object
 * @param {number} rightNodeIdx - Right node index
 * @param {number} leftNodeIdx - Left node index
 * @return {number[]} - Shift between nodes*/
function getShiftBetweenNodes(
  molGraph: MolGraph, rightNodeIdx: number, leftNodeIdx: number
): number[] {
  return [
    keepPrecision(molGraph.atoms.x[rightNodeIdx] - molGraph.atoms.x[leftNodeIdx]),
    keepPrecision(molGraph.atoms.y[rightNodeIdx] - molGraph.atoms.y[leftNodeIdx]),
  ];
}


/** Helper function necessary to build a correct V3000 molfile out of V2000 with
 * specified r-groups
 * @param {string} molfileV2K - V2000 molfile
 * @return {string} - V2000 molfile without R-group lines*/
function removeRGroupLines(molfileV2K: string): string {
  let begin = molfileV2K.indexOf(C.V2K_A_LINE, 0);
  if (begin === -1)
    begin = molfileV2K.indexOf(C.V2K_RGP_LINE);
  const end = molfileV2K.indexOf(C.V3K_END, begin);
  return molfileV2K.substring(0, begin) + molfileV2K.substring(end);
}

/** V2000 to V3000 converter
 * @param {string} molfileV2K - V2000 molfile
 * @param {any} moduleRdkit - RDKit module
 * @return {string} - V3000 molfile*/
export function convertMolfileToV3K(molfileV2K: string, moduleRdkit: any): string {
  // The standard Chem converter is not used here because it relies on creation of moduleRdkit on each iteration
  const molObj = moduleRdkit.get_mol(fixV2000MolfileRAtomLines(molfileV2K));
  const molfileV3K = molObj.get_v3Kmolblock();
  molObj.delete();
  return molfileV3K;
}

/** fixes the r lines of the molblock. rdkit sometimes gives v2000 with weird numbers in r atom lines,
 * that can be interpreted as isotopes or cause downstream weirdness
 * @param {string} molfileV2K
 * @return {string} */
export function fixV2000MolfileRAtomLines(molfileV2K: string): string {
  const lines = molfileV2K.split('\n');
  const fixedLines = lines.map((line) => {
    const rIndex1 = line.indexOf(' R#  ');
    const rIndex2 = line.indexOf(' R   ');
    if (rIndex1 === -1 && rIndex2 === -1)
      return line;
    const rIndex = rIndex1 !== -1 ? rIndex1 : rIndex2;
    // all numbers after R should be 0, otherwise they get interpreted as isotopes or smth weird
    const lineArray = line.split('');
    for (let i = rIndex + 5; i < lineArray.length; i++) {
      if (lineArray[i] === ' ')
        continue;
      lineArray[i] = '0';
    }
    return lineArray.join('');
  });
  return fixedLines.join('\n');
}

/** Parse V3000 bond block and construct the Bonds object
 * @param {string} molfileV3K - V3000 molfile
 * @param {number} bondCount - Number of bonds
 * @return {Bonds} - Bonds object*/
export function parseBondBlock(molfileV3K: string, bondCount: number): Bonds {
  const bondTypes: Uint32Array = new Uint32Array(bondCount);
  const atomPairs: number[][] = new Array(bondCount);
  const bondConfiguration = new Map<number, number>();
  const kwargs = new Map<number, string>();

  let begin = molfileV3K.indexOf(C.V3K_BEGIN_BOND_BLOCK);
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;
  for (let i = 0; i < bondCount; ++i) {
    // parse bond type and atom pair
    const parsedValues: number[] = new Array(3);
    begin = molfileV3K.indexOf(C.V3K_BEGIN_DATA_LINE, end) + C.V3K_IDX_SHIFT;
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
    let beginCfg = lineRemainder.indexOf(C.V3K_BOND_CONFIG);
    if (beginCfg !== -1) {
      beginCfg = lineRemainder.indexOf('=', beginCfg) + 1;
      let endCfg = lineRemainder.indexOf(' ', beginCfg);
      if (endCfg === -1)
        endCfg = lineRemainder.length;
      const bondConfig = parseInt(lineRemainder.slice(beginCfg, endCfg));
      bondConfiguration.set(i, bondConfig);
      const removedSubstring = C.V3K_BOND_CONFIG + bondConfig.toString();
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

/** Constructs mapping of r-group nodes to default capGroups, all numeration starting from 1.
 * According to https://pubs.acs.org/doi/10.1021/ci3001925, R1 and R2 are the chain extending attachment points,
 * while R3 is the branching attachment point.
 * @param {string} molfile - V2000 or V3000 molfile
 * @return {Map<number, number>} - Map of r-group nodes to default capGroups*/
export function parseCapGroupIdxMap(molfile: string): Map<number, number> {
  const isMolV3K = molfile.includes('V3000');
  return isMolV3K ? parseCapGroupIdxMapV3K(molfile) : parseCapGroupIdxMapV2K(molfile);
}

export function parseCapGroupIdxMapV2K(molfileV2K: string): Map<number, number> {
  const capGroupIdxMap = new Map<number, number>();

  // parse A-lines (RNA)
  let begin = molfileV2K.indexOf(C.V2K_A_LINE, 0);
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

    begin = molfileV2K.indexOf(C.V2K_A_LINE, end);
  }

  // parse RGP lines (may be more than one in RNA monomers)
  begin = molfileV2K.indexOf(C.V2K_RGP_LINE, 0);
  end = molfileV2K.indexOf('\n', begin);
  while (begin !== -1) {
    begin += C.V2K_RGP_SHIFT;
    end = molfileV2K.indexOf('\n', begin);
    const rgpStringParsed = molfileV2K.substring(begin, end)
      .replaceAll(/\s+/g, ' ')
      .split(' ');
    const rgpIndicesArray = rgpStringParsed.map((el) => parseInt(el))
      .slice(1); // slice from 1 because the 1st value is the number of pairs in the line
    for (let i = 0; i < rgpIndicesArray.length; i += 2) {
      // there may be conflicting cap group definitions, like 3-O-Methylribose (2,5 connectivity) in HELMCoreLibrary
      if (capGroupIdxMap.has(rgpIndicesArray[i]) && capGroupIdxMap.get(rgpIndicesArray[i]) !== rgpIndicesArray[i + 1])
        throw new Error(`r-group index ${rgpIndicesArray[i]} has already been added with a different value`);
      else
        capGroupIdxMap.set(rgpIndicesArray[i], rgpIndicesArray[i + 1]);
    }
    begin = molfileV2K.indexOf(C.V2K_RGP_LINE, end);
  }
  return capGroupIdxMap;
}

export function parseCapGroupIdxMapV3K(molfileV3K: string): Map<number, number> {
  const capGroupIdxMap = new Map<number, number>();
  const regex = /M  V30 (\d+) R#.+RGROUPS=\((\d+) (\d+)\).*/gm;
  let res;
  while ((res = regex.exec(molfileV3K)) !== null) {
    // This is necessary to avoid infinite loops with zero-width matches
    if (res.index === regex.lastIndex)
      regex.lastIndex++;

    // index 1 matches array is the atom number, index 3 is R group number
    capGroupIdxMap.set(parseInt(res[1]), parseInt(res[3]));
  }

  return capGroupIdxMap;
}

export function parseAtomAndBondCounts(molfileV3K: string): { atomCount: number, bondCount: number } {
  molfileV3K = molfileV3K.replaceAll('\r', ''); // to handle old and new sdf standards

  // parse atom count
  let begin = molfileV3K.indexOf(C.V3K_BEGIN_COUNTS_LINE) + C.V3K_COUNTS_SHIFT;
  let end = molfileV3K.indexOf(' ', begin + 1);
  const numOfAtoms = parseInt(molfileV3K.substring(begin, end));

  // parse bond count
  begin = end + 1;
  end = molfileV3K.indexOf(' ', begin + 1);
  const numOfBonds = parseInt(molfileV3K.substring(begin, end));

  return {atomCount: numOfAtoms, bondCount: numOfBonds};
}

/** Parse V3000 atom block and return Atoms object. NOTICE: only atomTypes, x, y
 * and kwargs fields are set in the return value, with other fields dummy initialized.
 * @param {string} molfileV3K - V3000 molfile
 * @param {number} atomCount - number of atoms in the molecule
 * @return {Atoms} - Atoms object */
function parseAtomBlock(molfileV3K: string, atomCount: number): Atoms {
  const atomTypes: string[] = new Array(atomCount);
  const x: Float32Array = new Float32Array(atomCount);
  const y: Float32Array = new Float32Array(atomCount);
  const kwargs: string[] = new Array(atomCount);

  let begin = molfileV3K.indexOf(C.V3K_BEGIN_ATOM_BLOCK); // V3000 atoms block
  begin = molfileV3K.indexOf('\n', begin);
  let end = begin;

  for (let i = 0; i < atomCount; i++) {
    begin = molfileV3K.indexOf(C.V3K_BEGIN_DATA_LINE, begin) + C.V3K_IDX_SHIFT;
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

/** Remove hydrogen nodes
 * @param {MolGraph} monomerGraph - monomer graph*/
function removeHydrogen(monomerGraph: MolGraph): void {
  let i = 0;
  while (i < monomerGraph.atoms.atomTypes.length) {
    if (monomerGraph.atoms.atomTypes[i] === C.HYDROGEN) {
      removeNodeAndBonds(monomerGraph, i + 1); // i + 1 because molfile node indexing starts from 1
      --i;
      // monomerGraph.atoms.atomTypes[i] = 'Li';
    }
    ++i;
  }
}

/** Replaces wrongly registered r-groups with corrections. for example, some monomers have [*:1][OH] or smth similar
 * @param {MolGraph} monomerGraph - monomer graph*/
function replaceWrongfulRGroups(monomerGraph: MolGraph): void {
  let i = 0;
  while (i < monomerGraph.atoms.atomTypes.length) {
    if (monomerGraph.atoms.atomTypes[i]?.toLowerCase() === 'oh')
      monomerGraph.atoms.atomTypes[i] = 'O';
    if (monomerGraph.atoms.atomTypes[i] === '?')
      monomerGraph.atoms.atomTypes[i] = 'H'; //replace ? with H
    ++i;
  }
}

/** After processing and removing r-groups, the KWARGS property might still contain RGROUPS lines.
 * this needs to be removed
 * @param {MolGraph} monomerGraph - monomer graph*/
function removeRgroupKwargs(monomerGraph: MolGraph): void {
  const molv3kRGroupKwargLine = ' RGROUPS=(1 1)';
  for (let i = 0; i < (monomerGraph.atoms.kwargs?.length ?? 0); i++) {
    const kwarg = monomerGraph.atoms.kwargs[i];
    if (kwarg && kwarg.includes(molv3kRGroupKwargLine))
      monomerGraph.atoms.kwargs[i] = kwarg.replace(molv3kRGroupKwargLine, ''); // remove the line
  }
}

/** Remove node 'removedNode' and the associated bonds. Notice, numeration of
 * nodes in molfiles starts from 1, not 0
 * @param {MolGraph} monomerGraph - monomer graph
 * @param {number} removedNode - node to be removed
 * @return {{removedAtom: string} | undefined} - removed atom type, if any
 * */
function removeNodeAndBonds(monomerGraph: MolGraph, removedNode?: number): {removedAtom: string} | undefined {
  if (typeof removedNode !== 'undefined') {
    const removedNodeIdx = removedNode - 1;
    const atoms = monomerGraph.atoms;
    const bonds = monomerGraph.bonds;
    const meta = monomerGraph.meta;

    // remove the node from atoms
    const removedAtomType = atoms.atomTypes.splice(removedNodeIdx, 1)[0];
    atoms.x = spliceTypedArray<Float32Array>(Float32Array, atoms.x, removedNodeIdx, 1);
    atoms.y = spliceTypedArray<Float32Array>(Float32Array, atoms.y, removedNodeIdx, 1);
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
        bonds.bondTypes = spliceTypedArray<Uint32Array>(Uint32Array, bonds.bondTypes, i, 1);
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
    return removedAtomType ? {removedAtom: removedAtomType} : undefined;
  }
}

// todo: rewrite the following two functions using templates,
function spliceTypedArray<T extends ITypedArray>(
  TConstructor: { new(length: number): T; }, typedArray: T, start: number, count: number
) {
  const result = new TConstructor(typedArray.length - count);
  let i = 0;
  let k = 0;
  while (i < typedArray.length) {
    if (i === start)
      i += count;
    result[k] = typedArray[i];
    ++k;
    ++i;
  }
  return result;
}


/** Adjust the peptide MolGraph to default/standardized position
 * @param {MolGraph} monomer - monomer graph*/
function adjustPeptideMonomerGraph(monomer: MolGraph): void {
  const centeredNode = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  const rotatedNode = monomer.meta.rNodes[0] - 1;
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // place nodeOne at origin
  shiftCoordinates(monomer, -x[centeredNode], -y[centeredNode]);

  // angle is measured between OY and the rotated node
  const angle = findAngleWithOY(x[rotatedNode], y[rotatedNode]);

  // rotate the centered graph, so that 'nodeTwo' ends up on the positive ray of OY
  rotateCenteredGraph(monomer.atoms, -angle);

  if (x[monomer.meta.rNodes[1] - 1] < 0)
    flipMonomerAroundOY(monomer);

  const doubleBondedOxygen = findDoubleBondedCarbonylOxygen(monomer);
  if (doubleBondedOxygen != null) {
    // flip carboxyl and R if necessary
    flipCarboxylAndRadical(monomer, doubleBondedOxygen);

    // flip hydroxyl group with double-bound O inside carboxyl group if necessary
    flipHydroxilGroup(monomer, doubleBondedOxygen);
  }
}

function adjustPhosphateMonomerGraph(monomer: MolGraph): void {
  // Orient the linker so that the atoms connected to its R-groups lie on a
  // horizontal line through the origin:
  //   - terminalNodes[0] (atom connected to R1) sits at the origin
  //   - terminalNodes[1] (atom connected to R2) sits on the +X ray
  //   - rNodes[0] (R1 placeholder) ends up on the left (-X side)
  //     and rNodes[1] (R2 placeholder) on the right (+X side)
  //
  // For a "small" linker (e.g. canonical phosphate) where R1 and R2 attach
  // to the SAME central atom, terminalNodes[0] == terminalNodes[1] and we
  // can't form a horizontal segment from a single atom — fall back to
  // orienting the rNodes[0]→rNodes[1] vector along +X instead.
  //
  // No atom-type assumptions: this works for phosphate, phosphorothioate,
  // sulfonate, carbonate or any other linker the library may carry as long
  // as the R-groups are present.
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  if (monomer.meta.terminalNodes.length === 0) return;

  const t0Idx = monomer.meta.terminalNodes[0] - 1;
  const t1Idx = monomer.meta.terminalNodes.length > 1 ? monomer.meta.terminalNodes[1] - 1 : -1;

  // 1. Center on terminalNodes[0]
  shiftCoordinates(monomer, -x[t0Idx], -y[t0Idx]);

  // 2. Pick the rotation reference. Prefer terminalNodes[1] when distinct
  //    from terminalNodes[0]; otherwise fall back to rNodes[1] / rNodes[0].
  let refX: number = 0; let refY: number = 0;
  if (t1Idx >= 0 && t1Idx !== t0Idx) {
    refX = x[t1Idx]; refY = y[t1Idx];
  } else if (monomer.meta.rNodes.length > 1 && monomer.meta.rNodes[1] - 1 >= 0) {
    refX = x[monomer.meta.rNodes[1] - 1]; refY = y[monomer.meta.rNodes[1] - 1];
  } else if (monomer.meta.rNodes.length > 0 && monomer.meta.rNodes[0] - 1 >= 0) {
    // single-R linker — treat the R1 placeholder direction as the "incoming"
    // (left) side, so flip its vector to land R2-ish geometry on +X.
    refX = -x[monomer.meta.rNodes[0] - 1]; refY = -y[monomer.meta.rNodes[0] - 1];
  }

  if (refX !== 0 || refY !== 0) {
    const angleFromOY = findAngleWithOY(refX, refY);
    rotateCenteredGraph(monomer.atoms, -angleFromOY - Math.PI / 2);
  }

  // 3. Make sure R1 placeholder is on the left of R2 placeholder.
  if (monomer.meta.rNodes.length >= 2) {
    const r1Idx = monomer.meta.rNodes[0] - 1;
    const r2Idx = monomer.meta.rNodes[1] - 1;
    if (r1Idx >= 0 && r2Idx >= 0 && x[r1Idx] > x[r2Idx])
      flipMonomerAroundOY(monomer);
  }
}

function adjustSugarMonomerGraph(monomer: MolGraph, pointerToBranchAngle: NumberWrapper): void {
  // Sugars connect via R1 (to the previous linker), R2 (to the next linker),
  // and R3 (to the branch — the nucleobase). Goal:
  //   - terminalNodes[0] (atom connected to R1) at the origin
  //   - terminalNodes[1] (atom connected to R2) on the +X ray
  //     → atoms connected to R1/R2 are on a horizontal line
  //   - rNodes[2] (R3 placeholder) lies in the upper half plane
  //     → R3 (and the nucleobase that will attach there) points up
  //
  // The previous algorithm rotated rNodes[1] (the R2 placeholder atom) onto
  // +X axis, which only worked when the molfile happened to place R2
  // directly opposite terminalNodes[0]. For non-canonical sugars (modified
  // riboses, sugar surrogates) this produced skewed layouts. Rotating
  // terminalNodes[1] is robust to whatever atom types or bond directions
  // the library uses.
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  if (monomer.meta.terminalNodes.length === 0) return;

  const t0Idx = monomer.meta.terminalNodes[0] - 1;
  const t1Idx = monomer.meta.terminalNodes.length > 1 ? monomer.meta.terminalNodes[1] - 1 : -1;

  // 1. Center on terminalNodes[0]
  shiftCoordinates(monomer, -x[t0Idx], -y[t0Idx]);

  // 2. Rotate so terminalNodes[1] (the atom bonded to R2) lies on +X.
  //    Skip when terminalNodes[1] is missing or coincides with t0
  //    (degenerate sugar — fall back to rNodes[1] for orientation).
  let refX: number = 0; let refY: number = 0;
  if (t1Idx >= 0 && t1Idx !== t0Idx) {
    refX = x[t1Idx]; refY = y[t1Idx];
  } else if (monomer.meta.rNodes.length > 1 && monomer.meta.rNodes[1] - 1 >= 0) {
    refX = x[monomer.meta.rNodes[1] - 1]; refY = y[monomer.meta.rNodes[1] - 1];
  }
  if (refX !== 0 || refY !== 0) {
    const angleFromOY = findAngleWithOY(refX, refY);
    rotateCenteredGraph(monomer.atoms, -angleFromOY - Math.PI / 2);
  }

  // 3. Flip around OX so R3 points up. Prefer the R3 placeholder itself,
  //    fall back to terminalNodes[2] (atom connected to R3) when missing.
  let r3y: number | null = null;
  if (monomer.meta.rNodes.length > 2 && monomer.meta.rNodes[2] - 1 >= 0)
    r3y = y[monomer.meta.rNodes[2] - 1];
  else if (monomer.meta.terminalNodes.length > 2 && monomer.meta.terminalNodes[2] - 1 >= 0)
    r3y = y[monomer.meta.terminalNodes[2] - 1];
  if (r3y !== null && r3y < 0)
    flipMonomerAroundOX(monomer);

  // 4. Abnormal-sugar override: when the sugar's R3-attachment atom
  //    (terminalNodes[2]) is NOT the topmost atom of the molecule, the
  //    natural branch direction would push the base sideways or even down
  //    (LNA, with its 2,4-O-CH2 bridge sitting above C1', is the canonical
  //    example). In that case reposition the R3 placeholder directly above
  //    the topmost real atom, so the branch shift puts the base's
  //    terminalNodes[0] above the entire sugar with a vertical bond going
  //    straight up. This is purely a 2D-depiction adjustment — it preserves
  //    connectivity and only changes where the base atoms land in space.
  repositionR3IfAbnormal(monomer);

  // 5. Compute branch angle for the base alignment step. Default to 0
  //    (straight up) when the sugar lacks the branch attachment — the base
  //    placement code only runs when a base actually follows. The override
  //    above, when triggered, places R3 directly above terminalNodes[2]
  //    so this naturally evaluates to ~0.
  if (monomer.meta.rNodes.length > 2 && monomer.meta.terminalNodes.length > 2 &&
      monomer.meta.rNodes[2] - 1 >= 0 && monomer.meta.terminalNodes[2] - 1 >= 0)
    pointerToBranchAngle.value = getAngleBetweenSugarBranchAndOY(monomer);
  else
    pointerToBranchAngle.value = 0;

  // 6. Re-anchor terminalNodes[0] at the origin (rotation/flip preserve it,
  //    but keep this as a safety net for floating-point drift).
  const t0IdxFinal = monomer.meta.terminalNodes[0] - 1;
  shiftCoordinates(monomer, -x[t0IdxFinal], -y[t0IdxFinal]);
}

/** When the sugar's R3-attachment atom (terminalNodes[2]) is not the
 * topmost atom, move the R3 placeholder so it sits directly above the
 * sugar's topmost real atom (with one bond length of vertical clearance).
 * Sugars like LNA — with a methylene/oxygen bridge spanning over the
 * face that carries C1' — fall into this branch.
 *
 * Only the R3 placeholder atom moves; the rest of the sugar (and the
 * bond from terminalNodes[2] to rNodes[2]) is preserved. The bond ends
 * up elongated and vertical, which is fine for SMILES (length is
 * irrelevant to connectivity) and gives the chain a clean visual layout
 * with the base above any bridging atoms. */
function repositionR3IfAbnormal(monomer: MolGraph): void {
  if (monomer.meta.rNodes.length <= 2 || monomer.meta.terminalNodes.length <= 2)
    return;
  const r3 = monomer.meta.rNodes[2] - 1;
  const t2 = monomer.meta.terminalNodes[2] - 1;
  if (r3 < 0 || t2 < 0) return;

  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  // Topmost Y among atoms that will remain in the assembly: exclude the
  // R-group placeholder positions (rNodes[*]), since those are either
  // removed or reassigned during setShiftsAndTerminalNodes.
  const rGroupSet = new Set<number>();
  for (const r of monomer.meta.rNodes) {
    if (r >= 1) rGroupSet.add(r - 1);
  }
  let topY = -Infinity;
  for (let i = 0; i < y.length; ++i) {
    if (!rGroupSet.has(i) && y[i] > topY) topY = y[i];
  }
  if (!Number.isFinite(topY)) return;

  // Tolerance: don't trigger for the normal case where terminalNodes[2]
  // is at (or near) the top — small floating-point or geometry wobble
  // shouldn't kick the base into "abnormal" placement.
  const epsilon = 0.5;
  if (y[t2] >= topY - epsilon) return;

  // Place R3 directly above terminalNodes[2] (so the bond from
  // terminalNodes[2] to base's eventual terminalNodes[0] is purely
  // vertical), at least one bond length above the topmost atom.
  const verticalClearance = 1.5;
  x[r3] = x[t2];
  y[r3] = topY + verticalClearance;
}

function adjustBaseMonomerGraph(monomer: MolGraph, pointerToBranchAngle: NumberWrapper): void {
  const x = monomer.atoms.x;
  const y = monomer.atoms.y;

  const centeredNode = monomer.meta.terminalNodes[0] - 1; // node indexing in molfiles starts from 1
  const rotatedNode = monomer.meta.rNodes[0] - 1;

  // center graph at centeredNode
  shiftCoordinates(monomer, -x[centeredNode], -y[centeredNode]);

  // rotate so that the branch bond is aligned with the sugar's branch
  const baseBranchToOYAngle = findAngleWithOY(x[rotatedNode], y[rotatedNode]);
  const sugarBranchToOYAngle = pointerToBranchAngle.value;
  if (sugarBranchToOYAngle === null)
    throw new Error('The value of sugarBranchToOYAngle is null');
  // Note: 0 is a valid value (R3 directly along +Y after the sugar fix);
  // the explicit null check above lets us treat it correctly.
  rotateCenteredGraph(monomer.atoms,
    Math.PI - baseBranchToOYAngle + sugarBranchToOYAngle);

  // scale graph in case its size does not fit the scale of phosphate and sugar
  // todo: consider extending to other monomer types
  const p1 = {
    x: x[monomer.meta.rNodes[0] - 1],
    y: y[monomer.meta.rNodes[0] - 1],
  };
  const p2 = {
    x: x[monomer.meta.terminalNodes[0] - 1],
    y: y[monomer.meta.terminalNodes[0] - 1],
  };
  const bondLength = getEuclideanDistance(p1, p2);
  if (bondLength != 1) {
    for (let i = 0; i < x.length; ++i) {
      x[i] = keepPrecision(x[i] / bondLength);
      y[i] = keepPrecision(y[i] / bondLength);
    }
  }
}

function getEuclideanDistance(p1: Point, p2: Point): number {
  return keepPrecision(Math.sqrt(
    (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2
  ));
}

/** Flip carboxyl group with the radical in a peptide monomer in case the
 * carboxyl group is in the lower half-plane
 * @param {MolGraph}monomer - peptide monomer
 * @param {number}doubleBondedOxygen - index of the double-bonded oxygen atom*/
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


/** Finds angle between OY and the ray joining origin with (x, y)
 * @param {number}x
 * @param {number}y
 * @return {number} angle in radians*/
function findAngleWithOY(x: number, y: number): number {
  let angle;
  if (x === 0)
    angle = y > 0 ? 0 : Math.PI;
  else if (y === 0)
    angle = x > 0 ? -Math.PI / 2 : Math.PI / 2;
  else {
    const tan = y / x;
    const atan = Math.atan(tan);
    angle = (x < 0) ? Math.PI / 2 + atan : -Math.PI / 2 + atan;
  }
  return angle;
}

/** Finds angle between OX and the ray joining origin with (x, y)
 * @param {number}x
 * @param {number}y
 * @return {number} angle in radians
 */
function findAngleWithOX(x: number, y: number): number {
  return findAngleWithOY(x, y) + Math.PI / 2;
}

/**  Rotate the graph around the origin by 'angle'
 * @param {Atoms}atoms - atoms of the graph
 * @param {number}angle - angle in radians*/
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

/** Flip monomer graph around OX axis preserving stereometry
 * @param {MolGraph}monomer - monomer graph*/
function flipMonomerAroundOX(monomer: MolGraph): void {
  flipMolGraph(monomer, true);
}

/** Flip monomer graph around OY axis preserving stereometry
 * @param {MolGraph}monomer - monomer graph*/
function flipMonomerAroundOY(monomer: MolGraph): void {
  flipMolGraph(monomer, false);
}

/** Flip graph around a specified axis: 'true' corresponds to OX, 'false' to OY
 * @param {MolGraph}molGraph - graph to flip
 * @param {boolean}axis - axis to flip around*/
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

function getAngleBetweenSugarBranchAndOY(molGraph: MolGraph): number {
  const x = molGraph.atoms.x;
  const y = molGraph.atoms.y;
  const rNode = molGraph.meta.rNodes[2] - 1;
  const terminalNode = molGraph.meta.terminalNodes[2] - 1;

  const xShift = x[rNode] - x[terminalNode];
  const yShift = y[rNode] - y[terminalNode];

  return Math.atan2(xShift, yShift); // notice the fliped order of arguments
}

/** Flips double-bonded 'O' in carbonyl group with 'OH' in order for the monomers
 * to have standard representation simplifying their concatenation. The
 * monomer must already be adjusted with adjustPeptideMonomerGraph in order for this function to be implemented
 * @param {MolGraph}monomer - peptide monomer
 * @param {number}doubleBondedOxygen - index of the double-bonded oxygen atom*/
function flipHydroxilGroup(monomer: MolGraph, doubleBondedOxygen: number): void {
  const x = monomer.atoms.x;
  // -1 below because indexing of nodes in molfiles starts from 1, unlike arrays
  if (x[monomer.meta.rNodes[1] - 1] > x[doubleBondedOxygen - 1])
    swapNodes(monomer, doubleBondedOxygen, monomer.meta.rNodes[1]);
}

/** Determine the number of node (starting from 1) corresponding to the
 * double-bonded oxygen of the carbonyl group
 * @param {MolGraph}monomer - peptide monomer
 * @return {number} index of the double-bonded oxygen atom*/
function findDoubleBondedCarbonylOxygen(monomer: MolGraph): number | null {
  const bondsMap = constructBondsMap(monomer);
  let doubleBondedOxygen = 0;
  const sentinel = monomer.atoms.atomTypes.length;
  let i = 0;
  if (monomer.meta.terminalNodes.length < 2)
    return null;
  // iterate over the nodes bonded to the carbon and find the double one
  while (doubleBondedOxygen === 0) {
    const node = bondsMap.get(monomer.meta.terminalNodes[1])![i];
    if (monomer.atoms.atomTypes[node - 1] === C.OXYGEN && node !== monomer.meta.rNodes[1])
      doubleBondedOxygen = node;
    i++;
    if (i > sentinel) {
      // throw new Error(`Search for double-bonded Oxygen in ${monomer} has exceeded the limit of ${sentinel}`);
      return null;
    }
  }
  return doubleBondedOxygen;
}

/** Swap the Cartesian coordinates of the two specified nodes in MolGraph
 * @param {MolGraph}monomer - monomer graph
 * @param {number}nodeOne - index of the first node
 * @param {number}nodeTwo - index of the second node*/
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

/** Maps a node to the list of nodes bound to it
 * @param {MolGraph}monomer - monomer graph
 * @return {Map<number, Array<number>>} map of nodes to the list of nodes bound to them*/
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

/** Shift molGraph in the XOY plane
 * @param {MolGraph}molGraph - graph to shift
 * @param {number}xShift - shift along X axis
 * @param {number}yShift - shift along Y axis*/
function shiftCoordinates(molGraph: MolGraph, xShift: number, yShift?: number): void {
  const x = molGraph.atoms.x;
  const y = molGraph.atoms.y;
  for (let i = 0; i < x.length; ++i) {
    x[i] = keepPrecision(x[i] + xShift);
    if (typeof yShift !== 'undefined')
      y[i] = keepPrecision(y[i] + yShift);
  }
}


export function convertMolGraphToMolfileV3K(molGraph: MolGraph): string {
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
  const molfileCountsLine = C.V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + C.V3K_COUNTS_LINE_ENDING;

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

    const atomLine = C.V3K_BEGIN_DATA_LINE + atomIdx + ' ' + atomType[i] + ' ' +
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
    const bondLine = C.V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
      firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
    molfileBondBlock += bondLine;
  }

  const molfileParts = [
    C.V3K_HEADER_FIRST_LINE,
    C.V3K_HEADER_SECOND_LINE,
    C.V3K_BEGIN_CTAB_BLOCK,
    molfileCountsLine,
    C.V3K_BEGIN_ATOM_BLOCK,
    molfileAtomBlock,
    C.V3K_END_ATOM_BLOCK,
    C.V3K_BEGIN_BOND_BLOCK,
    molfileBondBlock,
    C.V3K_END_BOND_BLOCK,
    C.V3K_END_CTAB_BLOCK,
    C.V3K_END,
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

    const molfileV3K = monomerLibObject[HELM_FIELDS.MOLFILE].includes('V3000') ? monomerLibObject[HELM_FIELDS.MOLFILE] :
      convertMolfileToV3K(removeRGroupLines(monomerLibObject[HELM_FIELDS.MOLFILE]), moduleRdkit);
    const counts = parseAtomAndBondCounts(molfileV3K);

    const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
    const bonds = parseBondBlock(molfileV3K, counts.bondCount);
    const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);

    const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

    removeHydrogen(monomerGraph);
    replaceWrongfulRGroups(monomerGraph);
    const molfile = convertMolGraphToMolfileV3K(monomerGraph);
    symbolToCappedMolfileMap.set(monomerSymbol, molfile);
  }
  return symbolToCappedMolfileMap;
}

/** Get the V3K molfile corresponding to the capped Monomer (default cap groups)
 * @param {Monomer} monomer
 * @return {string} V3K molfile*/
export function capPeptideMonomer(monomer: Monomer): string {
  const funcList: DG.Func[] = DG.Func.find({package: 'Chem', name: 'getRdKitModule'});
  const moduleRdkit = funcList[0].apply();

  const capGroups = parseCapGroups(monomer[HELM_FIELDS.RGROUPS]);
  const capGroupIdxMap = parseCapGroupIdxMap(monomer[HELM_FIELDS.MOLFILE]);
  const molfileV3K = monomer[HELM_FIELDS.MOLFILE].includes('V3000') ? monomer[HELM_FIELDS.MOLFILE] :
    convertMolfileToV3K(removeRGroupLines(monomer[HELM_FIELDS.MOLFILE]), moduleRdkit);
  const counts = parseAtomAndBondCounts(molfileV3K);

  const atoms = parseAtomBlock(molfileV3K, counts.atomCount);
  const bonds = parseBondBlock(molfileV3K, counts.bondCount);
  const meta = getMonomerMetadata(atoms, bonds, capGroups, capGroupIdxMap);

  const monomerGraph: MolGraph = {atoms: atoms, bonds: bonds, meta: meta};

  adjustPeptideMonomerGraph(monomerGraph);

  const molfile = convertMolGraphToMolfileV3K(monomerGraph);
  return molfile;
}
