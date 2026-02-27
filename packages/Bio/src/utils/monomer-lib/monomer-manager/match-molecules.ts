/* eslint-disable max-lines-per-function */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';

import {STANDRARD_R_GROUPS} from './const';
import {standardiseMonomers, capSmiles, getCorrectedSmiles} from './monomer-manager';

/** Represents a single monomer match result — one monomer that matched a molecule */
interface MonomerMatch {
  /** monomer symbol from the library */
  symbol: string;
  /** canonical (possibly capped) SMILES used for the match */
  smiles: string;
  /** original SMILES from the monomer definition */
  original: string;
  /** library source name */
  source: string;
}

/** Maps keyed by canonical SMILES, where each key can map to multiple monomers */
type MonomerSmilesMap = {[smiles: string]: MonomerMatch[]};

const MATCH_SEPARATOR = ' | ';

/**
 * Builds lookup maps from standardized monomers:
 * - uncappedMap: maps raw canonical monomer SMILES -> MonomerMatch[]
 * - cappedMap: maps capped (R-groups replaced with cap groups) canonical SMILES -> MonomerMatch[]
 * Both maps store arrays so that duplicate monomers (same structure, different symbols/libs) are preserved.
 */
async function buildMonomerSmilesMaps(
  fixedMonomers: Monomer[], originalMonomers: Monomer[], converterFunc: DG.Func,
): Promise<{cappedMap: MonomerSmilesMap; uncappedMap: MonomerSmilesMap}> {
  // build uncapped map from raw monomer SMILES
  const uncappedMap: MonomerSmilesMap = {};
  for (const m of fixedMonomers) {
    if (!m.smiles) continue;
    const match: MonomerMatch = {symbol: m.symbol, smiles: m.smiles, original: m.smiles, source: m.lib?.source ?? ''};
    if (!uncappedMap[m.smiles]) uncappedMap[m.smiles] = [];
    uncappedMap[m.smiles].push(match);
  }

  // build capped monomer entries: replace R-groups with cap group atoms
  const cappedEntries = fixedMonomers
    .map((m, i) => ({
      symbol: m.symbol,
      smiles: capSmiles(m.smiles ?? '', m.rgroups ?? []),
      original: m.smiles,
      source: originalMonomers[i]?.lib?.source ?? '',
    }))
    .filter((e) => !!e.smiles && !e.smiles.includes('[*:'));

  // canonicalize all capped SMILES in bulk
  const cappedSmilesCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CappedSmiles', cappedEntries.map((e) => e.smiles));
  cappedSmilesCol.semType = DG.SEMTYPE.MOLECULE;
  const canonicalCappedCol: DG.Column = await converterFunc.apply({molecule: cappedSmilesCol, targetNotation: DG.chem.Notation.Smiles});
  if (!canonicalCappedCol || canonicalCappedCol.length !== cappedSmilesCol.length)
    throw new Error('Error canonicalizing capped monomer SMILES');

  // build capped map with canonicalized SMILES as keys
  const cappedMap: MonomerSmilesMap = {};
  const canonicalCappedList = canonicalCappedCol.toList();
  for (let i = 0; i < canonicalCappedList.length; i++) {
    const smi = canonicalCappedList[i];
    if (!smi) continue;
    cappedEntries[i].smiles = smi;
    const match: MonomerMatch = cappedEntries[i];
    if (!cappedMap[smi]) cappedMap[smi] = [];
    cappedMap[smi].push(match);
  }

  return {cappedMap, uncappedMap};
}

/**
 * Corrects and canonicalizes the input molecule column.
 * Handles both SMILES and molblock inputs.
 * Returns the list of canonical SMILES strings (null for invalid molecules).
 */
async function canonicalizeMolecules(
  molDf: DG.DataFrame, molColName: string, converterFunc: DG.Func,
): Promise<(string | null)[]> {
  const moleculesOriginalCol = molDf.col(molColName)!;
  const correctedList = moleculesOriginalCol.toList().map((s) => {
    if (!s) return s;
    try {
      const isMolBlock = s.includes('\n');
      return getCorrectedSmiles([], isMolBlock ? undefined : s, isMolBlock ? s : undefined);
    } catch (_e) {
      return s;
    }
  });

  const correctedCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MoleculesCorrected', correctedList);
  correctedCol.semType = DG.SEMTYPE.MOLECULE;
  // dummy df needed for semtype detection by converterFunc
  const _ddf = DG.DataFrame.fromColumns([correctedCol]);

  const canonicalCol: DG.Column = await converterFunc.apply({molecule: correctedCol, targetNotation: DG.chem.Notation.Smiles});
  if (!canonicalCol || canonicalCol.length !== correctedCol.length)
    throw new Error('Error canonicalizing molecules');

  return canonicalCol.toList();
}

/**
 * Attempts to match a single canonical molecule SMILES against the lookup maps.
 * Tries in order: capped map -> uncapped map -> cap the molecule with standard R-groups and retry.
 * Returns all matching monomers (can be multiple from different libraries).
 */
function matchBySmiles(
  canonicalMol: string, cappedMap: MonomerSmilesMap, uncappedMap: MonomerSmilesMap,
): MonomerMatch[] {
  // try direct lookup in capped and uncapped maps
  let matches = cappedMap[canonicalMol] ?? uncappedMap[canonicalMol];
  if (matches && matches.length > 0) return matches;

  // fallback: cap the molecule with standard R-groups and try again
  const cappedMol = capSmiles(canonicalMol, STANDRARD_R_GROUPS);
  if (cappedMol !== canonicalMol) {
    const correctedMol = grok.chem.convert(cappedMol, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
    matches = cappedMap[correctedMol] ?? uncappedMap[correctedMol];
    if (matches && matches.length > 0) return matches;
  }

  return [];
}

/**
 * Builds a Morgan fingerprint lookup map for all capped monomer SMILES.
 * The map keys are fingerprint binary strings (via DG.BitSet.toBinaryString()),
 * which allows fast exact matching that is tolerant of explicit hydrogen
 * and minor stereochemistry differences.
 */
async function buildMonomerFingerprintMap(
  cappedMap: MonomerSmilesMap,
): Promise<{fpMap: {[fpString: string]: MonomerMatch[]}; cappedSmilesList: string[]}> {
  const cappedSmilesList = Object.keys(cappedMap);
  if (cappedSmilesList.length === 0)
    return {fpMap: {}, cappedSmilesList: []};

  const monomerCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MonomerSmiles', cappedSmilesList);
  monomerCol.semType = DG.SEMTYPE.MOLECULE;

  const fpCol: DG.Column = await grok.functions.call('Chem:getMorganFingerprints', {molColumn: monomerCol});

  const fpMap: {[fpString: string]: MonomerMatch[]} = {};
  for (let i = 0; i < fpCol.length; i++) {
    const fp: DG.BitSet | null = fpCol.get(i);
    if (!fp) continue;
    const fpStr = fp.toBinaryString();
    // merge monomer matches from the SMILES map into the fingerprint map
    const smilesMatches = cappedMap[cappedSmilesList[i]] ?? [];
    if (!fpMap[fpStr]) fpMap[fpStr] = [];
    fpMap[fpStr].push(...smilesMatches);
  }

  return {fpMap, cappedSmilesList};
}

/**
 * For molecules that were not matched by exact SMILES, attempts matching via
 * Morgan fingerprints. Computes fingerprints for unmatched molecules and looks
 * them up in the monomer fingerprint map. Also tries capping with standard R-groups.
 */
async function matchByFingerprint(
  unmatchedIndices: number[],
  canonicalizedMolecules: (string | null)[],
  monomerFpMap: {[fpString: string]: MonomerMatch[]},
): Promise<Map<number, MonomerMatch[]>> {
  const results = new Map<number, MonomerMatch[]>();
  if (unmatchedIndices.length === 0 || Object.keys(monomerFpMap).length === 0)
    return results;

  // collect SMILES for unmatched molecules (uncapped first)
  const uncappedSmiles: string[] = unmatchedIndices.map((idx) => canonicalizedMolecules[idx] ?? '');

  // also prepare capped versions
  const cappedSmiles: string[] = uncappedSmiles.map((s) => {
    if (!s) return '';
    const capped = capSmiles(s, STANDRARD_R_GROUPS);
    return capped !== s ? grok.chem.convert(capped, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) : s;
  });

  // compute fingerprints for both uncapped and capped molecules
  const uncappedCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'UnmatchedMols', uncappedSmiles);
  uncappedCol.semType = DG.SEMTYPE.MOLECULE;
  const cappedCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'UnmatchedMolsCapped', cappedSmiles);
  cappedCol.semType = DG.SEMTYPE.MOLECULE;

  const [uncappedFpCol, cappedFpCol]: [DG.Column, DG.Column] = await Promise.all([
    grok.functions.call('Chem:getMorganFingerprints', {molColumn: uncappedCol}),
    grok.functions.call('Chem:getMorganFingerprints', {molColumn: cappedCol}),
  ]);

  for (let i = 0; i < unmatchedIndices.length; i++) {
    const molIdx = unmatchedIndices[i];

    // try uncapped fingerprint first
    const uncappedFp: DG.BitSet | null = uncappedFpCol.get(i);
    if (uncappedFp) {
      const fpStr = uncappedFp.toBinaryString();
      const matches = monomerFpMap[fpStr];
      if (matches && matches.length > 0) {
        results.set(molIdx, matches);
        continue;
      }
    }

    // fallback: try capped fingerprint
    const cappedFp: DG.BitSet | null = cappedFpCol.get(i);
    if (cappedFp) {
      const fpStr = cappedFp.toBinaryString();
      const matches = monomerFpMap[fpStr];
      if (matches && matches.length > 0)
        results.set(molIdx, matches);
    }
  }

  return results;
}

/** Deduplicates matches by symbol, keeping one entry per unique monomer symbol */
function deduplicateMatches(matches: MonomerMatch[]): MonomerMatch[] {
  const seen = new Set<string>();
  return matches.filter((m) => {
    if (seen.has(m.symbol)) return false;
    seen.add(m.symbol);
    return true;
  });
}

/** Collects all source library names for matched monomers, including known duplicates */
function collectSources(
  matches: MonomerMatch[], duplicates: {[symbol: string]: Monomer[]},
): string {
  const sources = new Set<string>();
  for (const m of matches) {
    // check if monomerLib knows about duplicates for this symbol across libraries
    const dups = duplicates[m.symbol];
    if (dups && dups.length > 0) {
      for (const dup of dups) {
        const s = dup?.lib?.source;
        if (s) sources.add(s);
      }
    } else if (m.source)
      sources.add(m.source);
  }
  return Array.from(sources).join(', ');
}

/**
 * Matches molecules in a dataframe with monomers from a monomer library.
 *
 * Matching pipeline:
 *   1. Standardize monomers and build SMILES lookup maps (capped & uncapped)
 *   2. Canonicalize input molecules
 *   3. Phase 1: exact canonical SMILES matching (capped, uncapped, and fallback-capped molecule)
 *   4. Phase 2: Morgan fingerprint fallback for molecules that didn't match by SMILES
 *   5. Populate result columns (supports multiple matches per molecule via pipe-delimited values)
 *
 * @returns cloned input DataFrame with added match columns
 */
export async function matchMoleculesWithMonomers(
  molDf: DG.DataFrame, molColName: string, monomerLib: IMonomerLib, polymerType: PolymerType = 'PEPTIDE',
): Promise<DG.DataFrame> {
  const duplicates = monomerLib.duplicateMonomers?.[polymerType] ?? {};
  const converterFunc = DG.Func.find({package: 'Chem', name: 'convertMoleculeNotation'})[0];
  if (!converterFunc)
    throw new Error('Function convertMoleculeNotation not found, please install Chem package');

  // === Step 1: Standardize monomers and build SMILES lookup maps ===
  const monomers = monomerLib.getMonomerSymbolsByType(polymerType)
    .map((s) => monomerLib.getMonomer(polymerType, s)!)
    .filter((m) => m && (m.smiles || m.molfile));

  const fixedMonomers = await standardiseMonomers(monomers);
  // preserve library reference from original monomers (lost during standardization)
  fixedMonomers.forEach((m, i) => { m.lib = monomers[i].lib; });

  const {cappedMap, uncappedMap} = await buildMonomerSmilesMaps(fixedMonomers, monomers, converterFunc);

  // === Step 2: Canonicalize input molecules ===
  const canonicalizedMolecules = await canonicalizeMolecules(molDf, molColName, converterFunc);

  // === Step 3: Phase 1 — Exact canonical SMILES matching ===
  // matchResults[i] holds all MonomerMatch entries for molecule i (empty array if unmatched)
  const matchResults: MonomerMatch[][] = new Array(canonicalizedMolecules.length).fill(null).map(() => []);
  const unmatchedIndices: number[] = [];

  for (let i = 0; i < canonicalizedMolecules.length; i++) {
    const mol = canonicalizedMolecules[i];
    if (!mol) continue;
    const smilesMatches = matchBySmiles(mol, cappedMap, uncappedMap);
    if (smilesMatches.length > 0)
      matchResults[i] = smilesMatches;
    else
      unmatchedIndices.push(i);
  }

  // === Step 4: Phase 2 — Morgan fingerprint fallback for unmatched molecules ===
  if (unmatchedIndices.length > 0) {
    try {
      const {fpMap} = await buildMonomerFingerprintMap(cappedMap);
      const fpMatches = await matchByFingerprint(unmatchedIndices, canonicalizedMolecules, fpMap);
      for (const [idx, matches] of fpMatches)
        matchResults[idx] = matches;
    } catch (e) {
      console.warn('Fingerprint fallback matching failed, continuing with SMILES matches only:', e);
    }
  }

  // === Step 5: Populate result columns ===
  const resultDf = molDf.clone();
  const symbolCol = resultDf.columns.addNewString(resultDf.columns.getUnusedName('Matched monomer symbol'));
  const smilesCol = resultDf.columns.addNewString(resultDf.columns.getUnusedName('Matched monomer smiles'));
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const sourceCol = resultDf.columns.addNewString(resultDf.columns.getUnusedName('Matched monomer source'));
  const matchCountCol = resultDf.columns.addNewInt(resultDf.columns.getUnusedName('Match count'));
  const matchMethodCol = resultDf.columns.addNewString(resultDf.columns.getUnusedName('Match method'));
  resultDf.columns.setOrder([molColName, symbolCol.name, smilesCol.name, sourceCol.name, matchCountCol.name, matchMethodCol.name]);

  for (let i = 0; i < matchResults.length; i++) {
    const matches = matchResults[i];
    if (matches.length === 0) continue;

    // deduplicate matches by symbol (same monomer can appear from multiple lookup paths)
    const uniqueMatches = deduplicateMatches(matches);

    // collect all sources, including duplicates from the monomer library
    const allSources = collectSources(uniqueMatches, duplicates);

    symbolCol.set(i, uniqueMatches.map((m) => m.symbol).join(MATCH_SEPARATOR), false);
    smilesCol.set(i, uniqueMatches[0].original ?? uniqueMatches[0].smiles, false);
    sourceCol.set(i, allSources, false);
    matchCountCol.set(i, uniqueMatches.length, false);
    // fingerprint matches are those from phase 2 (indices that were in unmatchedIndices)
    const method = unmatchedIndices.includes(i) ? 'fingerprint' : 'exact';
    matchMethodCol.set(i, method, false);
  }

  return resultDf;
}
