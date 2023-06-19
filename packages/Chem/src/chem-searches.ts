import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitService} from './utils/chem-common-rdkit';
import {
  chemBeginCriticalSection,
  chemEndCriticalSection,
  defaultMorganFpLength,
  defaultMorganFpRadius,
  Fingerprint,
  rdKitFingerprintToBitArray,
} from './utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {assure} from '@datagrok-libraries/utils/src/test';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import { getMolSafe, getQueryMolSafe } from './utils/mol-creation_rdkit';
import { _package } from './package';

const enum FING_COL_TAGS {
  invalidatedForVersion = '.invalideted.for.version',
  molsCreatedForVersion = '.mols.created.for.version',
  substrLibCreatedForVersion = 'substr.lib.created.for.version',
}

let lastColumnInvalidated: string = '';

function _chemFindSimilar(molStringsColumn: DG.Column, fingerprints: (BitArray | null)[],
  queryMolString: string, settings: { [name: string]: any }): DG.DataFrame {
  const len = molStringsColumn.length;
  const distances = _chemGetSimilarities(queryMolString, fingerprints);
  const limit = Math.min((settings.hasOwnProperty('limit') ? settings.limit : len), len);
  const minScore = settings.hasOwnProperty('minScore') ? settings.minScore : 0.0;
  const sortedIndices = Array.from(Array(len).keys()).sort((i1, i2) => {
    const a1 = distances[i1];
    const a2 = distances[i2];
    if (a2 < a1) return -1;
    if (a2 > a1) return +1;
    return 0;
  });
  const sortedMolStrings = DG.Column.fromType(DG.TYPE.STRING, 'molecule', limit);
  const sortedMolInd = DG.Column.fromType(DG.TYPE.INT, 'index', limit);
  sortedMolStrings.semType = DG.SEMTYPE.MOLECULE;
  const sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', limit);
  for (let n = 0; n < limit; n++) {
    const idx = sortedIndices[n];
    const score = distances[idx];
    if (score < minScore) {
      sortedMolStrings.dataFrame.rows.removeAt(n, limit - n);
      sortedScores.dataFrame.rows.removeAt(n, limit - n);
      break;
    }
    sortedMolStrings.set(n, molStringsColumn.get(idx));
    sortedScores.set(n, score);
    sortedMolInd.set(n, idx);
  }
  return DG.DataFrame.fromColumns([sortedMolStrings, sortedScores, sortedMolInd]);
}

// Only this function receives {sorted} in settings
export function _chemGetSimilarities(queryMolString: string, fingerprints: (BitArray | null)[]): number[] {
  const distances = new Array(fingerprints.length).fill(0.0);
  try {
    const sample = chemGetFingerprint(queryMolString, Fingerprint.Morgan);
    for (let i = 0; i < fingerprints.length; ++i)
      distances[i] = tanimotoSimilarity(fingerprints[i] ?? BitArray.fromBytes(new Uint8Array()), sample);
  } catch {
    return distances;
  }
  return distances;
}

function _chemGetDiversities(limit: number, molStringsColumn: DG.Column, fingerprints: (BitArray | null)[]): string[] {
  limit = Math.min(limit, fingerprints.length);
  const indexes = ArrayUtils.indexesOf(fingerprints, (f) => f != null);
  const diverseIndexes = getDiverseSubset(indexes.length, limit,
    (i1: number, i2: number) => 1 - tanimotoSimilarity(fingerprints[indexes[i1]]!, fingerprints[indexes[i2]]!));

  const diversities = new Array(limit).fill('');

  for (let i = 0; i < limit; i++)
    diversities[i] = molStringsColumn.get(indexes[diverseIndexes[i]]);

  return diversities;
}

function colInvalidated(col: DG.Column, createMols: boolean, useSubstructLib?: boolean): Boolean {
  return (lastColumnInvalidated == col.name &&
    col.getTag(FING_COL_TAGS.invalidatedForVersion) == String(col.version) &&
    (!createMols || useSubstructLib ?
      col.getTag(FING_COL_TAGS.substrLibCreatedForVersion) == String(col.version) :
      col.getTag(FING_COL_TAGS.molsCreatedForVersion) == String(col.version)));
}

async function _invalidate(molCol: DG.Column, createMols: boolean, useSubstructLib?: boolean) {
  if (!colInvalidated(molCol, createMols, useSubstructLib)) {
    if (createMols) {
      await (await getRdKitService()).initMoleculesStructures(molCol.toList(), useSubstructLib);
      useSubstructLib ? molCol.setTag(FING_COL_TAGS.substrLibCreatedForVersion, String(molCol.version + 2)) :
        molCol.setTag(FING_COL_TAGS.molsCreatedForVersion, String(molCol.version + 2));
    }
    lastColumnInvalidated = molCol.name;
    molCol.setTag(FING_COL_TAGS.invalidatedForVersion, String(molCol.version + 1));
  }
}

function checkForFingerprintsColumn(col: DG.Column, fingerprintsType: Fingerprint): Uint8Array[] | null {
  const colNameTag = '.' + fingerprintsType + '.Column';
  const colVerTag = '.' + fingerprintsType + '.Version';

  if (col.getTag(colNameTag) &&
      col.getTag(colVerTag) == String(col.version) &&
      col.dataFrame) {
    const fingCol = col.dataFrame.columns.byName(col.getTag(colNameTag));
    if (fingCol)
      return fingCol.toList();
  }
  return null;
}

function saveFingerprintsToCol(col: DG.Column, fgs: (Uint8Array | null)[],
  fingerprintsType: Fingerprint, createdMols: boolean): void {
  if (!col.dataFrame)
    throw new Error('Column has no parent dataframe');

  const colNameTag = '.' + fingerprintsType + '.Column';
  const colVerTag = '.' + fingerprintsType + '.Version';

  const fingerprintColumnName = col.getTag(colNameTag) ??
    '~' + col.name + fingerprintsType + 'Fingerprints';
  const df = col.dataFrame;
  const newCol: DG.Column<Uint8Array> = df.columns.getOrCreate(
    fingerprintColumnName, DG.COLUMN_TYPE.BYTE_ARRAY, fgs.length);

  newCol.init((i) => fgs[i]);

  col.setTag(colNameTag, fingerprintColumnName);
  if (createdMols)
    col.setTag(FING_COL_TAGS.molsCreatedForVersion, String(col.version + 3));
  col.setTag(colVerTag, String(col.version + 2));
  col.setTag(FING_COL_TAGS.invalidatedForVersion, String(col.version + 1));
}

async function getUint8ArrayFingerprints(
  molCol: DG.Column, fingerprintsType: Fingerprint = Fingerprint.Morgan,
  useSection = true, createMols = true): Promise<(Uint8Array | null)[]> {
  if (useSection)
    await chemBeginCriticalSection();
  try {
    const fgsCheck = checkForFingerprintsColumn(molCol, fingerprintsType);
    if (fgsCheck)
      return fgsCheck;
    else {
      await _invalidate(molCol, createMols);
      const molecules = createMols ? undefined : molCol.toList();
      const fingerprints = await (await getRdKitService()).getFingerprints(fingerprintsType, molecules);
      saveFingerprintsToCol(molCol, fingerprints, fingerprintsType, createMols);
      chemEndCriticalSection();
      return fingerprints;
    }
  } finally {
    if (useSection)
      chemEndCriticalSection();
  }
}

function substructureSearchPatternsMatch(molString: string, querySmarts: string, fgs: (Uint8Array | null)[]): BitArray {
  const patternFpUint8Length = 256;
  const result = new BitArray(fgs.length, false);
  const queryMol = getQueryMolSafe(molString, querySmarts, getRdKitModule());

  if (queryMol) {
    try {
      const fpRdKit = queryMol.get_pattern_fp_as_uint8array();
      checkEl:
      for (let i = 0; i < fgs.length; ++i) {
        if (fgs[i]) {
          for (let j = 0; j < patternFpUint8Length; ++j) {
            if ((fgs[i]![j] & fpRdKit[j]) != fpRdKit[j])
              continue checkEl;
          }
          result.setBit(i, true, false);
        }
      }
      return result;
    } catch (e: any) { 
      _package.logger.error(`substructureSearchPatternsMatch failed with error: ${e.toString()}`);
    }
  }
  return result;
}

export async function chemGetFingerprints(...args: [DG.Column, Fingerprint?, boolean?, boolean?]):
  Promise<(BitArray | null)[]> {
  return (await getUint8ArrayFingerprints(...args)).map((el) => el ? rdKitFingerprintToBitArray(el) : null);
}

export async function chemGetSimilarities(molStringsColumn: DG.Column, queryMolString = '')
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, true, false)!;

  return queryMolString.length != 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
}

export async function chemGetDiversities(molStringsColumn: DG.Column, limit: number)
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(limit, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, true, false)!;

  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molecule',
    _chemGetDiversities(limit, molStringsColumn, fingerprints));
}

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, true, false)!;
  return queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molBlockFailover: string, usePatternFingerprints = false,
  useSubstructLib?: boolean)
  : Promise<DG.BitSet> {
  await chemBeginCriticalSection();
  try {
    let matches = DG.BitSet.create(molStringsColumn.length);
    if (molString.length != 0) {
      if (usePatternFingerprints) {
        const fgs: (Uint8Array | null)[] = await getUint8ArrayFingerprints(molStringsColumn, Fingerprint.Pattern, false, false);
        const filteredMolsIdxs: BitArray = substructureSearchPatternsMatch(molString, molBlockFailover, fgs);
        const filteredMolecules: string[] = getMoleculesFilteredByPatternFp(molStringsColumn, filteredMolsIdxs);
        const matchesBitArray: DG.BitSet = await (await getRdKitService()).searchSubstructure(molString, molBlockFailover, filteredMolecules, false);
        matches = restoreMatchesByFilteredIdxs(filteredMolsIdxs, matchesBitArray);
      } else {
        await _invalidate(molStringsColumn, true, useSubstructLib);
        matches = await (await getRdKitService()).searchSubstructure(molString, molBlockFailover, undefined, useSubstructLib);
      }
    }
    return matches;
  } finally {
    chemEndCriticalSection();
  }
}

function getMoleculesFilteredByPatternFp(molStringsColumn: DG.Column, filteredMolsIdxs: BitArray): string[] {
  const molStrings = molStringsColumn.toList();
  const filteredMolecules = Array<string>(filteredMolsIdxs.trueCount());
  let counter = 0;
  for (let i = filteredMolsIdxs.getBit(0) ? 0 : filteredMolsIdxs.findNext(0); i !== -1; i = filteredMolsIdxs.findNext(i)) {
      filteredMolecules[counter] = molStrings[i];
      counter++;
  }
  return filteredMolecules;
}

function restoreMatchesByFilteredIdxs(filteredMolsIdxs: BitArray, matchesBitSet: DG.BitSet): DG.BitSet {
  const matches = DG.BitSet.create(filteredMolsIdxs.length);
  let matchesCounter = 0;
  for (let i = filteredMolsIdxs.getBit(0) ? 0 : filteredMolsIdxs.findNext(0); i !== -1; i = filteredMolsIdxs.findNext(i)) {
    if (matchesBitSet.get(matchesCounter))
        matches.set(i, true);
    matchesCounter++;
  }
  return matches;
}

export function chemGetFingerprint(molString: string, fingerprint: Fingerprint): BitArray {
  let mol = null;
  try {
    mol = getMolSafe(molString, {}, getRdKitModule()).mol;
    if (mol) {
      let fp;
      if (fingerprint == Fingerprint.Morgan) {
        fp = mol.get_morgan_fp_as_uint8array(JSON.stringify({
          radius: defaultMorganFpRadius,
          nBits: defaultMorganFpLength,
        }));
      } else if (fingerprint == Fingerprint.Pattern)
        fp = mol.get_pattern_fp_as_uint8array();
      else
        throw new Error(`${fingerprint} does not match any fingerprint`);
  
      return rdKitFingerprintToBitArray(fp) as BitArray;
    } else
      throw new Error(`Chem | Possibly a malformed molString: ${molString}`);
  } catch {
    throw new Error(`Chem | Possibly a malformed molString: ${molString}`);
  } finally {
    mol?.delete();
  }
}

export async function geMolNotationConversions(molCol: DG.Column, targetNotation: string): Promise<string[]> {
  await chemBeginCriticalSection();
  try {
    await _invalidate(molCol, true);
    const conversions = await (await getRdKitService()).convertMolNotation(targetNotation);
    //saveFingerprintsToCol(molCol, fingerprints, fingerprintsType);
    chemEndCriticalSection();
    return conversions;
  } finally {
    chemEndCriticalSection();
  }
}
