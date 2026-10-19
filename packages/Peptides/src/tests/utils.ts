import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {scaleActivity} from '../utils/misc';
import {SCALING_METHODS, TAGS} from '../utils/constants';
import {MCLSettings} from '../utils/types';

export enum TEST_COLUMN_NAMES {
  SEQUENCE = 'sequence',
  ACTIVITY = 'activity',
  CLUSTER = 'cluster',
}

export enum DEMO_COLUMN_NAMES {
  ID = 'ID',
  SEQUENCE = 'AlignedSequence',
  ACTIVITY = 'IC50',
}

/** The demo peptides: 647 aligned separator sequences of 17 positions with ID and IC50, or their first rows. */
export async function loadDemoPeptides(rowCount?: number): Promise<DG.DataFrame> {
  const csv = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
  const df = rowCount === undefined ? csv : csv.clone(DG.BitSet.create(csv.rowCount, (i) => i < rowCount));
  df.name = 'peptides';
  await grok.data.detectSemanticTypes(df);
  return df;
}

/** Starts SAR on unscaled IC50; with a similarity threshold, MCL clusters feed the Logo Summary Table. */
export async function analyzeDemoPeptides(df: DG.DataFrame, mclThreshold?: number): Promise<PeptidesModel> {
  const activityCol = df.getCol(DEMO_COLUMN_NAMES.ACTIVITY);
  const mclSettings = new MCLSettings();
  mclSettings.threshold = mclThreshold ?? mclSettings.threshold;
  const addMCL = mclThreshold !== undefined;
  const model = await startAnalysis(activityCol, df.getCol(DEMO_COLUMN_NAMES.SEQUENCE), null, df,
    scaleActivity(activityCol, SCALING_METHODS.NONE), SCALING_METHODS.NONE,
    {addMCL, useEmbeddingsClusters: addMCL, mclSettings});
  if (model === null)
    throw new Error('Model is null');
  return model;
}

/** Monomers of every separator sequence, gaps as empty strings. */
export function splitSequences(sequences: DG.Column<string>): string[][] {
  const separator = sequences.getTag(TAGS.SEPARATOR);
  return sequences.toList().map((s) => s.split(separator));
}

/** Row pairs `a:b` (a < b) of equally aligned sequences that differ at exactly one position. */
export function singleMutationPairs(split: string[][]): Set<string> {
  const pairs = new Set<string>();
  for (let a = 0; a < split.length; a++) {
    for (let b = a + 1; b < split.length; b++) {
      if (split[a].filter((monomer, position) => monomer !== split[b][position]).length === 1)
        pairs.add(`${a}:${b}`);
    }
  }
  return pairs;
}
