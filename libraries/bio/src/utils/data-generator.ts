import * as DG from 'datagrok-api/dg';
import {NOTATION} from './macromolecule';
//@ts-ignore: no types
import * as jStat from 'jstat';

export type DistributionType = 'lognormal' | 'normal';
export type DistributionParams = {type: DistributionType, mean: number, std: number};

export const colNames = {sequences: 'sequences', clusters: 'clusters', activity: 'activity'};

const sampleFrom: {[key in DistributionType]: (mean: number, std: number) => number} = {
  'lognormal': jStat.lognormal.sample,
  'normal': jStat.normal.sample,
};

export function generatePeptidesDataFrame(type: NOTATION, length: number, alphabet: string[],
  clusterSequenceLengths: number[] = [15], distParams: DistributionParams = {type: 'lognormal', mean: 0, std: 1},
  ): DG.DataFrame {
  const df = DG.DataFrame.create(length);
  const sequencesCol = df.columns.addNewString(colNames.sequences);
  const clustersCol = df.columns.addNewInt(colNames.clusters);
  const clustersColRawData = clustersCol.getRawData();
  const activityCol = df.columns.addNewInt(colNames.activity);
  const activityColRawData = activityCol.getRawData();

  for (let rowIdx = 0; rowIdx < length; ++rowIdx) {
    const sequenceLength = clusterSequenceLengths[Math.floor(Math.random() * clusterSequenceLengths.length)];
    clustersColRawData[rowIdx] = sequenceLength;
    activityColRawData[rowIdx] = sampleFrom[distParams.type](distParams.mean, distParams.std);
    const sequence = generateSequence(type, alphabet, sequenceLength);
    sequencesCol.set(rowIdx, sequence);
  }

  return df;
}

function generateSequence(type: NOTATION, alphabet: string[], length: number): string {
  const sequence = [];
  for (let i = 0; i < length; ++i)
    sequence.push(alphabet[Math.floor(Math.random() * alphabet.length)]);

  return buildSequence(sequence, type);
}

function buildSequence(sequence: string[], type: NOTATION = NOTATION.HELM, separator: string = '/'): string {
  switch (type) {
    case NOTATION.SEPARATOR:
      return sequence.join(separator);
    case NOTATION.FASTA:
      return sequence.join('');
    case NOTATION.HELM:
      return `PEPTIDE1{${sequence.join('.')}}$$$$`;
    default:
      throw new Error('Unknown notation type');
  }
}
