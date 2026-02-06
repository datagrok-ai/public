import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DescriptorStatistics, SOURCE_PATH, SYNTHETIC_DRUG_NAME} from './pmpo-defs';
import {getDescriptorStatistics, getDesiredTables} from './stat-tools';

//@ts-ignore: no types
import * as jStat from 'jstat';

/** Generates synthetic data for pMPO model training and testing
 * @param samplesCount Number of samples to generate
 * @returns DataFrame with generated data */
export async function getSynteticPmpoData(samplesCount: number): Promise<DG.DataFrame> {
  const df = await grok.dapi.files.readCsv(SOURCE_PATH);
  const generator = new PmpoDataGenerator(df, 'Drug', 'CNS', 'Smiles');

  return generator.getGenerated(samplesCount);
}

/** Class for generating synthetic data for pMPO model training and testing */
export class PmpoDataGenerator {
  private sourceDf: DG.DataFrame;
  private drugName: string;
  private desirabilityColName: string;
  private smilesColName: string;
  private desiredProbability: number;
  private descriptorStats: Map<string, DescriptorStatistics>;

  constructor(df: DG.DataFrame, drugName: string, desirabilityColName: string, smilesColName: string) {
    this.sourceDf = df;
    this.drugName = drugName;
    this.desirabilityColName = desirabilityColName;
    this.smilesColName = smilesColName;

    const descriptorNames = df.columns.toList().filter((col) => col.isNumerical).map((col) => col.name);
    const {desired, nonDesired} = getDesiredTables(df, df.col(desirabilityColName)!);

    // Compute descriptors' statistics
    this.descriptorStats = new Map<string, DescriptorStatistics>();
    descriptorNames.forEach((name) => {
      this.descriptorStats.set(name, getDescriptorStatistics(desired.col(name)!, nonDesired.col(name)!));
    });

    // Probability of desired class
    this.desiredProbability = desired.rowCount / df.rowCount;
  } // constructor

  /** Generates synthetic data for pMPO model training and testing
   * @param samplesCount Number of samples to generate
   * @returns DataFrame with generated data */
  public getGenerated(samplesCount: number): DG.DataFrame {
    if (samplesCount <= 1)
      throw new Error('Failed to generate pMPO data: sample count must be greater than 1.');

    let result: DG.DataFrame;

    /* Use rows from the source dataframe if the requested sample count
       is less than or equal to the source dataframe row count */
    if (samplesCount <= this.sourceDf.rowCount) {
      const rowMask = DG.BitSet.create(this.sourceDf.rowCount);

      for (let i = 0; i < samplesCount; ++i)
        rowMask.set(i, true);

      result = this.sourceDf.clone(rowMask);
    } else {
      const cloneDf = this.getClonedSourceDfWithFloatNumericCols();
      result = cloneDf.append(this.getSyntheticTable(samplesCount - this.sourceDf.rowCount));
    }

    // Check boolean columns and ensure non-zero stdev
    for (const col of result.columns) {
      if (col.type === DG.COLUMN_TYPE.BOOL && col.stats.stdev === 0) {
        // All values are the same, flip the first value
        let value = col.get(0);
        col.set(0, !value);

        value = col.get(1);
        col.set(1, !value);
      }
    }

    return result;
  } // getGenerated

  /** Generates a synthetic data table
   * @param samplesCount Number of samples to generate
   * @returns DataFrame with synthetic data */
  private getSyntheticTable(samplesCount: number): DG.DataFrame {
    const desirabilityRaw = new Array<boolean>(samplesCount);

    for (let i = 0; i < samplesCount; ++i)
      desirabilityRaw[i] = (Math.random() < this.desiredProbability);


    const cols = [
      this.getDrugColumn(samplesCount),
      this.getSmilesColumn(samplesCount),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, this.desirabilityColName, desirabilityRaw),
    ];

    this.descriptorStats.forEach((stat, name) => {
      const arr = new Float32Array(samplesCount);

      for (let i = 0; i < samplesCount; ++i) {
        if (desirabilityRaw[i])
          arr[i] = jStat.normal.sample(stat.desAvg, stat.desStd);
        else
          arr[i] = jStat.normal.sample(stat.nonDesAvg, stat.nonDesStd);
      }

      // @ts-ignore
      cols.push(DG.Column.fromFloat32Array(name, arr));
    });

    return DG.DataFrame.fromColumns(cols);
  } // getSyntheticTable

  /** Generates a column with synthetic drug names
   * @param samplesCount Number of samples to generate
   * @returns Column with synthetic drug names */
  private getDrugColumn(samplesCount: number): DG.Column<string> {
    return DG.Column.fromList(
      DG.COLUMN_TYPE.STRING,
      this.drugName,
      Array.from({length: samplesCount}, (_, i) => `${SYNTHETIC_DRUG_NAME} ${i + 1}`));
  }

  /** Generates a column with synthetic SMILES strings
   * @param samplesCount Number of samples to generate
   * @returns Column with synthetic SMILES strings */
  private getSmilesColumn(samplesCount: number): DG.Column<string> {
    return DG.Column.fromList(
      DG.COLUMN_TYPE.STRING,
      this.smilesColName,
      Array.from({length: samplesCount}, () => 'C'));
  }

  /** Clones the source dataframe converting numerical columns to Float type
   * @returns Cloned dataframe */
  private getClonedSourceDfWithFloatNumericCols(): DG.DataFrame {
    const cols: DG.Column[] = [];

    this.sourceDf.columns.toList().forEach((col) => {
      if (col.isNumerical)
        cols.push(col.clone().convertTo(DG.COLUMN_TYPE.FLOAT));
      else
        cols.push(col.clone());
    });

    const clone = DG.DataFrame.fromColumns(cols);
    clone.name = this.sourceDf.name;

    return clone;
  }
} // PmpoDataGenerator
