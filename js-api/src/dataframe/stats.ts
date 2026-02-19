/**
 * Statistics and aggregation classes.
 * @module dataframe/stats
 */

import {AGG, AggregationType} from "../const";
import {toDart} from "../wrappers";
import type {Column} from "./column";
import type {BitSet} from "./bit-set";
import type {DataFrame} from "./data-frame";
import {IDartApi} from "../api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** Represents basic descriptive statistics calculated for a {@link Column}.
 *  See samples: {@link https://public.datagrok.ai/js/samples/data-frame/stats} */
export class Stats {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Calculates statistics for the specified column, optionally filtered by a mask.
   * @param {BitSet} mask */
  static fromColumn(col: Column, mask: BitSet | null = null): Stats {
    return new Stats(api.grok_Stats_FromColumn(col.dart, toDart(mask)));
  }

  /** Calculates statistics for the array of values. */
  static fromValues(values: number[] | Int8Array | Int16Array | Int32Array | Float32Array | Float64Array): Stats {
    return new Stats(api.grok_Stats_FromValues(values));
  }

  /** Total number of values (including missing values). */
  get totalCount(): number {
    return api.grok_Stats_Get_TotalCount(this.dart);
  }

  /** Number of missing (empty) values. */
  get missingValueCount(): number {
    return api.grok_Stats_Get_MissingValueCount(this.dart);
  }

  /** Number of unique values. */
  get uniqueCount(): number {
    return api.grok_Stats_Get_UniqueCount(this.dart);
  }

  /** Number of non-empty values. */
  get valueCount(): number {
    return api.grok_Stats_Get_ValueCount(this.dart);
  }

  /** @returns {number} - minimum */
  get min(): number {
    return api.grok_Stats_Get_Min(this.dart);
  }

  /** @returns {number} - maximum */
  get max(): number {
    return api.grok_Stats_Get_Max(this.dart);
  }

  /** @returns {number} - sum */
  get sum(): number {
    return api.grok_Stats_Get_Sum(this.dart);
  }

  /** @returns {number} - average */
  get avg(): number {
    return api.grok_Stats_Get_Avg(this.dart);
  }

  /** @returns {number} - standard deviation */
  get stdev(): number {
    return api.grok_Stats_Get_Stdev(this.dart);
  }

  /** @returns {number} - variance */
  get variance(): number {
    return api.grok_Stats_Get_Variance(this.dart);
  }

  /** @returns {number} - skewness */
  get skew(): number {
    return api.grok_Stats_Get_Skew(this.dart);
  }

  /** @returns {number} - kurtosis */
  get kurt(): number {
    return api.grok_Stats_Get_Kurt(this.dart);
  }

  /** @returns {number} - median value */
  get med(): number {
    return api.grok_Stats_Get_Med(this.dart);
  }

  /** @returns {number} - first quartile */
  get q1(): number {
    return api.grok_Stats_Get_Q1(this.dart);
  }

  /** @returns {number} - second quartile */
  get q2(): number {
    return api.grok_Stats_Get_Q2(this.dart);
  }

  /** @returns {number} - third quartile */
  get q3(): number {
    return api.grok_Stats_Get_Q3(this.dart);
  }

  /** Pearson correlation between this column and another column. */
  corr(otherColumn: Column): number { return api.grok_Stats_Corr(this.dart, otherColumn.dart); }

  /** Spearman correlation between this column and another column. */
  spearmanCorr(otherColumn: Column): number { return api.grok_Stats_SpearmanCorr(this.dart, otherColumn.dart); }

  /** Returns distributions of [valueColumn] for each category in [catColumn]. */
  static histogramsByCategories(valueColumn: Column, catColumn: Column): Int32Array[] {
    return api.grok_Stats_HistogramsByCategories(valueColumn.dart, catColumn.dart);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}

/**
 * Fluid API for building an aggregation query against a {@link DataFrame}.
 * Build a query by calling the following methods: {@link key}, {@link pivot}, {@link count},
 * {@link uniqueCount}, {@link missingValueCount}, {@link valueCount}, {@link min}, {@link max}, {@link sum},
 * {@link avg}, {@link stdev}, {@link variance}, {@link q1}, {@link q2}, {@link q3}.
 *
 * When the query is constructed, execute it by calling {@link aggregate}, which will
 * produce a {@link DataFrame}.
 *
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 *
 * @example
 * let avgAgesByRaceAndSex = demographicsTable
 *   .groupBy(['race', 'sex'])
 *   .avg('age')
 *   .aggregate();
 */
export class GroupByBuilder {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Performs the aggregation
   *  @returns {DataFrame} */
  aggregate(options?: {autoName?: boolean}): DataFrame {
    // Import dynamically to avoid circular dependency
    const {DataFrame} = require('./data-frame');
    return new DataFrame(api.grok_GroupByBuilder_Aggregate(this.dart, options?.autoName ?? false));
  }

  /**
   * Adds an aggregation to the query.
   * @param {AggregationType} agg - Aggregation type.
   * @param {string} colName - Column name.
   * @param {string} resultColName - Name of the resulting column. Default value is agg(colName).
   * @returns {GroupByBuilder} - this for chaining
   * */
  add(agg: AggregationType, colName?: string | null, resultColName?: string | null): GroupByBuilder {
    api.grok_GroupByBuilder_Add(this.dart, agg, colName, resultColName);
    return this;
  }

  /** Adds a key column to group values on.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  key(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.KEY, srcColName, resultColName);
  }

  /** Adds a column to pivot values on.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  pivot(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.PIVOT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts rows, including these will null values.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  count(resultColName: string = 'count'): GroupByBuilder {
    return this.add(AGG.TOTAL_COUNT, null, resultColName);
  }

  /** Adds an aggregation that counts number of unique values in the specified column.
   * See also {@link count}, {@link valueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  uniqueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.UNIQUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts number of missing values in the specified column.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  missingValueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MISSING_VALUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts rows, including these will null values.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  valueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.VALUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates minimum value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  min(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MIN, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates maximum value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  max(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MAX, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates sum of the values for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * */
  sum(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.SUM, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates median value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} - this for chaining */
  med(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MED, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates average value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * */
  avg(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.AVG, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates standard deviation for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  stdev(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.STDEV, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates varians for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  variance(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.VARIANCE, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates first quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q1(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q1, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates second quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q2(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q2, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates third quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q3(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q3, srcColName, resultColName);
  }

  /** Adds an aggregation that takes first value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  first(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.FIRST, srcColName, resultColName);
  }

  /** Gets groups of DataFrames
   * @returns {Map} - where keys are stings in format 'columnName=value' and values are DataFrames */
  getGroups(): Map<string, DataFrame> {
    return api.grok_GroupByBuilder_GetGroups(this.dart);
  }

  /**
   * Specifies the filter for the source rows.
   * @input {String|Object} pattern
   * @returns {GroupByBuilder}
   **/
  where(pattern: string | object): GroupByBuilder {
    api.grok_GroupByBuilder_Where(this.dart, pattern);
    return this;
  }

  /**
   * */
  whereRowMask(bitset: BitSet): GroupByBuilder {
    if (bitset != null)
      api.grok_GroupByBuilder_WhereBitSet(this.dart, bitset.dart);
    return this;
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}
