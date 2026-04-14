/**
 * Core types for the time-series feature extraction module.
 */

/** Accepted typed numeric array types for time-series data. */
export type NumericArray =
  | Int32Array<ArrayBufferLike>
  | Uint32Array<ArrayBufferLike>
  | Float32Array<ArrayBufferLike>
  | Float64Array<ArrayBufferLike>;

/** Single value column within a time-series DataFrame. */
export interface TimeSeriesColumn {
  /** Column name, used as prefix in feature naming (e.g. "pH"). */
  name: string;
  /** Raw data values; length must equal the parent DataFrame's rowCount. */
  data: NumericArray;
}

/**
 * Columnar input format for feature extraction.
 *
 * All rows for the same `id` must be contiguous (sorted by id).
 * Multiple value columns can be provided — each produces its own feature set.
 */
export interface TimeSeriesDataFrame {
  /** Sample identifier for each row. Groups rows into individual time series. */
  ids: Uint32Array<ArrayBufferLike>;
  /** Time axis values for each row. */
  time: NumericArray;
  /** Total number of rows across all samples. */
  rowCount: number;
  /** Value columns to extract features from. */
  columns: TimeSeriesColumn[];
}

/** Single feature column in the output matrix. */
export interface FeatureColumn {
  /** Feature name following tsfresh convention (e.g. "pH__mean"). */
  name: string;
  /** One value per sample. */
  data: Float64Array;
}

/** Output of feature extraction: one row per sample, one column per feature. */
export interface FeatureMatrix {
  /** Unique sample ids, one per row, in extraction order. */
  sampleIds: Uint32Array<ArrayBufferLike>;
  /** Feature columns; each has `nSamples` entries. */
  columns: FeatureColumn[];
  /** Number of samples (rows). */
  nSamples: number;
}

/** Options for {@link extractFeatures}. */
export interface ExtractOptions {
  /** If true, validate input for NaN, Inf, length mismatches, and contiguous ids (default: false). */
  validate?: boolean;
}
