/**
 * Type definitions for dataframe module.
 * @module dataframe/types
 */

import type {Row} from './row';
import type {Column} from './column';

export type RowPredicate = (row: Row) => boolean;
export type Comparer = (a: any, b: any) => number;
export type IndexSetter = (index: number, value: any) => void;
export type ColumnId = number | string | Column;
export type MatchType = 'contains' | 'starts with' | 'ends with' | 'equals' | 'regex';

/** Column name -> format */
export interface ColumnsFormatCsvExportOptions {
  [index: string]: string;
}

/** Csv export options to be used in {@link DataFrame.toCsv} */
export interface CsvExportOptions {

  /** Field delimiter; comma if not specified */
  delimiter?: string;

  /** New line character; \n if not specified */
  newLine?: string;

  /** Textual representation of the missing value. Empty string if not specified. */
  missingValue?: string;

  /** Whether the header row containing column names is included. True if not specified. */
  includeHeader?: boolean;

  /** Whether only selected columns are included. False if not specified. */
  selectedColumnsOnly?: boolean;

  /** Whether only visible columns are included. False if not specified. */
  visibleColumnsOnly?: boolean;

  /** Whether only filtered rows are included. Will be combined with [selectedRowsOnly]. */
  filteredRowsOnly?: boolean;

  /** Whether only selected rows are included. Will be combined with [filteredRowsOnly]. */
  selectedRowsOnly?: boolean;

  /** Explicitly defined order of rows. Mutually exclusive with [selectedRowsOnly] or [filteredRowsOnly]. */
  rowIndexes?: number[];

  /** Column order */
  columns?: string[];

  /** Expands qualified numbers into two columns: `qual(column)` and `column` */
  qualifierAsColumn?: boolean;

  /// Saves MOLBLOCKS as SMILES.
  moleculesAsSmiles?: boolean;

  /** Column-specific formats (column name -> format).
      For format examples, see [dateTimeFormatters]. */
  columnFormats?: ColumnsFormatCsvExportOptions;

  /**
   * Whether to include UTF-8 BOM marker at the beginning of the file.
   *   Helps Excel correctly detect UTF-8 encoding for multibyte characters.
   */
  includeUtf8Bom?: boolean;
}

export type GroupDescription = {
  columns?: string[]
  description?: string;
  color?: string;
}

export type GroupsDescription = {
  [index: string]: GroupDescription;
}
