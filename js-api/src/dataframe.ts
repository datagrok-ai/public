// noinspection JSUnusedGlobalSymbols

/**
 * Core data structures for the Datagrok platform.
 *
 * This module contains the fundamental classes for working with tabular data:
 *
 * ## Main Classes
 *
 * - **{@link DataFrame}** - High-performance tabular data structure with strongly-typed columns
 * - **{@link Column}** - Strongly-typed column with efficient storage (plus subclasses: FloatColumn, BigIntColumn, DateTimeColumn, ObjectColumn, DataFrameColumn)
 * - **{@link ColumnList}** - Collection of columns in a DataFrame
 * - **{@link Row}** - Single row accessor with property-style column access (uses Proxy)
 * - **{@link RowList}** - Collection of rows with iteration support
 * - **{@link Cell}** - Single cell reference
 *
 * ## Selection and Filtering
 *
 * - **{@link BitSet}** - Efficient bit storage for selection/filter masks
 * - **{@link RowMatcher}** / **{@link ValueMatcher}** - Pattern matching for row/value filtering
 *
 * ## Aggregation and Statistics
 *
 * - **{@link Stats}** - Statistical calculations on columns
 * - **{@link GroupByBuilder}** - Fluent API for grouping and aggregation
 * - **{@link Qnum}** - Qualified numbers with comparison operators
 *
 * ## Helpers
 *
 * - **{@link DataFrameMetaHelper}** - Metadata access (formulaLines, annotationRegions)
 * - **{@link DataFramePlotHelper}** - Quick plot creation
 * - **{@link DataFrameDialogHelper}** - Common dialogs
 * - **{@link ColumnMetaHelper}** - Column metadata
 * - **{@link ColumnColorHelper}** / **{@link ColumnMarkerHelper}** - Column visual styling
 *
 * ## Interfaces
 *
 * - **{@link CsvExportOptions}** - Options for CSV export
 * - **{@link ColumnsFormatCsvExportOptions}** - Column-specific format options
 *
 * @module dataframe
 */

export * from './dataframe/index';
