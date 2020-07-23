/**
 * Cheminformatics support
 */

import {SIMILARITY_METRIC} from './const';
import {BitSet, Column, DataFrame} from './dataframe';

/**
 * Returns molecules similar to the reference one.
 See example: {@link https://public.datagrok.ai/js/samples/domains/chem/similarity-search}
 * @param column - Molecule column to search in.
 * @param molecule - Reference molecule in SMILES format.
 * @param metric - Metric to use.
 * @param limit - Maximum number of results to return.
 * @param minScore - Minimum similarity score for a molecule to be included.
 */
export function similaritySearch(column: Column, molecule: string, metric?: SIMILARITY_METRIC, limit?: number, minScore?: number): Promise<DataFrame>;

/**
 * Returns the specified number of most diverse molecules in the column.
 See example: {@link https://datagrok.ai/help/domains/chem/diversity-search}
 * @param column - Column with molecules in which to search.
 * @param metric - Metric to use.
 * @param limit - Number of molecules to return.
 */
export function diversitySearch(column: Column, metric?: SIMILARITY_METRIC, limit?: number): Promise<DataFrame>;

/**
 * Searches for a molecular pattern in a given column, returning a bitset with hits.
 See example: {@link substructure-search}
 * @param column - Column with molecules to search.
 * @param pattern - Pattern, either SMARTS or SMILES.
 * @param isSmarts - Whether the pattern is SMARTS.
 */
export function substructureSearch(column: Column, pattern: string, isSmarts?: boolean): Promise<BitSet>;

/**
 * Performs R-group analysis.
 See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 * @param table - Table.
 * @param column - Column name with SMILES to analyze.
 * @param core - Core in the SMILES format.
 */
export function rGroup(table: DataFrame, column: string, core: string): Promise<DataFrame>;

/**
 * Finds Most Common Substructure in the specified column.
 See example: {@link https://public.datagrok.ai/js/samples/domains/chem/mcs}
 * @param column - Column with SMILES to analyze.
 */
export function mcs(column: Column): Promise<string>;

/**
 * Calculates specified descriptors for the molecular column.
 See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 * @param table - Table.
 * @param column - Column name with SMILES to calculate descriptors for.
 * @param descriptors - RDKit descriptors to calculate.
 */
export function descriptors(table: DataFrame, column: string, descriptors: string[]): Promise<DataFrame>;

/**
 * Returns available descriptors tree.
 See example: {@link https://public.datagrok.ai/js/samples/domains/chem/descriptors}
 */
export function descriptorsTree(): Promise<object>;

/**
 * Renders a molecule to SVG
 https://public.datagrok.ai/js/samples/domains/chem/mol-rendering
 */
export function svgMol(smiles: string, width?: number, height?: number): HTMLDivElement;

/**
 * Molecule sketcher.
 * @param handler - Molecule on change handler, SMILES
 * @param smiles - Initial molecule
 */
export function sketcher(handler: (...params: any[]) => any, smiles: string): HTMLElement;