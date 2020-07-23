/**
 * Machine learning-related routines
 * @module ml
 * */

import {DataFrame} from "./dataframe";

/** Applies predictive model to the specified table.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/predictive-model}
 * @async
 * @param {string} name - Model namespace path.
 * @param {DataFrame} table - Data table.
 * @param {Object} columnNamesMap - Columns map
 * @param {boolean} showProgress - Maximum number of results to return.
 * @returns {Promise<DataFrame>}
 * */
export function applyModel(name: string, table: DataFrame, columnNamesMap?: Object, showProgress?: boolean): Promise<DataFrame>

/** Imputes missing values.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} impute - List of column names to impute missing values.
 * @param {string[]} data - List of column names contains data.
 * @param {number} nearestNeighbours - Number of nearest neighbours.
 * @returns {Promise<DataFrame>}
 * */
export function missingValuesImputation(table: DataFrame, impute: string[], data: string[], nearestNeighbours: number): Promise<DataFrame>

/** Clusters data.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} features - List of column names contains features.
 * @param {number} clusters - Number of clusters.
 * @returns {Promise<DataFrame>}
 * */
export function cluster(table: DataFrame, features: string[], clusters: number): Promise<DataFrame>

/** Principal component analysis.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/pca}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} features - List of column names contains features.
 * @param {number} components - Number of clusters.
 * @param {boolean} center - Center features data before PCA.
 * @param {boolean} scale - Scale features data before PCA.
 * @returns {Promise<DataFrame>}
 * */
export function pca(table: DataFrame, features: string[], components: number, center: boolean, scale: boolean): DataFrame

/** Creates a table with random values from the specified distribution.
 * Documentation: {@link https://datagrok.ai/help/transform/random-data}
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/random-data}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string} distribution - Distribution name.
 * @param {Object} params - Distribution parameters.
 * @param {number} seed - Initial seed.
 * @returns {Promise<DataFrame>}
 * */
export function randomData(table: DataFrame, distribution: string, params: Object, seed: number): Promise<DataFrame>