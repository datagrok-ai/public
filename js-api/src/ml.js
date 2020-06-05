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
export function applyModel(name, table, columnNamesMap = {}, showProgress = true) {
    return new Promise((resolve, reject) =>
        grok_ML_ApplyModel(name, table.d, (t) => resolve(new DataFrame(t)), columnNamesMap, showProgress));
}

/** Imputes missing values.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} impute - List of column names to impute missing values.
 * @param {string[]} data - List of column names contains data.
 * @param {number} nearestNeighbours - Number of nearest neighbours.
 * @returns {Promise<DataFrame>}
 * */
export function missingValuesImputation(table, impute, data, nearestNeighbours) {
    return new Promise((resolve, reject) =>
        grok_ML_MissingValuesImputation(table.d, impute, data, nearestNeighbours, () => resolve(table)));
}

/** Clusters data.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} features - List of column names contains features.
 * @param {number} clusters - Number of clusters.
 * @returns {Promise<DataFrame>}
 * */
export function cluster(table, features, clusters) {
    return new Promise((resolve, reject) =>
        grok_ML_Cluster(table.d, features, clusters, () => resolve(table)));
}

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
export function pca(table, features, components, center, scale) {
    return new Promise((resolve, reject) =>
        grok_ML_PCA(table.d, features, components, center, scale, () => resolve(table)));
}

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
export function randomData(table, distribution, params, seed) {
    return new Promise((resolve, reject) =>
        grok_ML_RandomData(table.d, distribution, params, seed, () => resolve(table)));
}
