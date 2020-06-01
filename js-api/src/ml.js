/** Machine learning-related routines
 * @module ml
*/
import {DataFrame} from "./dataframe";


/** Applies predictive model to the specified table.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/predictive-model}
 * @async
 * @param {string} name - Model namespace path.
 * @param {DataFrame} table - Data table.
 * @param {Object} columnNamesMap - Columns map
 * @param {boolean} showProgress - Maximum number of results to return.
 * @returns {DataFrame}
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
 * @returns {DataFrame}
 * */
export async function missingValuesImputation(table, impute, data, nearestNeighbours) {
    return await grok_ML_MissingValuesImputation(table.d, impute, data, nearestNeighbours);
}

/** Clusters data.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} features - List of column names contains features.
 * @param {number} clusters - Number of clusters.
 * @returns {DataFrame}
 * */
export async function cluster(table, features, clusters) {
    return await grok_ML_Cluster(table.d, features, clusters);
}

/** Principal component analysis.
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/pca}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string[]} features - List of column names contains features.
 * @param {number} components - Number of clusters.
 * @param {boolean} center - Center features data before PCA.
 * @param {boolean} scale - Scale features data before PCA.
 * @returns {DataFrame}
 * */
export async function pca(table, features, components, center, scale) {
    return await grok_ML_PCA(table.d, features, components, center, scale);
}

/** Creates a table with random values from the specified distribution.
 * Documentation: {@link https://datagrok.ai/help/transform/random-data}
 * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/random-data}
 * @async
 * @param {DataFrame} table - Data table.
 * @param {string} distribution - Distribution name.
 * @param {Object} params - Distribution parameters.
 * @param {number} seed - Initial seed.
 * @returns {DataFrame}
 * */
export async function randomData(table, distribution, params, seed) {
    return await grok_ML_RandomData(table.d, distribution, params, seed);
}
