/**
 * Machine learning-related routines
 * @module ml
 * */
define(['require', 'exports', './dataframe'], function (require, exports, dataframe_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.randomData = exports.pca = exports.cluster = exports.missingValuesImputation = exports.applyModel = void 0;
  let api = window;
  /** Applies predictive model to the specified table.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/predictive-model}
     * @async
     * @param {string} name - Model namespace path.
     * @param {DataFrame} table - Data table.
     * @param {object} columnNamesMap - Columns map.
     * @param {boolean} showProgress - Maximum number of results to return.
     * @returns {Promise<DataFrame>}
     * */
  function applyModel(name, table, columnNamesMap = {}, showProgress = true) {
    return new Promise((resolve, reject) => api.grok_ML_ApplyModel(name, table.d, (t) => resolve(new dataframe_1.DataFrame(t)), (e) => reject(e), columnNamesMap, showProgress));
  }
  exports.applyModel = applyModel;
  /** Imputes missing values.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation}
     * @async
     * @param {DataFrame} table - Data table.
     * @param {string[]} impute - List of column names to impute missing values.
     * @param {string[]} data - List of column names containing data.
     * @param {number} nearestNeighbours - Number of nearest neighbours.
     * @returns {Promise<DataFrame>}
     * */
  function missingValuesImputation(table, impute, data, nearestNeighbours) {
    return new Promise((resolve, reject) => api.grok_ML_MissingValuesImputation(table.d, impute, data, nearestNeighbours, () => resolve(table), (e) => reject(e)));
  }
  exports.missingValuesImputation = missingValuesImputation;
  /** Clusters data.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
     * @async
     * @param {DataFrame} table - Data table.
     * @param {string[]} features - List of column names containing features.
     * @param {number} clusters - Number of clusters.
     * @returns {Promise<DataFrame>}
     * */
  function cluster(table, features, clusters) {
    return new Promise((resolve, reject) => api.grok_ML_Cluster(table.d, features, clusters, () => resolve(table), (e) => reject(e)));
  }
  exports.cluster = cluster;
  /** Principal component analysis.
     * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/pca}
     * @async
     * @param {DataFrame} table - Data table.
     * @param {string[]} features - List of column names containing features.
     * @param {number} components - Number of clusters.
     * @param {boolean} center - Center features data before PCA.
     * @param {boolean} scale - Scale features data before PCA.
     * @returns {Promise<DataFrame>}
     * */
  function pca(table, features, components, center, scale) {
    return new Promise((resolve, reject) => api.grok_ML_PCA(table.d, features, components, center, scale, () => resolve(table), (e) => reject(e)));
  }
  exports.pca = pca;
  /** Creates a table with random values from the specified distribution.
     * Documentation: {@link https://datagrok.ai/help/transform/random-data}
     * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/random-data}
     * @async
     * @param {DataFrame} table - Data table.
     * @param {string} distribution - Distribution name.
     * @param {object} params - Distribution parameters.
     * @param {number} seed - Initial seed.
     * @returns {Promise<DataFrame>}
     * */
  function randomData(table, distribution, params, seed) {
    return new Promise((resolve, reject) => api.grok_ML_RandomData(table.d, distribution, params, seed, () => resolve(table), (e) => reject(e)));
  }
  exports.randomData = randomData;
});
