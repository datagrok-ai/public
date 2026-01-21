/**
 * Machine learning-related routines
 * @module ml
 * */

import {DataFrame} from "./dataframe";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export namespace ml {
  /** Applies predictive model to the specified table.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/predictive-model}
   * @async
   * @param {string} name - Model namespace path.
   * @param {DataFrame} table - Data table.
   * @param {object} columnNamesMap - Columns map.
   * @param {boolean} showProgress - Maximum number of results to return.
   * */
  export async function applyModel(name: string, table: DataFrame, columnNamesMap: object = {}, showProgress: boolean = true): Promise<DataFrame> {
    await api.grok_ML_ApplyModel(name, table.dart, columnNamesMap, showProgress);
    return table;
  }

  /** Imputes missing values.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation}
   * @async
   * @param {DataFrame} table - Data table.
   * @param {string[]} impute - List of column names to impute missing values.
   * @param {string[]} data - List of column names containing data.
   * @param {number} nearestNeighbours - Number of nearest neighbours.
   * */
  export async function missingValuesImputation(table: DataFrame, impute: string[], data: string[], nearestNeighbours: number): Promise<DataFrame> {
    await api.grok_ML_MissingValuesImputation(table.dart, impute, data, nearestNeighbours);
    return table;
  }

  /** Clusters data.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
   * @async
   * @param {DataFrame} table - Data table.
   * @param {string[]} features - List of column names containing features.
   * @param {number} clusters - Number of clusters.
   * */
  export async function cluster(table: DataFrame, features: string[], clusters: number): Promise<DataFrame> {
    await api.grok_ML_Cluster(table.dart, features, clusters);
    return table;
  }

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
  export async function pca(table: DataFrame, features: string[], components: number, center: boolean, scale: boolean): Promise<DataFrame> {
    await api.grok_ML_PCA(table.dart, features, components, center, scale);
    return table;
  }

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
  export async function randomData(table: DataFrame, distribution: string, params: object, seed: number): Promise<DataFrame> {
    await api.grok_ML_RandomData(table.dart, distribution, params, seed);
    return table;
  }
}
