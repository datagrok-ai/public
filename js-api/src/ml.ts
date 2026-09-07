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
   * @param name - Model namespace path.
   * @param table - Data table.
   * @param columnNamesMap - Columns map.
   * @param showProgress - Whether to show a progress indicator. */
  export async function applyModel(name: string, table: DataFrame, columnNamesMap: object = {}, showProgress: boolean = true): Promise<DataFrame> {
    await api.grok_ML_ApplyModel(name, table.dart, columnNamesMap, showProgress);
    return table;
  }

  /** Imputes missing values.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/missing-values-imputation}
   * @param table - Data table.
   * @param impute - List of column names to impute missing values.
   * @param data - List of column names containing data.
   * @param nearestNeighbours - Number of nearest neighbours. */
  export async function missingValuesImputation(table: DataFrame, impute: string[], data: string[], nearestNeighbours: number): Promise<DataFrame> {
    await api.grok_ML_MissingValuesImputation(table.dart, impute, data, nearestNeighbours);
    return table;
  }

  /** Clusters data.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/cluster}
   * @param table - Data table.
   * @param features - List of column names containing features.
   * @param clusters - Number of clusters. */
  export async function cluster(table: DataFrame, features: string[], clusters: number): Promise<DataFrame> {
    await api.grok_ML_Cluster(table.dart, features, clusters);
    return table;
  }

  /** Principal component analysis.
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/pca}
   * @param table - Data table.
   * @param features - List of column names containing features.
   * @param components - Number of clusters.
   * @param center - Center features data before PCA.
   * @param scale - Scale features data before PCA. */
  export async function pca(table: DataFrame, features: string[], components: number, center: boolean, scale: boolean): Promise<DataFrame> {
    await api.grok_ML_PCA(table.dart, features, components, center, scale);
    return table;
  }

  /** Creates a table with random values from the specified distribution.
   * Documentation: {@link https://datagrok.ai/help/transform/random-data}
   * See example: {@link https://public.datagrok.ai/js/samples/domains/data-science/random-data}
   * @param table - Data table.
   * @param distribution - Distribution name.
   * @param params - Distribution parameters.
   * @param seed - Initial seed. */
  export async function randomData(table: DataFrame, distribution: string, params: object, seed: number): Promise<DataFrame> {
    await api.grok_ML_RandomData(table.dart, distribution, params, seed);
    return table;
  }
}
