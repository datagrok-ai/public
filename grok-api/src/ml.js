import {DataFrame} from "./dataframe";

/** Machine learning-related routines */
export class ml {
    /** Applies predictive modesl to the specified table. */
    static applyModel(name, table, columnNamesMap, showProgress) {
        if (columnNamesMap === undefined)
            columnNamesMap = new Map();
        return new Promise((resolve, reject) =>
            grok_ML_ApplyModel(name, table.d, (t) => resolve(new DataFrame(t)), columnNamesMap, showProgress));
    }

    /** Imputes missing values. */
    static missingValuesImputation(table, impute, data, nearestNeighbours) {
        return new Promise((resolve, reject) =>
            grok_ML_MissingValuesImputation(table.d, impute, data, nearestNeighbours, () => resolve(table)));
    }

    /** Clusters data. */
    static cluster(table, features, clusters) {
        return new Promise((resolve, reject) =>
            grok_ML_Cluster(table.d, features, clusters, () => resolve(table)));
    }

    /** Principal component analysis. */
    static pca(table, features, components, center, scale) {
        return new Promise((resolve, reject) =>
            grok_ML_PCA(table.d, features, components, center, scale, () => resolve(table)));
    }

    /** Creates a table with random values from the specified distribution. */
    static randomData(table, distribution, params, seed) {
        return new Promise((resolve, reject) =>
            grok_ML_RandomData(table.d, distribution, params, seed, () => resolve(table)));
    }
}
