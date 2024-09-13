// Tools for missing values imputation using the k-nearest neighbors method

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ERROR_MSG, COPY_SUFFIX} from './ui-constants';

/** Column types supported by the missing values imputer */
export const SUPPORTED_COLUMN_TYPES = [
  DG.COLUMN_TYPE.INT,
  DG.COLUMN_TYPE.FLOAT,
  DG.COLUMN_TYPE.STRING,
  DG.COLUMN_TYPE.DATE_TIME,
  DG.COLUMN_TYPE.QNUM,
] as string[];

/** Return null value with respect to the column type */
export function getNullValue(col: DG.Column): number {
  switch (col.type) {
  case DG.COLUMN_TYPE.INT:
    return DG.INT_NULL;

  case DG.COLUMN_TYPE.FLOAT:
    return DG.FLOAT_NULL;

  case DG.COLUMN_TYPE.QNUM:
    return DG.FLOAT_NULL;

  case DG.COLUMN_TYPE.DATE_TIME:
    return DG.FLOAT_NULL;

  case DG.COLUMN_TYPE.STRING:
    return col.max;

  default:
    throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** Metric types (between column elements) */
export enum METRIC_TYPE {
  ONE_HOT = 'One-hot',
  DIFFERENCE = 'Difference',
};

/** Distance types (over several columns). */
export enum DISTANCE_TYPE {
  EUCLIDEAN = 'Euclidean',
  MANHATTAN = 'Manhattan',
};

/** Metric specification. */
export type MetricInfo = {
  weight: number,
  type: METRIC_TYPE,
};

/** Default values */
export enum DEFAULT {
  WEIGHT = 1,
  NEIGHBORS = 4,
  IN_PLACE = 1,
  SELECTED = 1,
  KEEP_EMPTY = 0,
};

/** Min number of neighbors for KNN */
export const MIN_NEIGHBORS = 1;

/** Dataframe item: index - number of row,  dist - distance to the target element */
type Item = {
  index: number,
  dist: number,
};

/** Impute missing values using the KNN method and returns an array of items for which an imputation fails */
export function impute(df: DG.DataFrame, targetColNames: string[], featuresMetrics: Map<string, MetricInfo>,
  missingValsIndices: Map<string, number[]>, distance: DISTANCE_TYPE, neighbors: number,
  inPlace: boolean): Map<string, number[]> {
  // 1. Check inputs completness

  if (neighbors < MIN_NEIGHBORS)
    throw new Error(ERROR_MSG.INCORRECT_NEIGHBORS);

  if (df.rowCount < 2)
    throw new Error(ERROR_MSG.KNN_NOT_ENOUGH_OF_ROWS);

  if (targetColNames.length === 0)
    throw new Error(ERROR_MSG.KNN_NO_TARGET_COLUMNS);

  if (featuresMetrics.size === 0)
    throw new Error(ERROR_MSG.KNN_NO_FEATURE_COLUMNS);

  if (featuresMetrics.size === 1) {
    targetColNames.forEach((name) => {
      if (featuresMetrics.has(name))
        throw new Error(`${ERROR_MSG.KNN_NO_FEATURE_COLUMNS} can be used for the column '${name}'`);
    });
  }

  targetColNames.forEach((name) => {
    if (!missingValsIndices.has(name))
      throw new Error(`${ERROR_MSG.KNN_FAILS}: ${ERROR_MSG.WRONG_PREDICTIONS}`);
  });

  const columns = df.columns;

  // 2. Imputation

  targetColNames.forEach((name) => {
    if (!SUPPORTED_COLUMN_TYPES.includes(columns.byName(name).type))
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  });

  featuresMetrics.forEach((val, name) => {
    if (!SUPPORTED_COLUMN_TYPES.includes(df.getCol(name).type))
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  });

  /** Failed to impute items */
  const failedToImpute = new Map<string, number[]>();

  // 2. Missing values imputation in each target column
  targetColNames.forEach((name) => {
    const col = columns.byName(name);
    const nullValue = getNullValue(col);
    const len = col.length;
    const source = col.getRawData();
    const frequencies = new Uint16Array(col.categories.length);

    const featureSource = [] as Array<Int32Array | Uint32Array | Float32Array | Float64Array>;
    const featureNullVal = [] as number[];
    const metricFunc = [] as ((a: number, b: number) => number)[];

    const failedToImputeIndices = [] as number[];

    // create features tools
    featuresMetrics.forEach((metricInfo, name) => {
      if (name !== col.name) {
        const feature = columns.byName(name);
        featureSource.push(feature.getRawData());
        featureNullVal.push(getNullValue(feature));

        switch (metricInfo.type) {
        case METRIC_TYPE.DIFFERENCE:
          metricFunc.push((a: number, b: number) => metricInfo.weight * Math.abs(a - b));
          break;

        case METRIC_TYPE.ONE_HOT:
          metricFunc.push((a: number, b: number) => metricInfo.weight * ((a === b) ? 0 : 1));
          break;

        default:
          break;
        }
      }
    });

    const featuresCount = featureSource.length;
    const properIndices = new Uint32Array(featureSource.length);
    const bufferVector = new Float32Array(featureSource.length);
    let properIndicesCount = 0;

    // closest items
    const nearestItems = new Array<Item>(neighbors);
    let nearestItemsCount = 0;

    // auxiliry variables
    let maxInd = 0;
    let maxDist = 0;
    let sum = 0;
    let fillValue = 0;

    /** Obtain proper indices for KNN: features with missing vals are skipped */
    const getProperIndeces = (idx: number) => {
      properIndicesCount = 0;

      for (let i = 0; i < featuresCount; ++i) {
        if (featureSource[i][idx] !== featureNullVal[i]) {
          properIndices[properIndicesCount] = i;
          ++properIndicesCount;
        }
      }
    };

    /** Compute buffer vector */
    const computeBufferVector = (idx: number, cur: number) => {
      properIndices.forEach((properIndex, k) => {
        bufferVector[k] = metricFunc[properIndex](featureSource[properIndex][idx], featureSource[properIndex][cur]);
      });
    };

    /** Euclidean distance function */
    const euclideanDistFunc = () => {
      let sum = 0;

      for (let i = 0; i < properIndicesCount; ++i)
        sum +=bufferVector[i] * bufferVector[i];

      return Math.sqrt(sum);
    };

    /** Manhattan distance function */
    const manhattanDistFunc = () => {
      let sum = 0;

      for (let i = 0; i < properIndicesCount; ++i)
        sum += Math.abs(bufferVector[i]);

      return Math.sqrt(sum);
    };

    /** Return norm of the buffer vector (distance between i-th & j-th elements) */
    const dist = (distance === DISTANCE_TYPE.EUCLIDEAN) ? euclideanDistFunc : manhattanDistFunc;

    /** Check if the current item (i.e. table row) can be used */
    const canItemBeUsed = (cur: number) => {
      if (source[cur] === nullValue)
        return false;

      for (let i = 0; i < properIndicesCount; ++i) {
        if (featureSource[properIndices[i]][cur] === featureNullVal[properIndices[i]])
          return false;
      }

      return true;
    };

    /** Return the most frequent of the nearest items (for categorial data) */
    const mostFrequentOfTheNearestItems = () => {
      frequencies.forEach((v, i, arr) => arr[i] = 0);
      let i = 0;

      for (i = 0; i < nearestItemsCount; ++i)
        ++frequencies[source[nearestItems[i].index]];

      let maxFreq = frequencies[0];
      let maxFreqIdx = 0;

      frequencies.forEach((v, i) => {
        if (v > maxFreq) {
          maxFreq = v;
          maxFreqIdx = i;
        }
      });

      return maxFreqIdx;
    };

    /** Get imputation value */
    const getFillValue = (idx: number) => {
      getProperIndeces(idx);

      // check available features
      if (properIndicesCount === 0)
        throw new Error(`${ERROR_MSG.KNN_IMPOSSIBLE_IMPUTATION}: the column "${col.name}", row ${idx + 1}`);

      nearestItemsCount = 0;

      // search for the closest items
      for (let cur = 0; cur < len; ++cur) {
        if (canItemBeUsed(cur) && (cur !== idx)) {
        // 1) compute distance between cur-th and idx-th items
          computeBufferVector(idx, cur);
          const curDist = dist();

          // 2) insert the current item
          if (nearestItemsCount < neighbors) {
            nearestItems[nearestItemsCount] = {index: cur, dist: curDist};
            ++nearestItemsCount;
          } else {
          // 2.1) find the farest
            maxInd = 0;
            maxDist = nearestItems[0].dist;

            for (let i = 1; i < nearestItemsCount; ++i) {
              if (maxDist < nearestItems[i].dist) {
                maxDist = nearestItems[i].dist;
                maxInd = i;
              }
            }

            // 2.2) replace
            if (curDist < maxDist)
              nearestItems[maxInd] = {index: cur, dist: curDist};
          } // else
        }
      } // for cur

      // check found nearest items
      if (nearestItemsCount === 0)
        throw new Error(`${ERROR_MSG.KNN_IMPOSSIBLE_IMPUTATION}: the column "${col.name}", row ${idx + 1}`);

      if (col.type === DG.COLUMN_TYPE.STRING)
        return mostFrequentOfTheNearestItems();

      // compute fill value
      sum = 0;
      for (let i = 0; i < nearestItemsCount; ++i)
        sum += source[nearestItems[i].index];

      fillValue = sum / nearestItemsCount;

      if (col.type === DG.COLUMN_TYPE.INT)
        return Math.round(fillValue);

      return fillValue;
    }; // getFillValue

    if (inPlace) {
      // use indices found previousely
      for (const i of missingValsIndices.get(name)!) {
        try {
          source[i] = getFillValue(i);
        } catch (err) {
          failedToImputeIndices.push(i);

          if (!(err instanceof Error))
            grok.shell.error(ERROR_MSG.CORE_ISSUE);
        }
      }

      if (failedToImputeIndices.length > 0)
        failedToImpute.set(name, failedToImputeIndices);

      // to reset view
      col.set(0, col.get(0));
    } else {
      //@ts-ignore
      const copy = col.clone();

      let i = 1;
      let copyName = `${name}(${COPY_SUFFIX})`;

      // find an appropriate name
      while (df.columns.contains(copyName)) {
        copyName = `${name}(${COPY_SUFFIX} ${i})`;
        ++i;
      }

      copy.name = copyName;

      const copySource = copy.getRawData();

      // use indices found previousely
      for (const i of missingValsIndices.get(name)!) {
        try {
          copySource[i] = getFillValue(i);
        } catch (err) {
          failedToImputeIndices.push(i);

          if (!(err instanceof Error))
            grok.shell.error(ERROR_MSG.CORE_ISSUE);
        }
      }

      if (failedToImputeIndices.length > 0)
        failedToImpute.set(copyName, failedToImputeIndices);

      copy.set(0, copy.get(0));

      df.columns.add(copy);
    } // else
  });

  return failedToImpute;
} // impute

/** Return indices of missing values for each column */
export function getMissingValsIndices(columns: DG.Column[]): Map<string, number[]> {
  const misValsInds = new Map<string, number[]>();

  for (const col of columns) {
    if (!SUPPORTED_COLUMN_TYPES.includes(col.type))
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);

    if (col.stats.missingValueCount === 0)
      continue;

    const indices = [] as number[];
    const nullValue = getNullValue(col);

    col.getRawData().forEach((val, idx) => {
      if (val === nullValue)
        indices.push(idx);
    });

    misValsInds.set(col.name, indices);
  }

  return misValsInds;
}

/** Predict existence of missing values imputation fails */
export function areThereFails(targetColNames: string[], featureColNames: string[],
  misValsInds: Map<string, number[]>): boolean {
  // check feature columns
  for (const name of featureColNames) {
    if (!misValsInds.has(name))
      return false;
  }

  // check target columns
  for (const target of targetColNames) {
    const indices = misValsInds.get(target);

    if (indices === undefined)
      throw new Error(ERROR_MSG.FAILS_TO_PREDICT_IMPUTATION_FAILS);

    for (const idx of indices) {
      let failToImpute = true;

      for (const feature of featureColNames) {
        const featureInds = misValsInds.get(feature);

        if (featureInds === undefined)
          throw new Error(ERROR_MSG.FAILS_TO_PREDICT_IMPUTATION_FAILS);

        if (!featureInds.includes(idx)) {
          failToImpute = false;
          break;
        }
      }

      if (failToImpute)
        return true;
    }
  }

  return false;
} // predictFails

/** Returns first non-null value */
function getFirstNonNull<T>(col: DG.Column<T>): T {
  const nullValue = getNullValue(col);
  const raw = col.getRawData();
  const len = raw.length;

  for (let i = 0; i < len; ++i) {
    if (raw[i] !== nullValue)
      return col.get(i)!;
  }

  throw new Error(ERROR_MSG.EMPTY_COLUMN);
}

/** Return default fill value with respect to the column type */
function getDefaultFillValue<T>(col: DG.Column<T>): T {
  switch (col.type) {
  case DG.COLUMN_TYPE.STRING:
  case DG.COLUMN_TYPE.DATE_TIME:
    return getFirstNonNull(col); // TODO: replace by most frequent

  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.QNUM:
    return col.stats.avg as T;

  default:
    throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** Perform missing values imputation using the simple approach */
export function imputeFailed(df: DG.DataFrame, failedToImpute: Map<string, number[]>): void {
  failedToImpute.forEach((indices, colName) => {
    const col = df.col(colName);
    if (col !== null) {
      if (!SUPPORTED_COLUMN_TYPES.includes(col.type))
        throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);

      const fillVal = getDefaultFillValue(col);
      indices.forEach((idx) => col.set(idx, fillVal));
    }
  });
}
