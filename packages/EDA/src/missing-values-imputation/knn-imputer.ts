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

/** */
function getNullValue(col: DG.Column): number {
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

/** */
export function simpleImpute(col: DG.Column, strategy: string, inPlace: boolean): DG.Column {
  if (col.stats.missingValueCount === 0)
    throw new Error(ERROR_MSG.NO_MISSING_VALUES);

  const nullValue = getNullValue(col);
  const fillValue = 0;//getFillValue(col, strategy);

  const len = col.length;
  const source = col.getRawData();    

  if (inPlace) {
    for (let i = 0; i < len; ++i)
      if (source[i] === nullValue)
        source[i] = fillValue;

    const buf = col.get(0); 
    col.set(0, col.get(1));
    col.set(0, buf);

    return col;
  }

  //@ts-ignore
  const copy = col.clone();
  copy.name = `${col.name}${COPY_SUFFIX}`;

  const copyRaw = copy.getRawData();

  for (let i = 0; i < len; ++i)
    copyRaw[i] = (source[i] !== nullValue) ? source[i] : fillValue;

  return copy;
}

/** */
export enum METRIC_TYPE {
  ONE_HOT = 'One-hot',
  DIFFERENCE = 'Difference',
};

/** Distance types. */
export enum DISTANCE_TYPE {
  EUCLIDEAN = 'Euclidian',
  MANHATTAN = 'Manhattan',
};

/** Metric specification. */
export type MetricInfo = {
  weight: number,
  type: METRIC_TYPE,
};

/** */
export enum DEFAULT {
  WEIGHT = 1,
  NEIGHBORS = 2,
  IN_PLACE = 1,
};

/** */
export const MIN_NEIGHBORS = 1;

/** */
export function impute(df: DG.DataFrame, targetCols: DG.Column[], featuresMetrics: Map<string, MetricInfo>,
  distance: DISTANCE_TYPE, neighbors: number, inPlace: boolean) 
{
  // 1. Check inputs completness

  if (df.rowCount < 2)
    throw new Error(ERROR_MSG.KNN_NOT_ENOUGH_OF_ROWS);

  if (targetCols.length === 0)
    throw new Error(ERROR_MSG.KNN_NO_TARGET_COLUMNS);

  if (featuresMetrics.size === 0)
    throw new Error(ERROR_MSG.KNN_NO_FEATURE_COLUMNS);
  
  if (featuresMetrics.size === 1)
    targetCols.forEach((col) => {
      if (featuresMetrics.has(col.name))
        throw new Error(`${ERROR_MSG.KNN_NO_FEATURE_COLUMNS} can be used for the column '${col.name}'`);        
      });

  targetCols.forEach((col) => {
    if (!SUPPORTED_COLUMN_TYPES.includes(col.type))
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  });

  featuresMetrics.forEach((val, name) => {
    if (!SUPPORTED_COLUMN_TYPES.includes(df.getCol(name).type))
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  });

  // 2. Missing values imputation by KNN
  targetCols.forEach((col) => {
    const nullValue = getNullValue(col);
    const len = col.length;
    const source = col.getRawData();
    const columns = df.columns;

    const featureSource = [] as Array<Int32Array | Uint32Array | Float32Array | Float64Array>;
    const featureNullVal = [] as number[];
    const metricFunc = [] as ((a: number, b: number) => number)[];

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
    }});
    
    const featuresCount = featureSource.length;
    const properIndices = new Uint32Array(featureSource.length);
    const bufferVector = new Float32Array(featureSource.length);
    let properIndicesCount = 0;

    /** Obtain proper indices for KNN: features with missing vals are skipped */
    const getProperIndeces = (idx: number) => {
      properIndicesCount = 0;
      
      for (let i = 0; i < featuresCount; ++i)
        if (featureSource[i][idx] !== featureNullVal[i]) {
          properIndices[properIndicesCount] = i;
          ++properIndicesCount;
        }
    };

    /** Euclidean distance function */
    const euclideanDistFunc = (vector: Float32Array) => {
      let sum = 0;

      for (let i = 0; i < properIndicesCount; ++i)
        sum += vector[i] * vector[i];

      return Math.sqrt(sum);
    };

    /** Manhattan distance function */
    const manhattanDistFunc = (vector: Float32Array) => {
      let sum = 0;
    
      for (let i = 0; i < properIndicesCount; ++i)
        sum += Math.abs(vector[i]);
    
      return Math.sqrt(sum);
    };

    /** Return norm of the buffer vector (distance between i-th & j-th elements) */
    const dist = (distance === DISTANCE_TYPE.EUCLIDEAN) ? euclideanDistFunc : manhattanDistFunc;

    /** Check if the current item (i.e. table row) can be used */
    const canItemBeUsed = (cur: number) => {
      for (let i = 0; i < properIndicesCount; ++i)
        if (featureSource[i][cur] === featureNullVal[i])
          return false;

      return true;
    };

    /** Get imputation value */
    const getFillValue = (idx: number) => {      
      getProperIndeces(idx);

      // check available features
      if (properIndicesCount === 0) {
        grok.shell.error(`${ERROR_MSG.KNN_IMPOSSIBLE_IMPUTATION}: the column ${col.name}, row ${idx + 1}`);
        return;
      }

      // search for the closest items
      for (let cur = 0; cur < len; ++cur)
        if (canItemBeUsed(cur)) {

        }
      
      if (col.type === DG.COLUMN_TYPE.STRING)
        return "0";
      else 
        return 0;
    };
    
    if (inPlace) {
      for (let i = 0; i < len; ++i)
        if (source[i] === nullValue)
          source[i] = 0;

      const buf = col.get(0);
      col.set(0, col.get(1));
      col.set(1, buf);
    }
    else {
      //@ts-ignore
      const copy = col.clone();

      let i = 1;
      let name= `${col.name}(${COPY_SUFFIX})`;

      // find an appropriate name
      while (df.columns.contains(name)) {
        name = `${col.name}(${COPY_SUFFIX} ${i})`;
        ++i;
      }

      copy.name = name;
      
      const copySource = copy.getRawData();

      for (i = 0; i < len; ++i)
        if (copySource[i] === nullValue)
          copySource[i] = 0;

      df.columns.add(copy);
    }
  });
}
