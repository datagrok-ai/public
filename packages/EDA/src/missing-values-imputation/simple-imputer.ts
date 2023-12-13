import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMPUTATION_STRATEGY, ERROR_MSG, COPY_SUFFIX} from './constants';

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
function getFillValue(col: DG.Column, strategy: string): number {
  let res: number;

  switch (strategy) {
    case IMPUTATION_STRATEGY.MEAN:        
      res = col.stats.avg;
      break;
  
    case IMPUTATION_STRATEGY.MEDIAN:
      res = col.stats.med;
      break;
    
    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_IMPUTATION_STRATEGY);
  }

  return ((col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.STRING))? Math.round(res) : res;
}

/** */
function getTypedArrayOfTheSameSize(col: DG.Column): Int32Array | Float32Array | Float64Array {
  const size = col.getRawData().length;

  switch (col.type) {
    case DG.COLUMN_TYPE.INT:
    case DG.COLUMN_TYPE.STRING:  
      return new Int32Array(size);

    case DG.COLUMN_TYPE.FLOAT:
      return new Float32Array(size);

    case DG.COLUMN_TYPE.DATE_TIME:
    case DG.COLUMN_TYPE.QNUM:
      return new Float64Array(size);
  
    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** */
export function simpleImpute(col: DG.Column, strategy: string, inPlace: boolean, toMark: boolean): DG.Column {
  if (col.stats.missingValueCount === 0)
    throw new Error(ERROR_MSG.NO_MISSING_VALUES);

  const nullValue = getNullValue(col);
  const fillValue = getFillValue(col, strategy);

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
