import * as DG from 'datagrok-api/dg';
import * as moment from 'moment';
import { normalizeStrings } from './utils';

export const isNullVariable = (val: any) => { return val.column.isNone(val.index); };

export const startAfterEnd = (values: any[]) => { return values[1] < values[0] && values[1] !== DG.INT_NULL; }

export const negativeValue = (val: any[]) => { return val[0] <= 0 && val[0] !== DG.INT_NULL && val[0] !== ''; }

export const notISO8601 = (val: any[]) => { return !moment(val[0], moment.ISO_8601).isValid() && val[0] !== ''; }

export const seriousnessCriteriaNotIndicated = (val: any) => { 
  const isSerious = val[0];
  var criteria = val.filter(it => it !== isSerious);
  return !(isSerious === 'Y' && criteria.some(item => item === 'Y'));
 }

 export const identicalValues = (values: any[]) => {
   const valuesToCompare = normalizeStrings(values);
   return valuesToCompare[0] === valuesToCompare[1]; 
  }

 export const nonIdenticalValuesWithExternal = (values: any[], valueToCompare) => { 
   const valuesToCompare = normalizeStrings(values.concat(valueToCompare));
   return valuesToCompare[0] !== valuesToCompare[1]; 
  }

