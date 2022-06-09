import * as DG from 'datagrok-api/dg';

export function isNullValue(obValue : any) : boolean {
  if(typeof obValue === "number" && isNaN(obValue)) {
    return true;
  }

  return obValue === undefined || obValue === null || obValue === "" || obValue === DG.FLOAT_NULL || obValue === DG.INT_NULL;
}

export function isStrictInt(val : number) : boolean {
  const b = Math.floor(val) === val;
  return b;
}
