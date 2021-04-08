declare let DG: any;

export function paramsToJs(params: any): any {
  return DG.paramsToJs(params);
}

export function toJs(d: any, check: boolean = false): any {
  return DG.toJs(d);
}

export function toDart(x: any): any {
  return DG.toDart(x);
}
