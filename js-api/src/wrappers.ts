declare let DG: any;

export function paramsToJs(params: any): any {
  return DG.paramsToJs(params);
}

export function toJs(dart: any, check: boolean = false): any {
  return DG.toJs(dart);
}

export function toDart(x: any): any {
  return DG.toDart(x);
}
