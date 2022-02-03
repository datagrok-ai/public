export type pubChemSearchType = 'similarity' | 'substructure' | 'identity';
export type pubChemIdType = number | string | boolean;
export type anyObject = {[key: string]: any};
export type paramsType = {[key: string]: pubChemIdType};

export function urlParamsFromObject(obj: {[key: string]: pubChemIdType}): string {
  return Object.entries(obj).map(([key, value]) => `${key}=${value}`).join('&');
}
