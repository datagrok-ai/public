import * as DG from 'datagrok-api/dg';

import {sortByReverseLength} from '../helpers';

export function getSaltMass(
  saltNames: string[], saltsMolWeightList: number[], equivalentsCol: DG.Column, i: number, saltCol: DG.Column
): number {
  const saltRowIndex = saltNames.indexOf(saltCol.get(i));
  return (
    saltRowIndex == -1 || saltsMolWeightList[saltRowIndex] == DG.FLOAT_NULL || equivalentsCol.get(i) == DG.INT_NULL
  ) ?
    DG.FLOAT_NULL :
    saltsMolWeightList[saltRowIndex] * equivalentsCol.get(i);
}

export function getSaltMolWeigth(
  saltNamesList: string[], saltCol: DG.Column, saltsMolWeightList: number[], i: number
): number {
  const saltRowIndex = saltNamesList.indexOf(saltCol.get(i));
  return (saltRowIndex == -1) ? DG.FLOAT_NULL : saltsMolWeightList[saltRowIndex];
}

export function getBatchMolWeight(compoundMolWeightCol: DG.Column, saltMassCol: DG.Column, i: number): number {
  return (compoundMolWeightCol.getString(i) == '' || saltMassCol.getString(i) == '') ?
    DG.FLOAT_NULL :
    compoundMolWeightCol.get(i) + saltMassCol.get(i);
}

export function getMolWeight(sequence: string, codesToWeightsMap: Map<string, number>): number {
  const codes = sortByReverseLength(Array.from(codesToWeightsMap.keys()));
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    const code = sequence.slice(i, i + matchedCode.length);
    const value = codesToWeightsMap.get(code);
    weight += (value === undefined) ? 0 : value;
    i += matchedCode.length;
  }
  return weight - 61.97;
}
