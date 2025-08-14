import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() {
  return PackageFunctions.info();
}

//name: expectTable
//input: dataframe actual 
//input: dataframe expected 
//output: bool result
export function expectTable(actual: DG.DataFrame, expected: DG.DataFrame) {
  return PackageFunctions.expectTable(actual, expected);
}
