import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: init
//tags: init
//output: dynamic result
//top-menu: hello
export async function init() {
  return PackageFunctions.init();
}
