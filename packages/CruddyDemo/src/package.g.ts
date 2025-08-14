import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: northwindDemo
//tags: app
//meta.browsePath: Dev
export function northwindDemo() {
  return PackageFunctions.northwindDemo();
}

//name: chemblDemo
//tags: app
export function chemblDemo() {
  return PackageFunctions.chemblDemo();
}
