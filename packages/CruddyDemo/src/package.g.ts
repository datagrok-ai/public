import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: app
//meta.browsePath: Dev
export function northwindDemo() : void {
  PackageFunctions.northwindDemo();
}

//meta.role: app
export function chemblDemo() : void {
  PackageFunctions.chemblDemo();
}
