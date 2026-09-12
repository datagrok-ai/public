import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: App
//meta.browsePath: Dev
export function northwindDemo() : void {
  PackageFunctions.northwindDemo();
}

//meta.role: App
export function chemblDemo() : void {
  PackageFunctions.chemblDemo();
}
