import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Monomers DB Provider
//output: object result
//meta.role: monomerLibProvider
export async function getMonomerDBProvider() : Promise<any> {
  return await PackageFunctions.getMonomerDBProvider();
}

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}
