import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Monomers DB Provider
//tags: monomer-lib-provider
//output: object result
export async function getMonomerDBProvider() : Promise<any> {
  return await PackageFunctions.getMonomerDBProvider();
}

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}
