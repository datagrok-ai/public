import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Monomers Domain DB Provider
//output: object result
//meta.role: monomer-lib-provider
export async function getMonomerDomainDbProvider() : Promise<any> {
  return await PackageFunctions.getMonomerDomainDbProvider();
}
