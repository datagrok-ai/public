import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Consensus Pharmacophore
//output: view result
//meta.role: app
//meta.browsePath: Bio
export async function consensusPharmacophoreApp() : Promise<any> {
  return await PackageFunctions.consensusPharmacophoreApp();
}

//description: Demo: Consensus Pharmacophore on 5 EGFR kinase structures.
export async function consensusPharmacophoreDemo() : Promise<void> {
  await PackageFunctions.consensusPharmacophoreDemo();
}
