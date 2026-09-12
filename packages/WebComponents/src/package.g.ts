import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: Init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}
