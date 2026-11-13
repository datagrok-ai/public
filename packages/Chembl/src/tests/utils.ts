import * as grok from 'datagrok-api/grok';

import {skipCategoryWithoutSpawner} from '@datagrok-libraries/test/src/test-container-utils';

/** Skips a category's tests when the ChEMBL connection points at the package's container and the stand cannot start
 * it: every query would hang until its test's timeout. */
export async function skipWithoutChemblDatabase(categoryName: string): Promise<void> {
  const conn = (await grok.dapi.connections.filter('namespace = "Chembl:" and shortName = "Chembl"').list())[0];
  if (`${conn?.parameters?.server}`.includes('DockerContainer'))
    await skipCategoryWithoutSpawner(categoryName, 'The ChEMBL database');
}
