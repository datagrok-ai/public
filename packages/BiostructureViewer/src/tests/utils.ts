import {expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PdbEntry} from '../pdb-entry';

export async function requireText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

/**
 * Tests if RCSB REST API returns any live data for a sample PDB id.
 *
 * @param {string} [pdbID='2V0A']
 */
export async function _testRCSBAlive(pdbID = '2V0A') {
  const pdb = new PdbEntry(pdbID);
  let noException = true;

  try {
    await pdb.fetchInfo();
  } catch (error) {
    noException = false;
  }
  expect(noException, true);
  expect(pdb.body.length > 0, true);
  expect(pdb.entities.length > 0, true);
}
