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
  expect(pdb.entities[0].chains.length > 0, true);
}

/** Test if parsing of a sample PDB id is done properly. */
export async function _testParse2V0A() {
  const pdb = new PdbEntry('2V0A');

  await pdb.fetchInfo();

  expect(pdb.entities.length == 1, true);
  expect(pdb.entities[0].chains.length == 2, true);

  for (const chain of pdb.entities[0].chains) {
    const keys = Object.keys(chain.tracks);
    expect(keys.includes('SHEET'), true);
    expect(keys.includes('HELIX_P'), true);
    expect(chain.tracks['SHEET'].length > 0, true);
    expect(chain.tracks['HELIX_P'].length > 0, true);
  }
}
