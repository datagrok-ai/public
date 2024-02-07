import { after, before, category, test } from '@datagrok-libraries/utils/src/test';
import { _testSetMolfile, _testSetSmarts, _testSetSmiles } from './ketcher-utils';


category('ketcher', async () => {

  before(async () => {
  });

  test('setSmiles', async () => {
    await _testSetSmiles();
  });

  test('setMolfile', async () => {
    await _testSetMolfile();
  });

  test('setSmarts', async () => {
    await _testSetSmarts();
  });

  after(async () => {
  });

});