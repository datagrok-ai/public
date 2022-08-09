import { after, before, category, delay, expect, expectFloat, test } from '@datagrok-libraries/utils/src/test';
import { Func } from 'datagrok-api/src/entities';
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