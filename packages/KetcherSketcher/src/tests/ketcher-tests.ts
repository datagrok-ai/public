import { after, before, category, test } from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import { _testSetMolfile, _testSetSmarts, _testSetSmiles } from './ketcher-utils';


category('ketcher', async () => {
  let previousSketcherType: string | undefined;

  before(async () => {
    previousSketcherType = grok.chem.currentSketcherType;
  });

  after(async () => {
    if (previousSketcherType !== undefined)
      grok.chem.currentSketcherType = previousSketcherType;
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

});
