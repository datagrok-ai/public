import * as DG from 'datagrok-api/dg';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {aligned1} from './test-data';


category('Viewers: Basic', () => {
  const df = DG.DataFrame.fromCsv(aligned1);
  const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), true);
    }, {skipReason: 'GROK-11534'});
  }
});

category('Viewers: Monomer-Position', () => {
  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Modes', async () => {

  }, {skipReason: 'Not implemented yet'});
});

category('Viewers: Most Potent Residues', () => {
  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});
});

category('Viewers: Logo Summary Table', () => {
  test('Properties', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});
});
