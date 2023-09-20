import {test, category} from '@datagrok-libraries/utils/src/test';
import {app, pricesPanel, samplesPanel} from '../package';

category('Chemspace', () => {
  const mol = 'Oc1ccccc1';

  test('Prices panel', async () => {
    pricesPanel(mol);
  }, {skipReason: 'Requires API key'});

  test('Samples panel', async () => {
    samplesPanel(mol);
  }, {skipReason: 'Requires API key'});

  test('App', async () => {
    app();
  }, {skipReason: 'Requires API key'});
});
