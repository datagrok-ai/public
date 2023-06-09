import {test, category} from '@datagrok-libraries/utils/src/test';
import {app, pricesPanel, samplesPanel} from '../package';

category('Chemspace', () => {
  const mol = 'Oc1ccccc1';

  test('Prices panel', async () => {
    pricesPanel(mol);
  });

  test('Samples panel', async () => {
    samplesPanel(mol);
  });

  test('App', async () => {
    app();
  });
});
