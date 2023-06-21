import {test, category} from '@datagrok-libraries/utils/src/test';
import {enamineStoreApp, enamineStorePanel} from '../package';

category('Enamine Store', () => {
  const mol = 'Oc1ccccc1';

  test('Panel', async () => {
    enamineStorePanel(mol);
  });

  test('App', async () => {
    enamineStoreApp();
  });
});
