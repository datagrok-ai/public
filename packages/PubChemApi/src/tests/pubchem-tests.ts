import {category, test} from '@datagrok-libraries/utils/src/test';
import {pubChemSearchWidget} from '../widget';

category('PubChem API', () => {
  const molString = 'C';

  test('Substructure Search', async () => {
    await pubChemSearchWidget(molString, 'substructure');
  });

  test('Similarity Search', async () => {
    await pubChemSearchWidget(molString, 'similarity');
  });

  test('Identity Search', async () => {
    await pubChemSearchWidget(molString, 'identity');
  });
});
