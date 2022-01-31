import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {drugBankSimilaritySearchPanel, drugBankSubstructureSearchPanel, initDrugBank} from '../package';

category('DrugBank', () => {
  const molString = 'C';
  before(async () => {
    await initDrugBank();
  });

  test('similarity-search', async () => {
    await drugBankSimilaritySearchPanel(molString);
  });

  test('substructure-search', async () => {
    await drugBankSubstructureSearchPanel(molString);
  });
});
