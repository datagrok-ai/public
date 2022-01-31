import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {drugBankSearchWidget} from '../widgets';

category('DrugBank', () => {
  const molString = 'C';
  let dbdf: DG.DataFrame;

  before(async () => {
    dbdf = DG.DataFrame.fromCsv(await _package.files.readAsText('db.csv'));
  });

  test('similarity-search', async () => {
    await drugBankSearchWidget(molString, 'similarity', dbdf);
  });

  test('substructure-search', async () => {
    await drugBankSearchWidget(molString, 'substructure', dbdf);
  });
});
