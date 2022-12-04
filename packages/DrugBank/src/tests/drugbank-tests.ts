import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {searchWidget, drugNameMoleculeConvert} from '../widgets';

category('DrugBank', () => {
  const molString = 'C';
  let dbdf: DG.DataFrame;
  let synonymsCol: DG.Column<string>;
  let smilesCol: DG.Column<string>;
  let dbdfRowCount: number;

  before(async () => {
    dbdf = DG.DataFrame.fromCsv(await _package.files.readAsText('db.csv'));
    synonymsCol = dbdf.getCol('SYNONYMS');
    smilesCol = dbdf.getCol('Smiles');
    dbdfRowCount = dbdf.rowCount;
  });

  test('similarity-search', async () => {
    await searchWidget(molString, 'similarity', dbdf);
  });

  test('substructure-search', async () => {
    await searchWidget(molString, 'substructure', dbdf);
  });

  test('drugNameMolecule', async () => {
    expect(drugNameMoleculeConvert('db:aspirin', dbdfRowCount, synonymsCol, smilesCol), 'CC(Oc(cccc1)c1C(O)=O)=O');
    expect(drugNameMoleculeConvert('db:carbono', dbdfRowCount, synonymsCol, smilesCol), '[C]');
    expect(drugNameMoleculeConvert('db:gadolinio', dbdfRowCount, synonymsCol, smilesCol), '[Gd]');
  });
});
