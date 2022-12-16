import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {searchWidget, drugNameMoleculeConvert} from '../widgets';
import * as CONST from './const';

category('DrugBank', () => {
  const molStrings = [CONST.SMILES, CONST.SMARTS, CONST.MOL2000, CONST.MOL3000, CONST.EMPTY];
  let dbdf: DG.DataFrame;
  let synonymsCol: DG.Column<string>;
  let moleculeCol: DG.Column<string>;
  let dbdfRowCount: number;

  before(async () => {
    dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
    synonymsCol = dbdf.getCol('SYNONYMS');
    moleculeCol = dbdf.getCol('molecule');
    dbdfRowCount = dbdf.rowCount;
  });

  test('similarity-search', async () => {
    for (const molString of molStrings)
      await searchWidget(molString, 'similarity', dbdf);
  });

  test('substructure-search', async () => {
    for (const molString of molStrings)
      await searchWidget(molString, 'substructure', dbdf);
  });

  test('drugNameMolecule', async () => {
    drugNameMoleculeConvert('db:aspirin', dbdfRowCount, synonymsCol, moleculeCol)
    // expect(drugNameMoleculeConvert('db:aspirin', dbdfRowCount, synonymsCol, smilesCol), 'CC(Oc(cccc1)c1C(O)=O)=O');
    // expect(drugNameMoleculeConvert('db:carbono', dbdfRowCount, synonymsCol, smilesCol), '[C]');
    // expect(drugNameMoleculeConvert('db:gadolinio', dbdfRowCount, synonymsCol, smilesCol), '[Gd]');
  });
});
