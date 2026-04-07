import {category, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {searchWidget, drugNameMoleculeConvert, SEARCH_TYPE} from '../widgets';
import * as CONST from './const';

//NOTE: Skip reason here requires manual check of the test results. The tests are disabled for GH Actions.

category('DrugBank', () => {
  const molStrings = [CONST.SMILES, CONST.SMARTS, CONST.MOL2000, CONST.MOL3000, CONST.EMPTY];

  test('similarity-search', async () => {
    const dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
    for (const molString of molStrings)
      await searchWidget(molString, SEARCH_TYPE.SIMILARITY, dbdf);
  });

  test('substructure-search', async () => {
    const dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
    for (const molString of molStrings)
      await searchWidget(molString, SEARCH_TYPE.SUBSTRUCTURE, dbdf);
  });

  test('drugNameMolecule', async () => {
    const dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
    drugNameMoleculeConvert('db:aspirin', dbdf.rowCount, dbdf.getCol('SYNONYMS'), dbdf.getCol('molecule'));
    // expect(drugNameMoleculeConvert('db:aspirin', dbdfRowCount, synonymsCol, smilesCol), 'CC(Oc(cccc1)c1C(O)=O)=O');
    // expect(drugNameMoleculeConvert('db:carbono', dbdfRowCount, synonymsCol, smilesCol), '[C]');
    // expect(drugNameMoleculeConvert('db:gadolinio', dbdfRowCount, synonymsCol, smilesCol), '[Gd]');
  });
});
