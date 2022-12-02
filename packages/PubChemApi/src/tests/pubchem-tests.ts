import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {inchiKeysToSmiles, pubChemToSmiles} from '../package';
import {smilesToPubChem} from '../pubchem';
import {buildAccordion, getSearchWidget} from '../widget';
import * as CONST from './const';

category('Panels', () => {
  const molStrings = [CONST.SMILES, CONST.MOL2000, CONST.MOL3000, CONST.SMARTS, CONST.EMPTY];

  test('Info', async () => {
    for (const molString of molStrings) {
      const pubChemId = await smilesToPubChem(molString);
      await buildAccordion(pubChemId);
    }
  });

  test('Substructure Search', async () => {
    for (const molString of molStrings)
      await getSearchWidget(molString, 'substructure');
  });

  test('Similarity Search', async () => {
    for (const molString of molStrings)
      await getSearchWidget(molString, 'similarity');
  });

  test('Identity Search', async () => {
    for (const molString of molStrings)
      await getSearchWidget(molString, 'identity');
  });

  test('pubChemToSmiles', async () => {
    expect(await pubChemToSmiles('3334'), 'COC(=O)NC1=NC2=C(N1)C=C(C=C2)SC3=CC=CC=C3');
    expect(await pubChemToSmiles('44219'), 'CC(=O)OC1=CC=CC=C1C(=O)[O-].C(CC[NH3+])CC(C(=O)O)N');
    expect(await pubChemToSmiles('2244'), 'CC(=O)OC1=CC=CC=C1C(=O)O');
    expect(await pubChemToSmiles('12762'), 'C1=CC=C(C=C1)[As](Cl)Cl');
  });

  test('inchiKeysToSmiles', async () => {
    expect(await inchiKeysToSmiles('UDHDFEGCOJAVRE-UHFFFAOYSA-N'), 'C1=CC=C(C=C1)[As](Cl)Cl');
    expect(await inchiKeysToSmiles('QWVGKYWNOKOFNN-UHFFFAOYSA-N'), 'CC1=CC=CC=C1O');
    expect(await inchiKeysToSmiles('FEWJPZIEWOKRBE-UHFFFAOYSA-N'), 'C(C(C(=O)O)O)(C(=O)O)O');
    expect(await inchiKeysToSmiles('VFUGCQKESINERB-UHFFFAOYSA-N'), 'CCCC1(CCN(C1)C)C2=CC(=CC=C2)O');
  });
});
