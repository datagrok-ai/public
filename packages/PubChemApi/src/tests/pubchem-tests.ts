import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {inchiKeys, pubChem} from '../package';
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
  
  test('pubChem', async () => {
    expect(await pubChem('3334'), 'COC(=O)NC1=NC2=C(N1)C=C(C=C2)SC3=CC=CC=C3');
    expect(await pubChem('44219'), 'CC(=O)OC1=CC=CC=C1C(=O)[O-].C(CC[NH3+])CC(C(=O)O)N');
    expect(await pubChem('2244'), 'CC(=O)OC1=CC=CC=C1C(=O)O');
    expect(await pubChem('12762'), 'C1=CC=C(C=C1)[As](Cl)Cl');
  });

  test('inchiKeys', async () => {
    expect(await inchiKeys('UDHDFEGCOJAVRE-UHFFFAOYSA-N'), 'C1=CC=C(C=C1)[As](Cl)Cl');
    expect(await inchiKeys('QWVGKYWNOKOFNN-UHFFFAOYSA-N'), 'CC1=CC=CC=C1O');
    expect(await inchiKeys('FEWJPZIEWOKRBE-UHFFFAOYSA-N'), 'C(C(C(=O)O)O)(C(=O)O)O');
    expect(await inchiKeys('VFUGCQKESINERB-UHFFFAOYSA-N'), 'CCCC1(CCN(C1)C)C2=CC(=CC=C2)O');
  });
});
