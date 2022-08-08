import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
//import {findMonomers, helmToFasta, helmToPeptide, helmToRNA, initHelm} from '../package';
import {_package} from '../package-test';
import {loadDialog, manageFiles, monomerManager} from '../package';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


category('Helm', () => {
  //These tests require webservice that is not present on test stand
  // test('helmToFasta', async () => {
  //   expect(await helmToFasta('RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$'), '>RNA1UTGCA');
  //   expect(await helmToFasta('RNA1{P.R(U).P.R(T)}$$$$'), '>RNA1UT');
  //   expect(await helmToFasta('PEPTIDE1{A.G}$$$$V2.0'), '>PEPTIDE1AG');
  // });
  //
  // test('helmToRNA', async () => {
  //   expect(await helmToRNA('RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$'), 'UTGCA');
  //   expect(await helmToRNA('RNA1{P.R(U).P.R(T)}$$$$'), 'UT');
  //   // eslint-disable-next-line max-len
  //   expect(await helmToRNA('RNA1{R(U)P.R(T)P}|RNA2{P.R(A)P.R(A)}$RNA1,RNA2,2:pair-6:pair|RNA1,RNA2,5:pair-3:pair$$$'), 'UT AA');
  // });
  //
  // test('helmToPeptide', async () => {
  //   expect(await helmToPeptide('PEPTIDE1{A.G}$$$$V2.0'), 'AG');
  //   expect(await helmToPeptide('PEPTIDE1{L.V.A}|PEPTIDE2{L.V.A}$$$$'), 'LVA LVA');
  //   // eslint-disable-next-line max-len
  //   expect(await helmToPeptide('PEPTIDE1{A.R.C.A.A.K.T.C.D.A}$PEPTIDE1,PEPTIDE1,8:R3-3:R3$$$'), 'ARCAAKTCDA');
  // });

  // detectMacromolecule is a function of Bio package
  // test('detectMacromolecule', async () => {
  //   const file = await _package.files.readAsText('tests/test.csv');
  //   const df = DG.DataFrame.fromCsv(file);
  //   const col = df.columns.byName('HELM string');
  //   await grok.data.detectSemanticTypes(df);
  //   expect(col.semType, DG.SEMTYPE.MACROMOLECULE);
  // });

  // semType and units tags are detecting by Bio package function detectMacromolecule
  // test('detectHelm', async () => {
  //   const file = await _package.files.readAsText('tests/test.csv');
  //   const df = DG.DataFrame.fromCsv(file);
  //   const col = df.columns.byName('HELM string');
  //   await grok.data.detectSemanticTypes(df);
  //   expect(col.tags[DG.TAGS.UNITS], 'HELM');
  // });


  test('manageFiles', async () => {
    return await manageFiles();
  });

  test('loadDialog', async() => {
    return await loadDialog();
  });

  test('monomerManager', async() => {
    await grok.functions.call('Helm:monomerManager', {value: 'HELMCoreLibrary.json'});
    const checkName = '2-Chloroadenine';
    let flag = false;
    const types = Object.keys(org.helm.webeditor.monomerTypeList());
    const monomers: any = [];
    for (var i = 0; i < types.length; i++) {
      //@ts-ignore
      monomers.push(new scil.helm.Monomers.getMonomerSet(types[i]));
      Object.keys(monomers[i]).forEach(k => {
        if (monomers[i][k].n == checkName){
          flag = true;
        }
      });
    }
    expect(flag, true);
  });

});
