import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

async function compare(input: string, expectArr: Array<any>, callMethod: string, columnName: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(input);

  let dfResult = null;
  try {
    dfResult = await grok.functions.call(callMethod, {'ids': df});
  }
  catch (e) {
    console.log(e);
  }

  expect(dfResult!.rowCount, expectArr.length);
  for (let i = 0; i < dfResult!.rowCount; i++) {
    expect(dfResult!.get(columnName, i), expectArr[i]);
  }
}

category('Converting queries', () => {

  test('converter Id to Smiles', async () => {
    const smiles = await grok.functions.call(`${_package.name}:chemblIdToSmiles`, {'id': 'CHEMBL6302'});
    expect(smiles, 'COc1cc2nc(N3CCN(/C(S)=N/c4ccc(NC(=S)N5CC5)cc4)CC3)nc(N)c2cc1OC');

    const smiles2 = await grok.functions.call(`${_package.name}:chemblIdToSmiles`, {'id': 'CHEMBL6323'});
    expect(smiles2, 'COc1cc2nc(N3CCN(C(=O)c4ccc5ccccc5c4)CC3)nc(N)c2cc1OC');

    const smiles3 = await grok.functions.call(`${_package.name}:chemblIdToSmiles`, {'id': 'CHEMBL6365'});
    expect(smiles3, 'O=c1oc2ccccc2c2cc(O)c(O)cc12');
  });

  test('converter Molregno to Smiles', async () => {
    const smiles = await grok.functions.call(`${_package.name}:molregnoToSmiles`, {'molregno': '1824756'});
    expect(smiles, 'CCCCCCCCCOc1ccc(F)c(C(=O)NCCN)c1F');

    const smiles2 = await grok.functions.call(`${_package.name}:molregnoToSmiles`, {'molregno': '1824757'});
    expect(smiles2, 'CCCCCCCCCOc1ccc(F)c(C(=O)NCCNC(=O)c2c(F)ccc(OCCCCCCCCC)c2F)c1F');

    const smiles3 = await grok.functions.call(`${_package.name}:molregnoToSmiles`, {'molregno': '1824758'});
    expect(smiles3, 'CCCCCCCCCOc1ccc2c(c1F)C(=O)NCCN2');
  });

  test('converter Name to Smiles', async () => {
    const smiles = await grok.functions.call(`${_package.name}:nameToSmiles`,
      {'compoundName': 'Bis-quinolinium cyclophane analogue'});
    expect(smiles, 'c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1');

    const smiles2 = await grok.functions.call(`${_package.name}:nameToSmiles`,
      {'compoundName': '1-[4-(4-Amino-6,7-dimethoxy-quinazolin-2-yl)-piperazin-1-yl]-3-phenyl-propan-1-one'});
    expect(smiles2, 'COc1cc2nc(N3CCN(C(=O)CCc4ccccc4)CC3)nc(N)c2cc1OC');

    const smiles3 = await grok.functions.call(`${_package.name}:nameToSmiles`,
      {'compoundName': '2-{4-[(4-Chloro-phenyl)-hydroxy-methyl]-3,5-dimethyl-phenyl}-2H-[1,2,4]triazine-3,5-dione'});
    expect(smiles3, 'Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1');
  });

  test('converter InchiKey to Chembl', async () => {
    let csvDf1 = `key,order
    OWRSAHYFSSNENM-UHFFFAOYSA-N,1
    ZJYUMURGSZQFMH-UHFFFAOYSA-N,4
    YOMWDCALSDWFSV-UHFFFAOYSA-N,2
    PSOPUAQFGCRDIP-UHFFFAOYSA-N,5`;
    let expectArr = ['CHEMBL6329', 'CHEMBL265667', 'CHEMBL6328', 'CHEMBL6362'];

    await compare(csvDf1, expectArr, `${_package.name}:inchiKeyToChembl`, 'chembl_id');
  });

  test('converter InchiKey to Smiles', async () => {
    let csvDf1 = `key,order
    OWRSAHYFSSNENM-UHFFFAOYSA-N,1
    ZJYUMURGSZQFMH-UHFFFAOYSA-N,4
    YOMWDCALSDWFSV-UHFFFAOYSA-N,2
    PSOPUAQFGCRDIP-UHFFFAOYSA-N,5`;
    let expectArr = ['Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl', 'Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1',
      'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1', 'Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1'];

    await compare(csvDf1, expectArr, `${_package.name}:inchiKeyToSmiles`, 'canonical_smiles');
  });
  
  test('converter InchiKey to Inchi', async () => {
    let csvDf1 = `key,order
    OWRSAHYFSSNENM-UHFFFAOYSA-N,1
    ZJYUMURGSZQFMH-UHFFFAOYSA-N,4
    YOMWDCALSDWFSV-UHFFFAOYSA-N,2
    PSOPUAQFGCRDIP-UHFFFAOYSA-N,5`;
    let expectArr = ['InChI=1S/C17H12ClN3O3/c1-10-8-11(21-17(24)20-15(22)9-19-21)6-7-12(10)16(23)13-4-2-3-5-14(13)18/h2-9H,1H3,(H,20,22,24)',
      'InChI=1S/C18H16ClN3O3/c1-10-7-14(22-18(25)21-15(23)9-20-22)8-11(2)16(10)17(24)12-3-5-13(19)6-4-12/h3-9,17,24H,1-2H3,(H,21,23,25)',
      'InChI=1S/C18H12N4O3/c1-11-8-14(22-18(25)21-16(23)10-20-22)6-7-15(11)17(24)13-4-2-12(9-19)3-5-13/h2-8,10H,1H3,(H,21,23,25)',
      'InChI=1S/C17H13N3O3/c1-11-2-4-12(5-3-11)16(22)13-6-8-14(9-7-13)20-17(23)19-15(21)10-18-20/h2-10H,1H3,(H,19,21,23)'];

    await compare(csvDf1, expectArr, `${_package.name}:inchiKeyToInchi`, 'standard_inchi');
  });

  test('converter Chembl to Smiles', async () => {
    let csvDf1 = `key, order
                  CHEMBL6329,1
                  CHEMBL267864,5
                  CHEMBL6328,2
                  CHEMBL6362,4
                  CHEMBL265667,3`;

    let expectArr = ['Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl',
      'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(C#N)cc1',
      'Cc1cc(-n2ncc(=O)[nH]c2=O)cc(C)c1C(O)c1ccc(Cl)cc1',
      'Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1',
      'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccc(Cl)cc1'];

    await compare(csvDf1, expectArr, `${_package.name}:chemblToSmiles`, 'canonical_smiles');
  });

  test('converter Chembl to Inchi', async () => {
    let csvDf1 = `key, order
                  CHEMBL6329,1
                  CHEMBL267864,5
                  CHEMBL6328,2
                  CHEMBL6362,4
                  CHEMBL265667,3`;
    let expectArr = ['InChI=1S/C17H12ClN3O3/c1-10-8-11(21-17(24)20-15(22)9-19-21)6-7-12(10)16(23)13-4-2-3-5-14(13)18/h2-9H,1H3,(H,20,22,24)',
      'InChI=1S/C18H12N4O3/c1-11-8-14(22-18(25)21-16(23)10-20-22)6-7-15(11)17(24)13-4-2-12(9-19)3-5-13/h2-8,10H,1H3,(H,21,23,25)',
      'InChI=1S/C18H16ClN3O3/c1-10-7-14(22-18(25)21-15(23)9-20-22)8-11(2)16(10)17(24)12-3-5-13(19)6-4-12/h3-9,17,24H,1-2H3,(H,21,23,25)',
      'InChI=1S/C17H13N3O3/c1-11-2-4-12(5-3-11)16(22)13-6-8-14(9-7-13)20-17(23)19-15(21)10-18-20/h2-10H,1H3,(H,19,21,23)',
      'InChI=1S/C17H12ClN3O3/c1-10-8-13(21-17(24)20-15(22)9-19-21)6-7-14(10)16(23)11-2-4-12(18)5-3-11/h2-9H,1H3,(H,20,22,24)'];

    await compare(csvDf1, expectArr, `${_package.name}:chemblToInchi`, 'standard_inchi');
  });
  test('converter Chembl to InchiKey', async () => {
    let csvDf1 = `key, order
                  CHEMBL6329,1
                  CHEMBL267864,5
                  CHEMBL6328,2
                  CHEMBL6362,4
                  CHEMBL265667,3`;
    let expectArr = ['OWRSAHYFSSNENM-UHFFFAOYSA-N',
      'ZJYUMURGSZQFMH-UHFFFAOYSA-N',
      'YOMWDCALSDWFSV-UHFFFAOYSA-N',
      'PSOPUAQFGCRDIP-UHFFFAOYSA-N',
      'KEZNSCMBVRNOHO-UHFFFAOYSA-N'];


    await compare(csvDf1, expectArr, `${_package.name}:chemblToInchiKey`, 'standard_inchi_key');
  });

});
