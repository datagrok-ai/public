import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, test, delay, before, after} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure, _testSearchSubstructureAllParameters} from './utils';
import {_testFindSimilar, _testGetSimilarities} from './menu-tests-similarity-diversity';
import {testCsv, testSubstructure} from './substructure-search-tests';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import {_importSdf} from '../open-chem/sdf-importer';

category('server features', () => {
  test('descriptors', async () => {

    const tree = await grok.chem.descriptorsTree();
    expect(tree !== undefined, true);
    const df = DG.DataFrame.fromCsv(testCsv);
    const t: DG.DataFrame = await grok.chem.descriptors(df, 'smiles', 
      ['MolWt', 'NumAromaticCarbocycles','NumHAcceptors', 'NumHeteroatoms', 'NumRotatableBonds', 'RingCount']);
    grok.shell.addTableView(t);

    isColumnPresent(grok.shell.t.columns, 'MolWt');
    isColumnPresent(grok.shell.t.columns, 'NumAromaticCarbocycles');
    isColumnPresent(grok.shell.t.columns, 'NumHAcceptors');
    isColumnPresent(grok.shell.t.columns, 'NumHeteroatoms');
    isColumnPresent(grok.shell.t.columns, 'NumRotatableBonds');
    isColumnPresent(grok.shell.t.columns, 'RingCount');
  });

  test('sketcher', async () => {
    const result: HTMLElement = grok.chem.sketcher(()=>{}, 'CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3');
    expect(result !== null, true);
  });
});



category('chem exported', () => {
  test('findSimilar.api.sar-small', async () => {
    await _testFindSimilar(grok.chem.findSimilar);
  });

  test('getSimilarities.api.molecules', async () => {
    await _testGetSimilarities(grok.chem.getSimilarities);
  });

  test('substructureSearch', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const trueIndices = [0, 2];
    const bitset: DG.BitSet = await grok.chem.searchSubstructure(df.col('smiles')!, testSubstructure);
    const bitsetString = bitset.toBinaryString();
    const bitsetArray = [...bitsetString];
    for (let k = 0; k < trueIndices.length; k++) {
      expect(bitsetArray[trueIndices[k]] === '1', true);
      bitsetArray[trueIndices[k]] = '0';
    }
  });

  test('mcs', async () => {
    // const df = DG.DataFrame.fromCsv(testCsv);
    // const mcs = await grok.chem.mcs(df, 'smiles');
    // expect(mcs, 'C:CCC1:C:C:C:C:C:1');
  });
});
