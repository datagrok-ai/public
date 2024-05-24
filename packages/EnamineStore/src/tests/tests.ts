import * as grok from 'datagrok-api/grok';
import {test, category, expect} from '@datagrok-libraries/utils/src/test';
import {EnamineStoreSearchResult, SEARCH_MODE, enamineStoreApp, enamineStorePanel, searchModeToCommandMap} from '../package';

category('Enamine Store', () => {
  const mol = 'Oc1ccccc1';

  test('Panel', async () => {
    enamineStorePanel(mol);
  });
 
  test('App', async () => {
    enamineStoreApp();
  }); 

  test('Similarity Search', async () => {
    await testEnamineSearch(SEARCH_MODE.SIMILAR, 5);
  });

  test('Exact Search', async () => {
    await testEnamineSearch(SEARCH_MODE.EXACT, 1);
  });

  test('Substructure Search', async () => {
    await testEnamineSearch(SEARCH_MODE.SUBSTRUCTURE, 9);
  });
});

async function testEnamineSearch(searchMode: SEARCH_MODE, numResults: number) {
  const searchType = searchModeToCommandMap[searchMode];
  const options: {[key: string]: any} = {
    'q': 'Oc1ccccc1',
    'cat': 'SCR',
    'currency': 'USD',
    'type': 'SMARTS',
    'sstype': searchType
  };
  const fc = await grok.data.callQuery('EnamineStore:Search', options, true, 100);
  const res = JSON.parse(fc.getParamValue('stringResult'))['searchResults'] as EnamineStoreSearchResult[];
  //console.log(JSON.stringify(res));
  expect(res.length, numResults, `Incorrect ${searchType} search results, returned ${numResults}`);

   const opts: {[key: string]: any} = {
    'id': res[0]['code'],
    'cat': 'SCR',
    'cur': 'USD',
  };
  const fcPrice = await grok.data.callQuery('EnamineStore:Price', opts, true, 100);
  const resPrice = JSON.parse(fcPrice.getParamValue('stringResult'))['samples'];
  expect(resPrice.length > 0, true, 'Incorrect price results'); 
}

