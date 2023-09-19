// import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';

import {category, test, expect, after} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

category('Packages: Properties', () => {
  test('INT', async () => {
    await changeProp('INT', 10);
  }, {skipReason: 'skip'});

  test('BIG_INT', async () => {
    await changeProp('BIG_INT', 10);
  }, {skipReason: 'skip'});

  test('FLOAT', async () => {
    await changeProp('FLOAT', 10.567);
  }, {skipReason: 'skip'});

  test('NUM', async () => {
    await changeProp('NUM', 1.234);
  }, {skipReason: 'skip'});

  test('BOOL', async () => {
    await changeProp('BOOL', false);
  }, {skipReason: 'skip'});

  test('STRING', async () => {
    await changeProp('STRING', 'value2');
  }, {skipReason: 'skip'});

  after(async () => {
    await _package.setProperties({INT: 1, BIG_INT: 1, FLOAT: 1.234, NUM: 1, BOOL: true, STRING: 'value1'});
  });
});

async function changeProp(name: string, value: any): Promise<void> {
  await _package.setProperties({[name]: value});
  const props = await _package.getProperties();
  expect(props.get(name), value);
}
