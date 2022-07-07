import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import {mmSemType} from '../const';
// import {importFasta} from '../package';

category('converters', () => {
//   test('a', async () => {await _a();});
//   test('b', async () => {await _b();});
  test('testFastaToSeparator', async () => { await _testFastaToSeparator(); });
  test('testSeparatorToFasta', async () => { await _testSeparatorToFasta(); });
});

// export async function _a() {
//   expect(1, 1);
// }
// 
// export async function _b() {
//   expect(1, 2);
// }

export async function _testFastaToSeparator() {
  expect(1, 1);
}

export async function _testSeparatorToFasta() {
  expect(1, 2);
}
