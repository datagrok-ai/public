import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: problems', () => {
  test('a key no stack hash has is null', async () => {
    expect(await grok.dapi.admin.getProblem('error-incident:zzzzzzzzzzzz'), null);
  }, {owner: 'aparamonov@datagrok.ai'});

  test('an id that is neither a uuid nor kind:key is null', async () => {
    expect(await grok.dapi.admin.getProblem('nothing'), null);
  }, {owner: 'aparamonov@datagrok.ai'});
});
