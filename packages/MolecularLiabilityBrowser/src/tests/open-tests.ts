import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import {packageName} from '../package';
import {DataQuery} from 'datagrok-api/dg';
import {catchToLog, DataLoader} from '../utils/data-loader';

category('open', () => {
  test('open1', async () => { await testOpen1(); });

  test('loadData1', async () => { await testLoadData1(); });

  async function testOpen1(): Promise<void> {
    expect(1, 1);
  }

  async function testLoadData1(): Promise<void> {


    let k = 11;

    // const queryNameToQuery = {};
    // queries.forEach(e => queryNameToQuery[e.name.charAt(0).toLowerCase() + e.name.slice(1)] = e);
    //
    // // query
    // const call = queryNameToQuery['firstQuery'].prepare();
    // await call.call();
  }
});

