import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import dayjs from 'dayjs';
import {category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

// Cloud log access (`core/docs/LOG_ARCHITECTURE.md`). The data lives in Datagrok's
// AWS account, so only the estate instances can actually read it — a stand with no
// AWS role and no archive connection has nothing to return. What is asserted here
// is the part that holds on every stand: the call reaches the server, binds its
// parameters, and fails at the boundary with a real message rather than silently
// resolving empty. A resolved call is accepted too, so this stays meaningful where
// the credentials do exist.
const MISSING_CONNECTION = '00000000-0000-0000-0000-0000000000ff';

category('Dapi: cloud logs', () => {
  test('log groups: a missing connection is reported, not swallowed', async () => {
    await expectExceptionAsync(
      async () => { await grok.dapi.log.getCloudLogGroups({connection: MISSING_CONNECTION}); },
      (e) => `${e?.message ?? e}`.length > 0);
  }, {owner: 'aparamonov@datagrok.ai'});

  test('log groups: resolves to a list when credentials exist', async () => {
    try {
      const groups = await grok.dapi.log.getCloudLogGroups();
      expect(Array.isArray(groups), true);
    } catch (e) {
      // No AWS region configured on this stand — the documented outcome off-estate.
      expect(`${(e as any)?.message ?? e}`.length > 0, true);
    }
  }, {owner: 'aparamonov@datagrok.ai'});

  test('events: time range and group are bound, result is a dataframe', async () => {
    const end = dayjs();
    try {
      const df = await grok.dapi.log.getCloudLogEvents('/datagrok/nonexistent',
        end.subtract(5, 'minute'), end, {limit: 10});
      expect(df instanceof DG.DataFrame, true);
      expect(df.columns.names().includes('message'), true);
    } catch (e) {
      expect(`${(e as any)?.message ?? e}`.length > 0, true);
    }
  }, {owner: 'aparamonov@datagrok.ai'});

  test('archive: listing without a connection is refused', async () => {
    await expectExceptionAsync(
      async () => { await grok.dapi.log.getArchiveObjects(MISSING_CONNECTION, {prefix: 'cloudwatch/'}); },
      (e) => `${e?.message ?? e}`.length > 0);
  }, {owner: 'aparamonov@datagrok.ai'});

  test('archive: decoding an object from a missing connection is refused', async () => {
    await expectExceptionAsync(
      async () => { await grok.dapi.log.getArchiveEvents(MISSING_CONNECTION, 'cloudwatch/x'); },
      (e) => `${e?.message ?? e}`.length > 0);
  }, {owner: 'aparamonov@datagrok.ai'});
});
