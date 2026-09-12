import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {DomainFrameEditor} from '@datagrok-libraries/domain-ui';
import {thrown, withRestrictedUser} from './domain-lifecycle';

// The access surface: DomainTableClient.access() — ONE server-composed answer in the
// shape every permission-aware control consumes ({can, fields}) — and the per-row
// `withAccess` columns that replaced the per-row permissions endpoint. Real server,
// real grants; the restricted-user probes ride withRestrictedUser (the user is
// blocked in finally) and skip cleanly where self-signup is unavailable.
category('Dapi: domain access', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const hidden = () => grok.dapi.domains.table('apitests.hidden_item');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const ACCESS_KEYS = [...DG.DOMAIN_ACCESS_COLUMNS];
  // `~can_share` is per row only in row mode (apitests.item is): a boolean there, null elsewhere.
  const isAccessValue = (k: string, v: any) => typeof v === 'boolean' || (k === '~can_share' && v === null);

  test('shape: can, fields, identity; cache survives reads, drops on invalidation', async () => {
    const access = await items().access();
    const can = access.can;
    expect(can.view && can.insert && can.edit && can.delete && can.share, true, JSON.stringify(can));
    expect(access.securityMode, 'row');
    expect(access.securingTable, 'apitests.item', 'a row-mode table secures itself');
    expect(access.audit, true);
    expect(access.hasBusinessKey, true);
    for (const c of ['sku', 'name', 'quantity'])
      expect(access.fields[c], 'editable', `${c} must be editable for admin: ${JSON.stringify(access.fields)}`);
    for (const c of ['id', 'version'])
      expect(c in access.fields, true, `system column ${c} must be listed: ${JSON.stringify(access.fields)}`);
    expect(Object.keys(access.fields).some((c) => c.startsWith('~')), false,
      `service columns are not fields: ${JSON.stringify(access.fields)}`);
    expect(access.travelableRelations.includes('tags'), true,
      `the declared tags relation must be travelable for admin: ${JSON.stringify(access.travelableRelations)}`);
    grok.dapi.domains.invalidateUiCaches();
    expect((await items().access()).can.edit, true, 'recompute after invalidation must succeed');
  });

  test('master mode secures on the delegate target', async () => {
    const junction = await grok.dapi.domains.table('apitests.item_tag').access();
    expect(junction.securityMode, 'master');
    expect(junction.securingTable, 'apitests.item', JSON.stringify(junction));
    expect(junction.travelableRelations.length, 0, 'the junction declares no relations');
  });

  test('unknown table rejects with a typed validation error', async () => {
    const e = await thrown(() => grok.dapi.domains.table('apitests.nosuch').access());
    expect(e instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
  });

  test('withAccess: JSON rows and get() carry the boolean row keys, plain reads do not', async () => {
    const sku = `SKU-ACC-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Access probe'});
    try {
      const filter = `sku = "${sku}"`;
      const [plain] = await items().query({filter});
      for (const k of ACCESS_KEYS)
        expect(k in plain, false, `${k} must be absent without the flag: ${JSON.stringify(plain)}`);
      const [row] = await items().query({filter, withAccess: true});
      for (const k of ACCESS_KEYS)
        expect(isAccessValue(k, row[k]), true, `${k} must be an access flag: ${JSON.stringify(row)}`);
      expect(typeof row['~can_share'], 'boolean', `a row-mode table carries per-row Share: ${JSON.stringify(row)}`);
      expect(row['~can_edit'] && row['~can_delete'], true,
        `admin must edit and delete its own row: ${JSON.stringify(row)}`);
      const one = await items().get(ins.id, {withAccess: true});
      expect(one.id, ins.id, 'get(withAccess) must resolve the addressed row');
      for (const k of ACCESS_KEYS)
        expect(isAccessValue(k, one[k]), true, `get(withAccess) must carry ${k}: ${JSON.stringify(one)}`);
      const bare = await items().get(ins.id);
      for (const k of ACCESS_KEYS)
        expect(k in bare, false, `${k} must be absent from a plain get: ${JSON.stringify(bare)}`);
      // Off row mode (item_event is master-secured) the per-row Share is null, not false.
      const events = grok.dapi.domains.table('apitests.item_event');
      const [ev] = await events.insert({item_id: ins.id, kind: `acc-${stamp()}`, amount: 1});
      const [event] = await events.query({filter: `id = "${ev.id}"`, withAccess: true});
      expect(event['~can_share'], null, `~can_share must be null off row mode: ${JSON.stringify(event)}`);
      expect(typeof event['~can_edit'], 'boolean', `~can_edit must still be a boolean: ${JSON.stringify(event)}`);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('withAccess: queryDf carries the row keys as bool columns', async () => {
    const sku = `SKU-ACC-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Access probe'});
    try {
      const filter = `sku = "${sku}"`;
      const df = await items().queryDf({filter, withAccess: true});
      expect(df.rowCount, 1);
      for (const k of ACCESS_KEYS) {
        const col = df.col(k);
        expect(col != null, true, `${k} column missing: ${df.columns.names().join(', ')}`);
        expect(col!.type, DG.TYPE.BOOL, `${k} must be a bool column, got ${col!.type}`);
        expect(col!.meta.includeInBinaryExport, false, `${k} is not excluded from binary export`);
        expect(col!.meta.includeInCsvExport, false, `${k} is not excluded from csv export`);
      }
      expect(df.toCsv().includes('~can_'), false, 'the access columns leaked into toCsv()');
      const bare = await items().queryDf({filter});
      for (const k of ACCESS_KEYS)
        expect(bare.col(k), null, `${k} must be absent without the flag`);
      // Off row mode a bool column has no null slot: `~can_share` is ABSENT from the
      // d42 frame (null in the JSON row), the other two stay.
      const events = grok.dapi.domains.table('apitests.item_event');
      const [ev] = await events.insert({item_id: ins.id, kind: `acc-${stamp()}`, amount: 1});
      const evDf = await events.queryDf({filter: `id = "${ev.id}"`, withAccess: true});
      expect(evDf.rowCount, 1);
      expect(evDf.col('~can_share'), null, `~can_share must be absent off row mode: ${evDf.columns.names().join(', ')}`);
      for (const k of ['~can_edit', '~can_delete'])
        expect(evDf.col(k)?.type, DG.TYPE.BOOL, `${k} must stay a bool column off row mode`);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('withAccess frame: the editor adds and saves a row, the insert carries no ~ keys', async () => {
    const prefix = `SKU-ACCFE-${stamp()}`;
    await items().insert({sku: `${prefix}-0`, name: 'Access frame seed'});
    try {
      const query = {filter: {property: 'sku', operator: 'like', value: `${prefix}%`} as any, withAccess: true};
      const df = await items().queryDf(query);
      expect(df.columns.contains('~can_edit'), true, 'the fixture frame carries no access columns');
      const editor = await DomainFrameEditor.attach(df, items() as any, {query, quiet: true});
      try {
        expect(editor.quiet, true, 'the quiet option did not reach the editor');
        const row = editor.addRow({sku: `${prefix}-1`, name: 'Added on a withAccess frame', quantity: 1});
        expect(editor.errorCount, 0,
          `an added row on a withAccess frame must validate: ${JSON.stringify(editor.errorsOf(row))}`);
        expect(editor.isChanged(row, '~can_edit'), false, 'an access column reads as a changed cell');
        const ops = editor.buildOps();
        expect(ops.length, 1, `one insert expected: ${JSON.stringify(ops)}`);
        const keys = Object.keys((ops[0].op as any).values ?? {});
        expect(keys.some((k) => k.startsWith('~')), false, `the insert payload carries service keys: ${keys}`);
        const balloons = document.querySelectorAll('.d4-balloon').length;
        expect(await editor.save(), true, 'save() must succeed on a withAccess frame');
        // Earlier balloons may expire meanwhile — only a NEW one is the editor's.
        expect(document.querySelectorAll('.d4-balloon').length <= balloons, true,
          'a quiet editor raised its own Saved balloon');
        expect((await items().query({filter: `sku = "${prefix}-1"`})).length, 1, 'the added row was not inserted');
        // The post-save re-read lands the server's row: system columns it filled
        // in, and the per-row access — not the payload and not the frame defaults.
        const id = df.get('id', row);
        for (const c of ['version', 'created_on', 'author_id'])
          expect(df.get(c, row) != null, true, `${c} is blank in the frame after save`);
        const fresh = await items().get(id, {withAccess: true});
        expect(df.get('~can_edit', row), fresh['~can_edit'],
          'the inserted row must carry the server\'s per-row edit right after save');
      } finally {
        editor.detach();
      }
    } finally {
      await items().deleteWhere({property: 'sku', operator: 'like', value: `${prefix}%`} as any);
    }
  });

  test('a restricted user: flags flip on a real grant round-trip, a restricted column is absent', async () => {
    const outcome = await withRestrictedUser('wo4acc', async (probe) => {
      // access() under the restricted user's session — the probe is a per-call
      // server request, so it follows the swapped auth; the client-side cache is
      // dropped first because its user key tracks the SESSION user.
      const asUser = () => probe.asUser(async () => {
        grok.dapi.domains.invalidateUiCaches();
        return await items().access();
      });
      try {
        const before = await asUser();
        expect(before.can.insert, false, `no grant yet, can.insert must deny: ${JSON.stringify(before)}`);
        expect(before.can.edit, false, 'no grant yet, can.edit must deny');
        // Travel is gated like the reads it rides: a ROW-secured junction and target
        // pass for any authenticated caller and let the row predicate hide the links.
        expect(before.travelableRelations.includes('tags'), true,
          `row-secured tags must stay travelable for an ungranted user: ${JSON.stringify(before.travelableRelations)}`);
        expect(before.securingTable, 'apitests.item', 'securingTable is identity, not permission');
        await items().grant(probe.group, 'Edit'); // drops the caches automatically
        const after = await asUser();
        expect(after.can.insert, true, `Edit grant must flip can.insert: ${JSON.stringify(after)}`);
        expect(after.can.edit, true, 'Edit grant must flip can.edit');
        await items().restrictColumn('quantity');
        try {
          const restricted = await asUser();
          expect('quantity' in restricted.fields, false,
            `a restricted column must be absent from fields: ${JSON.stringify(restricted.fields)}`);
          expect('name' in restricted.fields, true, 'the other columns stay listed');
          grok.dapi.domains.invalidateUiCaches();
          expect('quantity' in (await items().access()).fields, true, 'admin bypasses column security');
        } finally {
          await items().restoreColumnVisibility('quantity');
        }
        await items().revoke(probe.group, 'Edit');
        expect((await asUser()).can.edit, false, 'revoke must flip can.edit back');
        return 'checked';
      } finally {
        try {
          await items().revoke(probe.group, 'Edit');
        } catch (_) { /* already revoked */ }
        grok.dapi.domains.invalidateUiCaches();
      }
    });
    expect(outcome, 'checked',
      'the grant flip and column restriction were NOT verified: no restricted session (see the console for why)');
  });

  test('a table grant does not reach the rows of a defaultRowVisibility:"none" table', async () => {
    // access() ANDs can.edit/can.delete with "a grant on the securing TABLE reaches
    // its ROWS" — false exactly for this fixture. What this pins is the SERVER
    // semantics the term mirrors, side by side with apitests.item whose rows
    // default to table visibility.
    const e = await thrown(() => hidden().access());
    if (e instanceof DG.DomainValidationError) {
      console.log('skipped: apitests.hidden_item is not deployed (schema.json < 1.2.0)');
      return;
    }
    expect(e, null, `apitests.hidden_item access failed: ${e?.message}`);
    // Every server call is labelled: this test touches two tables, two sessions
    // and four grants, and a bare 'forbidden' would not say which.
    const step = async <T>(label: string, action: () => Promise<T>): Promise<T> => {
      try {
        return await action();
      } catch (x: any) {
        throw new Error(`${label}: ${x?.message ?? x} (${JSON.stringify(x?.body ?? {})})`);
      }
    };
    const outcome = await withRestrictedUser('wo4fvis', async (probe) => {
      const group = probe.group;
      const sku = `SKU-VIS-${probe.login}`;
      let visibleId: string | undefined;
      let hiddenId: string | undefined;
      try {
        visibleId = (await step('insert into apitests.item', () => items().insert({sku, name: 'Visibility probe'})))[0].id;
        hiddenId = (await step('insert into apitests.hidden_item',
          () => hidden().insert({sku, name: 'Visibility probe'})))[0].id;
        await step('grant View on apitests.item', () => items().grant(group, 'View'));
        await step('grant View on apitests.hidden_item', () => hidden().grant(group, 'View'));
        // Same grant, same row-mode: only the table whose rows default to table
        // visibility lets it through.
        const seenVisible = await step('read apitests.item as the restricted user',
          () => probe.asUser(() => items().query({filter: `sku = "${sku}"`})));
        const seenHidden = await step('read apitests.hidden_item as the restricted user',
          () => probe.asUser(() => hidden().query({filter: `sku = "${sku}"`})));
        expect(seenVisible.length, 1, 'a table View grant must reach rows of a table-visibility table');
        expect(seenHidden.length, 0,
          `a table View grant must NOT reach rows of a defaultRowVisibility:"none" table: ${JSON.stringify(seenHidden)}`);
        return 'checked';
      } finally {
        // Best-effort, one report per failure: a throw here would mask the result
        // AND skip the rest of the cleanup.
        const cleanup: [string, () => Promise<any>][] = [
          ['revoke View on apitests.item', () => items().revoke(group, 'View')],
          ['revoke View on apitests.hidden_item', () => hidden().revoke(group, 'View')],
          ['delete the apitests.item probe row', async () => visibleId && await items().delete(visibleId)],
          ['delete the apitests.hidden_item probe row', async () => hiddenId && await hidden().delete(hiddenId)],
        ];
        for (const [label, action] of cleanup)
          try {
            await action();
          } catch (x) {
            console.error(`cleanup — ${label} failed: ${x}`);
          }
        grok.dapi.domains.invalidateUiCaches();
      }
    });
    expect(outcome, 'checked',
      'the reaches-rows rule was NOT verified: no restricted session (see the console for why)');
  });
}, {owner: 'askalkin@datagrok.ai'});
