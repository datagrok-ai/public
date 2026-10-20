import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, awaitCheck} from '@datagrok-libraries/test/src/test';

// Tests for grok.dapi.domains against the 'apitests' domain schema that this
// package declares in databases/apitests/schema.json (deployed on publish).
category('Dapi: domains', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const sku = () => `SKU-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  test('schemas listing', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    const s = schemas.find((x) => x.name === 'apitests');
    expect(s != null, true, 'apitests schema not registered');
    expect(s!.pgSchema, 'apitests');
    expect(s!.tables.some((t) => t.name === 'item'), true, 'item table not registered');
    expect(s!.tables.find((t) => t.name === 'item')!.securityMode, 'row');
  });

  test('insert, query, get, update, audit, delete', async () => {
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Widget', quantity: 5, note: 'first'});
    expect(ins.created, true);
    const id: string = ins.id;
    try {
      const rows = await items().query({filter: `sku = "${key}"`});
      expect(rows.length, 1);
      expect(rows[0].name, 'Widget');
      expect(rows[0].note, 'first', 'jsonb property schema key not returned');

      const row = await items().get(id);
      expect(row.sku, key);
      expect(row.version, 1);

      const upd = await items().update(id, {name: 'Widget 2'}, {version: 1});
      expect(upd.version, 2);

      const audit = await items().audit(id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].before.name, 'Widget');
      expect(audit[1].after.name, 'Widget 2');
    } finally {
      await items().delete(id);
    }
    expect(await items().get(id), null, 'soft-deleted row is still visible');
  });

  test('optimistic concurrency conflict', async () => {
    const [ins] = await items().insert({sku: sku(), name: 'A'});
    try {
      await items().update(ins.id, {name: 'B'}, {version: 1});
      let conflict = '';
      try {
        await items().update(ins.id, {name: 'C'}, {version: 1});
      } catch (e: any) {
        conflict = e.message ?? `${e}`;
      }
      expect(conflict.includes('Version conflict'), true, `unexpected error: '${conflict}'`);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('duplicate business key reported per row', async () => {
    const key = sku();
    const [first] = await items().insert({sku: key});
    try {
      const [dup] = await items().insert({sku: key});
      expect(dup.status, 'duplicate');
      expect(dup.existingId, first.id);
    } finally {
      await items().delete(first.id);
    }
  });

  test('promote makes a row an entity', async () => {
    const [ins] = await items().insert({sku: sku()});
    try {
      const res = await items().promote(ins.id);
      expect(res.promoted, true);
      const audit = await items().audit(ins.id);
      expect(audit[audit.length - 1].op, 'promote');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('dapi2 generated client: domains removed at parity, init export intact', async () => {
    // WO-4b: every dapi2 function belonged to the domains namespace — with it
    // removed at typed-surface parity the generated client has no value members
    // left (the chats namespace never emitted functions), so `dapi2` is
    // type-only now; dapi2Init survives as a value and the OpenAPI yaml keeps
    // every /domains/ route.
    expect(typeof grok.dapi2Init, 'function');
    // NB: compare against true — utils expect() treats a passed undefined
    // `expected` as its default (true), so expect(x, undefined) never passes.
    expect((grok as any).dapi2?.domains === undefined, true, 'dapi2.domains must be gone');
  });

  test('table name validation', async () => {
    let error = '';
    try {
      grok.dapi.domains.table('noseparator');
    } catch (e: any) {
      error = e.message;
    }
    expect(error.includes('<schema>.<table>'), true);
  });

  test('read scope: count/exists take one object, the positional form still answers it', async () => {
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Scoped'});
    const filter = `sku = "${key}"`;
    try {
      // Every form of the same scope must answer the same number — that is what makes a
      // paged list's total agree with its rows.
      expect(await items().count({filter}), 1, 'the scope object');
      expect(await items().count(filter), 1, 'the deprecated positional filter');
      expect(await items().count({search: 'Scoped', filter}), 1, 'search inside the scope');
      expect(await items().exists({filter}), true, 'exists over the scope object');
      // An empty scope and no scope are the whole table, and a filter object is a FILTER,
      // not a scope — the discriminator's two edges.
      const all = await items().count();
      expect(await items().count({}), all, 'an empty scope is the whole table');
      expect(await items().count({property: 'sku', operator: '=', value: key} as any), 1,
        'a condition object must be read as a filter, not as a scope');
      await items().delete(ins.id);
      expect(await items().count({filter}), 0, 'a deleted row leaves the default scope');
      expect(await items().count({filter, deleted: 'only'}), 1, 'the trash scope, as one object');
      expect(await items().count(filter, {deleted: 'only'}), 1, 'the deprecated pair agrees with it');
      expect(await items().exists({filter, deleted: 'only'}), true, 'exists follows the same scope');
      // aggregate() takes the same scope, so a summary and a count cannot disagree.
      const [row] = await items().aggregate<string, string>(
        {measures: [{fn: 'count'}], filter, deleted: 'only'});
      expect(Number(row.count), 1, `aggregate must honour the scope: ${JSON.stringify(row)}`);
    } finally {
      try {
        await items().delete(ins.id);
      } catch (_) { /* already deleted by the scope case */ }
    }
  });

  test('captions: a ref column projects the target row\'s display name', async () => {
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Caption parent'});
    const events = grok.dapi.domains.table('apitests.item_event');
    try {
      const [ev] = await events.insert({item_id: ins.id, kind: 'captioned', amount: 1});
      const filter = `id = "${ev.id}"`;
      const caption = DG.domainCaptionColumn('item_id');
      expect(caption, '~caption_item_id', 'the caption column name is spelled once, in the js-api');
      const [row] = await events.query({filter, captions: ['item_id']});
      // the name column wins over the business key (C3): apitests.item names itself by `name`
      expect(row[caption], 'Caption parent', `the caption must be the parent's display name: ${JSON.stringify(row)}`);
      // A caption is not an expand and does not drag its ref column in.
      const [narrow] = await events.query({filter, columns: ['kind'], captions: ['item_id']});
      expect(narrow[caption], 'Caption parent', 'a caption is independent of the projected columns');
      expect('item_id' in narrow, false,
        `asking for a caption must not project the ref column: ${JSON.stringify(narrow)}`);
      const [plain] = await events.query({filter});
      expect(caption in plain, false, `captions are never on by default: ${JSON.stringify(plain)}`);
    } finally {
      await items().delete(ins.id); // cascades item_event
    }
  });

  test('captions: the queryDf column is typed and export-tagged, and toCsv drops it', async () => {
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Caption frame parent'});
    const events = grok.dapi.domains.table('apitests.item_event');
    try {
      const [ev] = await events.insert({item_id: ins.id, kind: 'captioned-df', amount: 2});
      const df = await events.queryDf({filter: `id = "${ev.id}"`, captions: ['item_id']});
      expect(df.rowCount, 1);
      const col = df.col('~caption_item_id');
      expect(col != null, true, `no caption column: ${df.columns.names().join(', ')}`);
      expect(col!.type, DG.TYPE.STRING, `a caption is a string column, got ${col!.type}`);
      expect(col!.get(0), 'Caption frame parent', 'the caption column carries the parent display name');
      expect(col!.meta.includeInCsvExport, false, 'the caption leaked into csv export');
      expect(col!.meta.includeInBinaryExport, false, 'the caption leaked into binary export');
      expect(df.toCsv().includes('~caption_'), false, 'the caption column leaked into toCsv()');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('captions: no oracle — unknown, non-ref, nested and duplicate all refuse', async () => {
    const events = grok.dapi.domains.table('apitests.item_event');
    const refused = async (captions: string[]): Promise<any> => {
      try {
        await events.query({filter: 'amount > -1', limit: 1, captions: captions as any});
      } catch (e: any) {
        return e;
      }
      return null;
    };
    // An unknown column and a real-but-not-a-ref column must be indistinguishable: the
    // refusal may not tell a caller whether the column exists.
    const unknown = await refused(['nosuchcolumn']);
    const nonRef = await refused(['kind']);
    for (const [label, e] of [['unknown', unknown], ['non-ref', nonRef]] as [string, any][]) {
      expect(e instanceof DG.DomainFilterError, true,
        `${label} caption: expected DomainFilterError, got ${e?.constructor?.name}: ${e?.message}`);
      expect(e.message.includes('Unknown or inaccessible caption column'), true,
        `${label} caption: unexpected message '${e.message}'`);
    }
    expect(unknown.message.replace('nosuchcolumn', 'X'), nonRef.message.replace('kind', 'X'),
      'the two refusals differ, so the message is an oracle for column existence');
    const nested = await refused(['item_id.sku']);
    expect(nested?.message.includes('Nested caption'), true,
      `a dotted caption must refuse as nested: ${nested?.message}`);
    const duplicate = await refused(['item_id', 'item_id']);
    expect(duplicate?.message.includes('Duplicate caption'), true,
      `a repeated caption must refuse: ${duplicate?.message}`);
  });

}, {owner: 'askalkin@datagrok.ai'});

// Registry reflection (ui-js-api WO-2): grok.dapi.domains.registry — the runtime
// Property metadata, table info with FK-inverted child tables, and batched
// display-name resolution. Assertions run against this package's own 'apitests'
// schema; grit.issue assertions (the dogfood schema with refs/min/nameColumn)
// skip cleanly where Grit is not deployed.
category('Dapi: domain registry', () => {
  const registry = () => grok.dapi.domains.registry;
  const sku = () => `SKU-REG-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function gritDeployed(): Promise<boolean> {
    const schemas = await grok.dapi.domains.schemas.list();
    return schemas.some((s) => s.name === 'grit');
  }

  test('rowProperties: schema.json constraints round-trip', async () => {
    const props = await registry().rowProperties('apitests.item');
    const by = (n: string) => props.find((p) => p.name === n);
    expect(by('sku') != null, true, `sku property missing: ${props.map((p) => p.name).join(',')}`);
    expect(by('sku')!.nullable, false, 'required column must not be nullable');
    expect(by('name')!.nullable, true, 'optional column must be nullable');
    expect(by('quantity')!.min, 0, 'min constraint from schema.json');
    expect(by('quantity')!.propertyType, 'int');
  });

  test('rowProperties: grit.issue matches schema.json (refs/min/nullable)', async () => {
    if (!await gritDeployed()) {
      console.log('skipped: Grit is not deployed');
      return;
    }
    const props = await registry().rowProperties('grit.issue');
    const by = (n: string) => props.find((p) => p.name === n)!;
    expect(by('number').min, 1, 'min from schema.json');
    expect(by('number').nullable, true, 'auto-numbered column is nullable (server fills it)');
    expect(by('title').nullable, false);
    expect(by('description').nullable, true);
    // Lookup tables since Grit 2.0.0: status/priority/type are ref columns, so the
    // allowed values live in the target table and the Property carries its semType.
    for (const [column, target] of [['project_id', 'grit.project'], ['status_id', 'grit.status'],
      ['priority_id', 'grit.priority'], ['type_id', 'grit.issue_type']] as [string, string][])
      expect(by(column).semType, target, `${column} must carry the target row semType`);
    expect(by('project_id').nullable, false, 'a required ref must not be nullable');
    expect(by('status_id').nullable, true, 'an optional ref must be nullable');
    expect(by('project_id').friendlyName, 'Project',
      'ref columns must carry the label the platform renders, not the wire name');
    expect(by('status_id').friendlyName, 'Status', 'the schema.json friendlyName must win');
  });

  test('rowProperties: unknown table rejects with a typed validation error', async () => {
    let e: any = null;
    try {
      await registry().rowProperties('apitests.nosuch');
    } catch (x) {
      e = x;
    }
    expect(e instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
  });

  test('tableInfo: identity, security, and FK-inverted childTables', async () => {
    const info = await registry().tableInfo('apitests.item');
    expect(JSON.stringify(info.businessKey), JSON.stringify(['sku']));
    // NB: nameColumn is not asserted here — deployed registries may carry an
    // isName drift on apitests.item (same manifest version, no reapply); the
    // grit.issue test pins nameColumn against its committed schema.json.
    expect(info.securityMode, 'row');
    expect(info.audit, true);
    expect(info.singularName, 'item', 'effective singular derived from the table name');
    expect(info.pluralName, 'items');
    const child = info.childTables.find((c) => c.table === 'item_event');
    expect(child != null, true, `item_event missing from childTables: ${JSON.stringify(info.childTables)}`);
    expect(child!.schema, 'apitests');
    expect(child!.fkColumn, 'item_id');
    expect(child!.label, 'Item', 'friendly FK label (the _id suffix dropped)');
  });

  test('tableInfo: grit.issue childTables list comment', async () => {
    if (!await gritDeployed()) {
      console.log('skipped: Grit is not deployed');
      return;
    }
    const info = await registry().tableInfo('grit.issue');
    expect(info.nameColumn, 'title');
    expect(info.childTables.some((c) => c.table === 'comment' && c.fkColumn === 'issue_id'), true,
      `comment missing from childTables: ${JSON.stringify(info.childTables)}`);
  });

  test('resolveNames: display-identity chain, null for unresolvable ids', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const info = await registry().tableInfo('apitests.item');
    const named = sku();
    const bare = sku();
    const ghost = '00000000-0000-0000-0000-000000000000';
    let inserted: {id: string}[] = [];
    try {
      // One row with a name value, one without: the second proves the
      // business-key fallback regardless of whether the registry declares a
      // name column for apitests.item (deployed registries drift on isName).
      // A single insert inside the try so no partial pair can leak.
      inserted = await items.insert([{sku: named, name: 'Resolve probe'}, {sku: bare}]);
      const [insNamed, insBare] = inserted;
      const names = await registry().resolveNames('apitests.item', [insNamed.id, insBare.id, ghost]);
      expect(names[insNamed.id], info.nameColumn != null ? 'Resolve probe' : named,
        `display identity must follow the declared name column (${info.nameColumn})`);
      expect(names[insBare.id], bare, 'empty name — the business key is the display identity');
      expect(Object.keys(names).includes(ghost), true, 'every requested id must be a key');
      expect(names[ghost] == null, true, 'unresolvable ids must map to null');
    } finally {
      for (const r of inserted)
        await items.delete(r.id);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

// DomainQuery state class (ui-js-api WO-6): the single serializable representation
// of what a user is looking at — the `DomainQuery` function's parameters, a URL deep
// link, and a REST spec, all the same object. UI-only state (view mode, current
// entity) never enters it: the reserved 'view='/'entity=' URL params are ignored.
category('Dapi: domain query state', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  /** Key-order-independent comparison of two plain objects. NB the replacer is the
   * TOP-level key list, so it also masks keys nested objects do not share — enough for
   * the flat parameter/URL maps here, not a deep comparison. */
  const sorted = (o: any) => JSON.stringify(o, Object.keys(o).sort());
  /** The message of whatever [action] throws, or null when it succeeds. */
  const thrown = (action: () => any): string | null => {
    try {
      action();
      return null;
    } catch (e: any) {
      return e.message ?? `${e}`;
    }
  };

  test('URL parameters round-trip losslessly, reserved params ignored', async () => {
    // The platform's list-element binding: one URL key per element (filters[0]), so
    // a recorded run can be re-parameterized element by element.
    // Every list parameter of the function, so nothing is lost silently (a link that
    // drops joins or groupBy is a different query). The fixture deliberately mixes
    // select-mode (columns/offset) and aggregate-mode (aggregations/groupBy) parameters:
    // the URL layer is mode-agnostic and nothing is executed here — the modes are
    // reconciled by toSpec()/the function, covered below and in the malformed test.
    const url = {
      'columns[0]': 'kind',
      'columns[1]': 'amount',
      'filters[0]': 'amount > 1',
      'filters[1]': '{"property":"kind","operator":"=","value":"Widget"}',
      'joins[0]': 'item_id',
      'aggregations[0]': 'avg(amount) as avg_amount',
      'groupBy[0]': 'kind',
      'orderBy[0]': '!created_on',
      'limit': '25',
      'offset': '10',
    };
    const q = DG.DomainQuery.fromUrlParams('apitests', 'item_event',
      {...url, 'view': 'grid', 'entity': 'abc'});     // UI-only state rides along, stays out
    expect(q.schema, 'apitests');
    expect(q.table, 'item_event');
    expect(q.columns!.join(','), 'kind,amount');
    expect(q.filters!.length, 2);
    expect(q.filters![1].includes('Widget'), true, `JSON element mangled: ${q.filters![1]}`);
    expect(q.joins!.join(','), 'item_id');
    expect(q.aggregations!.join(','), 'avg(amount) as avg_amount');
    expect(q.groupBy!.join(','), 'kind');
    expect(q.orderBy!.join(','), '!created_on');
    expect(q.limit, 25);
    expect(q.offset, 10);
    expect(sorted(q.toUrlParams()), sorted(url), 'URL round trip is not lossless');
    // toParams() is the same content in the function's own shape.
    expect(sorted(DG.DomainQuery.fromParams(q.toParams()).toUrlParams()), sorted(url),
      'params round trip is not lossless');

    // Gaps close and indices are read in ascending order (a hand-edited link stays usable).
    const gapped = DG.DomainQuery.fromUrlParams('apitests', 'item',
      {'filters[2]': 'b = 2', 'filters[0]': 'a = 1'});
    expect(gapped.filters!.join('|'), 'a = 1|b = 2');

    // A lone smart-filter string is passed to the REST spec verbatim (the server parses
    // that grammar); only SEVERAL of them need the function itself.
    const one = DG.DomainQuery.fromUrlParams('apitests', 'item', {'filters[0]': 'quantity > 1'});
    expect(one.toSpec().filter, 'quantity > 1', 'a single grammar element must pass through');
  });

  test('view state -> DomainQuery -> params, and run() reproduces the subset', async () => {
    const mine = `SKU-QS-${stamp()}`;
    let inserted: {id: string}[] = [];
    let view: any = null;
    let df: _DG.DataFrame | null = null;
    try {
      inserted = await items().insert({sku: mine, name: mine});
      view = DG.DomainView.create({schema: 'apitests', table: 'item',
        permanentFilter: `sku = "${mine}"`, embedded: true});
      grok.shell.addView(view);
      await awaitCheck(() => view.root.textContent!.includes(mine),
        'the filtered row never appeared in the view', 15000);
      const params = view.query;
      const q = DG.DomainQuery.fromParams(params);
      expect(sorted(q.toParams()), sorted(params), 'DomainView.query did not round-trip');
      expect((q.filters ?? []).some((f) => f.includes(mine)), true,
        `the view's filter is missing from its query: ${JSON.stringify(params)}`);
      df = await q.run();
      expect(df!.rowCount, 1, 'run() must reproduce the subset the view shows');
      expect(df!.col('sku')!.get(0), mine);
    } finally {
      if (df != null)
        grok.shell.closeTable(df);
      if (view != null)
        view.close();
      for (const r of inserted)
        await items().delete(r.id);
    }
  });

  test('run() matches queryDf(toSpec()) and records a creation script', async () => {
    const prefix = `SKU-QR-${stamp()}`;
    // JSON condition elements: values are bound server-side, and toSpec() AND-joins
    // them into one REST condition tree.
    const q = new DG.DomainQuery({schema: 'apitests', table: 'item',
      filters: [`{"property":"sku","operator":"like","value":"${prefix}%"}`,
        '{"property":"quantity","operator":">","value":1}'],
      orderBy: ['!quantity'], limit: 10});
    // Recording is gated on the user's data-history setting (data_history.dart), so a
    // profile with it off would fail the creation-script assertions: force it on and
    // restore whatever the profile had.
    const dataHistory = grok.shell.settings.dataHistory;
    let inserted: {id: string}[] = [];
    let df: _DG.DataFrame | null = null;
    try {
      grok.shell.settings.dataHistory = true;
      inserted = await items().insert([
        {sku: `${prefix}-1`, name: 'one', quantity: 1},
        {sku: `${prefix}-2`, name: 'two', quantity: 2},
        {sku: `${prefix}-3`, name: 'three', quantity: 3}]);
      const direct = await items().queryDf(q.toSpec());
      df = await q.run();
      expect(df!.rowCount, 2, 'the AND-combined filter elements did not select 2 rows');
      expect(df!.rowCount, direct.rowCount, 'run() and queryDf(toSpec()) disagree on row count');
      expect(df!.col('sku')!.toList().join(','), direct.col('sku')!.toList().join(','),
        'run() and queryDf(toSpec()) returned different rows');
      // The function arranges its output (system columns behind '~') — proof the run
      // went through the function, not the raw REST path.
      expect(df!.col('~id') != null, true, `system columns are not hidden: ${df!.columns.names()}`);
      // Recorded: the frame carries a creation script, so it refreshes, data-syncs,
      // and takes URL parameters.
      const script = df!.tags['.script'] ?? '';
      expect(script.includes('DomainQuery'), true, `no DomainQuery creation script: '${script}'`);
      expect(script.includes(prefix), true, `the creation script lost the filter values: '${script}'`);
    } finally {
      grok.shell.settings.dataHistory = dataHistory;
      if (df != null)
        grok.shell.closeTable(df);
      for (const r of inserted)
        await items().delete(r.id);
    }
  });

  test('fromBuilder preserves where/orderBy/top and selects the same rows', async () => {
    const prefix = `SKU-QB-${stamp()}`;
    let inserted: {id: string}[] = [];
    try {
      inserted = await items().insert([
        {sku: `${prefix}-1`, quantity: 1},
        {sku: `${prefix}-2`, quantity: 2}]);
      const builder = items().query()
        .where('sku', 'like', `${prefix}%`)
        .where('quantity', '>', 1)
        .orderBy('sku')
        .top(5);
      const q = DG.DomainQuery.fromBuilder(builder);
      expect(q.schema, 'apitests');
      expect(q.table, 'item');
      expect(q.filters!.length, 2, 'AND conjuncts must become one filter element each');
      expect(q.orderBy!.join(','), 'sku');
      expect(q.limit, 5);
      const rows = await builder;
      const df = await items().queryDf(q.toSpec());
      expect(df.rowCount, rows.length, 'the converted query selects different rows');
      expect(df.rowCount, 1);
    } finally {
      for (const r of inserted)
        await items().delete(r.id);
    }
  });

  test('malformed input throws instead of degrading', async () => {
    const bad = (params: {[key: string]: string}) =>
      thrown(() => DG.DomainQuery.fromUrlParams('apitests', 'item', params));
    expect((bad({'filters[k]': 'a = 1'}) ?? '').includes('Malformed URL parameter'), true,
      'a non-numeric element index must throw');
    expect((bad({'limit': 'ten'}) ?? '').includes('Malformed URL parameter'), true,
      'a non-integer limit must throw, not resolve to NaN');
    expect((bad({'limit': '-5'}) ?? '').includes('Malformed URL parameter'), true,
      'a negative limit must throw — the server clamps it to 0 and returns nothing');
    expect((bad({'filters': 'a = 1'}) ?? '').includes('filters[0]'), true,
      'a list bound without an index must name the element form');
    expect(bad({'view': 'grid', 'unknown[3]': 'x'}), null, 'unknown keys must be ignored');

    const q = new DG.DomainQuery({schema: 'apitests', table: 'item'});
    q.filters = ['{not json'];
    expect((thrown(() => q.toSpec()) ?? '').includes('invalid JSON'), true);
    // Sub-group members must be conditions, nested groups, or connectors — the same
    // shape check the function applies, so run() and toSpec() reject the same input.
    q.filters = ['["a","b"]'];
    expect((thrown(() => q.toSpec()) ?? '').includes('condition sub-group'), true,
      'a bare-string sub-group member must be rejected');
    // Several smart-filter strings can only be parsed by the function itself (a bare
    // string inside a REST condition tree means a connector, not a filter).
    q.filters = ['a = 1', 'b = 2'];
    expect((thrown(() => q.toSpec()) ?? '').includes('run()'), true);
    q.filters = undefined;
    q.aggregations = ['count'];
    expect(q.isAggregate, true);
    expect((thrown(() => q.toSpec()) ?? '').includes('aggregate'), true);
  });
}, {owner: 'askalkin@datagrok.ai'});
