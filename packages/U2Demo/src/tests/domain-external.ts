/* The read UI over an externally bound table (EMS external bindings, phase C WO-C2): the same
   `domains.table(...).app()` Stockroom runs, over `northwind.order`, with every difference driven by
   the flags the server declares — no captions in the query, rows addressed by their canonical id
   (`10248`, `10248,11`), no New / Import / Bulk edit / Trash / History, no live probe, a `basic`
   filter offer — and every refusal that still reaches the UI named. Fixture: the `northwind`
   schema of the NorthwindBinding scratch package (later ApiSamples), read-only and never written;
   the category skips itself when it is not registered. */
import * as grok from 'datagrok-api/grok';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Filters, Rows} from '@datagrok-libraries/u2';
import type {DomainSource, FilterSchema, FilterValue} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const SCHEMA = 'northwind';
const ORDERS = `${SCHEMA}.order`;
const LINES = `${SCHEMA}.order_detail`;
const BASE = `/domains/${SCHEMA}/order`;

category('U2: domain external', () => {
  let skip: string | null = null;
  let orders: DomainTable;
  let app: DomainApp;

  async function settled(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
  }

  async function ready(src: DomainSource): Promise<void> {
    await settled(src);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  async function until(check: () => boolean, what: string): Promise<void> {
    for (let i = 0; i < 100 && !check(); i++)
      await delay(50);
    if (!check())
      throw new Error(`timed out waiting for ${what}`);
  }

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  before(async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === SCHEMA)) {
      skip = `the ${SCHEMA} schema is not registered`;
      return;
    }
    orders = await domains.table(ORDERS);
    app = domains.app({table: orders, base: BASE});
    document.body.append(app.root);
    await ready(app.listSource);
  });

  after(async () => {
    app?.dispose();
    app?.root.remove();
  });

  test('the handle declares the read-only external shape', async () => {
    if (skipped())
      return;
    const support = orders.table.support;
    expect(support.captions, false);
    expect(support.filters, 'basic');
    expect(support.transaction, false);
    expect(support.updateWhere, false);
    expect(support.audit, false);
    expect(support.probe, false);
    expect(orders.info.rowAddress, 'id');
    expect('batch' in support, true);
    expect(orders.table.probe === undefined, true, 'no probe member: a live source over it is not polled');
    expect(orders.table.audit === undefined, true);
    expect(orders.table.updateWhere === undefined, true);
  });

  test('the list reads without captions; the ref id is what a row carries', async () => {
    if (skipped())
      return;
    expect(app.page.value, 'list');
    expect(app.listSource.rows.items.value.length > 0, true);
    const df = app.listSource.df.value!;
    expect(df.columns.names().some((n) => n.startsWith(Rows.CAPTION_PREFIX)), false, 'no ~caption_ column');
    const first = app.listSource.rows.items.value[0];
    expect(first[Rows.caption('customer_id')] === undefined, true, 'no ~caption_ column on an external frame');
    expect(typeof first.customer_id, 'string');
    expect(app.listSource.readOnly.value, true, 'no transaction lands here');
  });

  test('New, Import, Bulk edit and Trash are absent, not disabled', async () => {
    if (skipped())
      return;
    const ribbon = app.ribbon().main;
    const el = (at: number): HTMLElement => {
      const item = ribbon[at];
      return item instanceof HTMLElement ? item : item.root;
    };
    await delay(50);
    expect(el(0).hidden, true, 'no New');
    expect(app.menuActions().map((a) => a.name).join(), '', 'no Import, no Bulk edit, no Trash');
    expect(el(3).hidden, true, 'and no ⋯ button');
    const row = app.listSource.rows.items.value[0];
    expect(app.list.actionsFor(row).some((a) => a.name === 'Delete'), false, 'no Delete on a row');
  });

  test('open 10248: the entity page by id, the address ends with /10248, the child tab lists its lines', async () => {
    if (skipped())
      return;
    expect(await app.open(`${BASE}/10248`), true);
    expect(app.page.value, 'entity');
    expect(app.entity.value, '10248');
    const source = app.entitySource.value!;
    await ready(source);
    await until(() => source.currentRow.value !== null, 'the row');
    expect(source.currentRow.value!.id, '10248');
    expect(app.keyOf(source.currentRow.value!), '10248');
    expect(app.path.value, `${BASE}/10248`);
    expect(app.panes.querySelector('[data-u2="domain-history"]'), null, 'no History pane');
    const children = app.panes.querySelector('[data-u2="domain-children"]');
    expect(children !== null, true, 'the child tab is there');
    // the child pane is a d4 grid on a canvas: its rows are not DOM; the empty marker hides once lines are in
    await until(() => children!.querySelector('.u2-domain-children-empty')?.hasAttribute('hidden') === true, 'the order lines');
    const lines = await grok.dapi.domains.table(LINES).query({filter: 'order_id = "10248"'});
    expect(lines.length > 0, true);
    expect(lines.every((line) => line.id.startsWith('10248,')), true, 'every line is 10248,<product>');
    expect(await app.goTo('list'), true);
  });

  test('order_detail 10248,11: the shared link is its id as one encoded segment, and it opens', async () => {
    if (skipped())
      return;
    const lines = await domains.table(LINES);
    expect(lines.info.rowAddress, 'id');
    const segment = encodeURIComponent('10248,11');
    const values = await grok.dapi.domains.table(LINES).get('10248,11');
    const link = lines.handler.deepLink(lines.handler.rowFrom(values));
    expect(link, `/domains/${SCHEMA}/order_detail/${segment}`, 'the link spells the id, not a dash-joined key');
    const a = domains.app({table: lines, base: `/domains/${SCHEMA}/order_detail`, children: false});
    document.body.append(a.root);
    try {
      await ready(a.listSource);
      expect(await a.open(link!), true);
      const source = a.entitySource.value!;
      await ready(source);
      await until(() => source.currentRow.value !== null, 'the line');
      expect(source.currentRow.value!.id, '10248,11');
      expect(source.currentRow.value!.product_id, 11);
      expect(a.keyOf(source.currentRow.value!), '10248,11');
      expect(a.path.value, `/domains/${SCHEMA}/order_detail/${segment}`);
    } finally {
      a.dispose();
      a.root.remove();
    }
  });

  test('a datetime key: the link encodes the canonical id once more, and it opens', async () => {
    if (skipped())
      return;
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'extwlive')) {
      console.log('skipped: the extwlive schema is not registered');
      return;
    }
    const id = '2026-09-18T07%3A08%3A09.123Z';
    const keyed = await domains.table('extwlive.dtkeyed');
    const values = await grok.dapi.domains.table('extwlive.dtkeyed').get(id);
    expect(values.id, id, 'the id as the server encodes it');
    const link = keyed.handler.deepLink(keyed.handler.rowFrom(values));
    expect(link, `/domains/extwlive/dtkeyed/${encodeURIComponent(id)}`);
    const a = domains.app({table: keyed, base: '/domains/extwlive/dtkeyed', children: false});
    document.body.append(a.root);
    try {
      await ready(a.listSource);
      expect(await a.open(link!), true);
      const source = a.entitySource.value!;
      await ready(source);
      await until(() => source.currentRow.value !== null, 'the row');
      expect(source.currentRow.value!.id, id);
      expect(a.keyOf(source.currentRow.value!), id);
    } finally {
      a.dispose();
      a.root.remove();
    }
  });

  test('live: true issues no probe — no version, no aggregate', async () => {
    if (skipped())
      return;
    const client = grok.dapi.domains.table(ORDERS);
    const proto = Object.getPrototypeOf(client);
    const calls: string[] = [];
    const spied = ['version', 'aggregate'].filter((name) => typeof proto[name] === 'function');
    const originals = Object.fromEntries(spied.map((name) => [name, proto[name]]));
    for (const name of spied)
      proto[name] = function(...args: unknown[]) { calls.push(name); return originals[name].apply(this, args); };
    const src = orders.source({pageSize: 5, live: true, liveMs: 100});
    try {
      await ready(src);
      await delay(350);
      expect(calls.join(), '', 'nothing was asked');
      expect(src.state.value, 'ready');
      expect(src.error.value === undefined, true);
    } finally {
      src.dispose();
      for (const name of spied)
        proto[name] = originals[name];
    }
  });

  test('a cold ?trash=1 is refused by name', async () => {
    if (skipped())
      return;
    expect(await app.open(`${BASE}?trash=1`), true);
    expect(app.mode.value, 'trash');
    await settled(app.listSource);
    expect(app.listSource.state.value, 'error');
    expect(String(app.listSource.error.value).includes('does not support restoring deleted rows'), true,
      String(app.listSource.error.value));
    expect(app.status.value.problem !== null, true, 'said on the status line');
    expect(await app.open(BASE), true);
    await ready(app.listSource);
  });

  test('the filter offer on ship_country lacks under and matches; a typed equality filters', async () => {
    if (skipped())
      return;
    await until(() => app.filters.input.value !== null, 'the filter box');
    const prop = app.listSource.schema.properties.find((p) => p.name === 'ship_country')!;
    const schema = (app.filters.input.value as unknown as {options: {schema: FilterSchema}}).options.schema;
    const offered = schema.operators!(prop, Filters.operators.for(prop)).map((o) => o.id);
    for (const id of ['under', 'matches', '!matches', '!like'])
      expect(offered.includes(id), false, `${id} is not offered`);
    expect(offered.includes('='), true);
    expect(offered.includes('like'), true);
    app.listSource.query.value = 'ship_country = "France"';
    await delay(100);
    await ready(app.listSource);
    const rows = app.listSource.rows.items.value;
    expect(rows.length > 0, true);
    expect(rows.every((r) => r.ship_country === 'France'), true);
    expect(app.status.value.problem, null);
    app.listSource.query.value = '';
    await delay(100);
    await ready(app.listSource);
  });

  test('every operator still offered on a string, a datetime and a float column runs', async () => {
    if (skipped())
      return;
    await until(() => app.filters.input.value !== null, 'the filter box');
    const schema = (app.filters.input.value as unknown as {options: {schema: FilterSchema}}).options.schema;
    const samples: Record<string, FilterValue> = {
      ship_country: 'France', ordered_on: new Date('1997-07-01T00:00:00Z'), freight: 32.38};
    const lists: Record<string, FilterValue> = {ship_country: ['France', 'Germany'], freight: [32.38, 11.61]};
    const valueFor = (column: string, op: string): FilterValue | undefined => {
      if (op === 'in' || op === 'not in')
        return lists[column];
      if (op === 'like' || op === 'starts' || op === 'ends')
        return 'Fr';
      return op.endsWith('null') ? undefined : samples[column];
    };
    const failed: string[] = [];
    for (const column of Object.keys(samples)) {
      const prop = app.listSource.schema.properties.find((p) => p.name === column)!;
      const offered = schema.operators!(prop, Filters.operators.for(prop)).map((o) => o.id);
      for (const op of offered) {
        app.listSource.query.value = Filters.group('and', [Filters.cond(column, op, valueFor(column, op))]);
        await delay(100);
        await settled(app.listSource);
        if (app.listSource.state.value !== 'ready' || app.status.value.problem !== null)
          failed.push(`${column} ${op}: ${String(app.listSource.error.value ?? app.status.value.problem)}`);
      }
      console.log(`offered on ${column}: ${offered.join(', ')}`);
    }
    expect(failed.join('; '), '', 'every offered operator runs');
    app.listSource.query.value = '';
    await delay(100);
    await ready(app.listSource);
  });
});
