/**
 * Integration tests for NodeDomainsDataSource (`grok s domains`). Requires a running Datagrok
 * server in ~/.grok/config.yaml whose dev key holds the CreateDomainSchema privilege (admin).
 * A throwaway user-managed schema is created and purged in afterAll.
 *
 * Run: npm run test:integration   (HOST=<alias> to target a specific server)
 */
import {afterAll, beforeAll, describe, expect, it} from 'vitest';
import {NodeApiClient, NodeDapi} from '../utils/node-dapi';
import {getDevKey} from '../utils/test-utils';

const HOST = process.env['HOST'] ?? '';
const SCHEMA = `clitest_${Date.now().toString(36)}`;
const TABLE = `${SCHEMA}.widget`;

let dapi: NodeDapi;
let offline = false;
let created = false;

beforeAll(async () => {
  try {
    const {url, key} = getDevKey(HOST);
    dapi = new NodeDapi(await NodeApiClient.login(url, key));
    await dapi.domains.createSchema(SCHEMA, 'CLI test', 'grok s domains integration test');
    created = true;
    await dapi.domains.applySchema(SCHEMA, {tables: {widget: {
      businessKey: ['sku'],
      columns: {
        sku: {type: 'string', required: true, unique: true},
        name: {type: 'string', isName: true},
        quantity: {type: 'int', min: 0},
      },
    }}});
  } catch (err) {
    offline = true;
    console.warn(`Domains integration tests skipped: ${(err as Error).message}`);
  }
});

afterAll(async () => {
  if (created)
    await dapi.domains.deleteSchema(SCHEMA);
});

describe('NodeDomainsDataSource (integration)', () => {
  it('lists the schema with its table', async () => {
    if (offline) return;
    const s = await dapi.domains.schema(SCHEMA);
    expect(s.managedBy).toBe('user');
    expect(s.tables.map((t: any) => t.name)).toEqual(['widget']);
    const manifest = await dapi.domains.manifest(SCHEMA);
    expect(Object.keys(manifest.tables.widget.columns)).toEqual(['sku', 'name', 'quantity']);
  });

  it('inserts, reads, updates with optimistic concurrency, and deletes a row', async () => {
    if (offline) return;
    const [ins] = await dapi.domains.insert(SCHEMA, 'widget', {sku: 'W-1', name: 'One', quantity: 1});
    expect(ins.created).toBe(true);
    const row = await dapi.domains.getRow(SCHEMA, 'widget', ins.id);
    expect(row).toMatchObject({sku: 'W-1', quantity: 1, version: 1});
    const upd = await dapi.domains.update(SCHEMA, 'widget', ins.id, {quantity: 2}, 1);
    expect(upd.version).toBe(2);
    await expect(dapi.domains.update(SCHEMA, 'widget', ins.id, {quantity: 3}, 1)).rejects.toThrow(/Version conflict/);
    await expect(dapi.domains.insert(SCHEMA, 'widget', {sku: 'W-2', quantity: -1})).rejects.toMatchObject({
      apiError: {body: {error: 'validation'}},
    });
    await dapi.domains.deleteRow(SCHEMA, 'widget', ins.id);
    expect(await dapi.domains.getRow(SCHEMA, 'widget', ins.id)).toBeNull();
  });

  it('uploads csv, upserts, queries, counts, aggregates, and bulk-deletes', async () => {
    if (offline) return;
    const report = await dapi.domains.batch(SCHEMA, 'widget', Buffer.from('sku,name,quantity\nW-10,Ten,10\nW-11,Eleven,11\n'), 'text/csv');
    expect(report).toMatchObject({inserted: 2, errorCount: 0});
    const upsert = await dapi.domains.batch(SCHEMA, 'widget',
      Buffer.from(JSON.stringify([{sku: 'W-10', quantity: 100}, {sku: 'W-12', quantity: 12}])), 'application/json', {mode: 'upsert'});
    expect(upsert).toMatchObject({inserted: 1, updated: 1});

    const rows = await dapi.domains.query(SCHEMA, 'widget', {filter: 'quantity > 11', sort: '!quantity', columns: ['sku', 'quantity']});
    expect(rows.map((r) => r.sku)).toEqual(['W-10', 'W-12']);
    expect(await dapi.domains.count(SCHEMA, 'widget')).toBe(3);
    const agg = await dapi.domains.aggregate(SCHEMA, 'widget', {measures: [{fn: 'sum', column: 'quantity', as: 'total'}]});
    expect(agg[0].total).toBe(123);
    const d42 = await dapi.domains.queryD42(SCHEMA, 'widget', {limit: 1});
    expect(d42.length).toBeGreaterThan(0);

    const del = await dapi.domains.deleteWhere(SCHEMA, 'widget', 'quantity > 11', 1);
    expect(del).toEqual({deleted: 1, hasMore: true});
    expect(await dapi.domains.count(SCHEMA, 'widget')).toBe(2);
    const audit = await dapi.domains.tableAudit(SCHEMA, 'widget', 5);
    expect(audit[0].op).toBe('delete');
  });

  it('runs a transaction with a $ref', async () => {
    if (offline) return;
    const res = await dapi.domains.transaction(SCHEMA, [
      {op: 'insert', table: 'widget', ref: 'a', values: {sku: 'W-30', name: 'Thirty'}},
      {op: 'update', table: 'widget', id: '$a', values: {quantity: 30}},
    ]);
    expect(res[1]).toMatchObject({id: res[0].id, version: 2});
  });

  it('grants and revokes on the table entity', async () => {
    if (offline) return;
    const [group] = await dapi.groups.by(1).filter('personal = false').list();
    if (!group) return;
    const entityId = await dapi.domains.entityId({schema: SCHEMA, table: 'widget'});
    await dapi.domains.grant(entityId, group.id, 'Edit');
    const grants = await dapi.domains.grants(entityId);
    expect(grants.some((g: any) => g.group.id === group.id && g.permission === 'Edit')).toBe(true);
    await dapi.domains.revoke(entityId, group.id);
    expect((await dapi.domains.grants(entityId)).some((g: any) => g.group.id === group.id)).toBe(false);
    const caps = await dapi.domains.capabilities(SCHEMA, 'widget');
    expect(caps).toMatchObject({canView: true, securityMode: 'table', hasBusinessKey: true});
  });
});
