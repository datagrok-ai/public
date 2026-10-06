import {describe, it, expect} from 'vitest';
import {parseKeyValues, parseMeasures, querySpec, applyBody, schemaRow, tableRow, columnRows, rowsToCsv} from '../commands/server-domains';
import {NodeDomainsDataSource, parseDomainAddress} from '../utils/node-dapi';

interface Call {method: string; path: string; body?: any; contentType?: string}

function makeMock(responder: (method: string, path: string, body?: any) => any) {
  const calls: Call[] = [];
  const client: any = {
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      return responder(method, path, body);
    },
    get(path: string) { return this.request('GET', path); },
    post(path: string, body?: any) { return this.request('POST', path, body); },
    del(path: string) { return this.request('DELETE', path); },
    async putBytes(path: string, bytes: Buffer, contentType: string) {
      calls.push({method: 'POST', path, body: bytes.toString(), contentType});
      return responder('POST', path, bytes);
    },
    async postForBytes(path: string, body: any) {
      calls.push({method: 'POST', path, body});
      return Buffer.from('d42');
    },
  };
  return {client, calls};
}

describe('parseDomainAddress', () => {
  it('splits schema.table on the first dot', () => {
    expect(parseDomainAddress('grit.issue')).toEqual({schema: 'grit', table: 'issue'});
    expect(parseDomainAddress('grit')).toEqual({schema: 'grit'});
  });

  it('enforces the level a verb needs', () => {
    expect(() => parseDomainAddress('grit', {table: true})).toThrow(/needs a table/);
    expect(() => parseDomainAddress('grit.issue', {table: false})).toThrow(/needs a schema/);
    expect(() => parseDomainAddress('grit.', {})).toThrow(/Invalid domain address/);
    expect(() => parseDomainAddress('', {})).toThrow(/Invalid domain address/);
  });
});

describe('parseKeyValues', () => {
  it('types values that parse as JSON and keeps the rest as strings', () => {
    expect(parseKeyValues(['title=Crash on save', 'quantity=5', 'done=true', 'note=null', 'tags=["a","b"]', 'sku=W-1']))
      .toEqual({title: 'Crash on save', quantity: 5, done: true, note: null, tags: ['a', 'b'], sku: 'W-1'});
  });

  it('keeps everything after the first = in the value', () => {
    expect(parseKeyValues(['expr=a=b'])).toEqual({expr: 'a=b'});
  });

  it('rejects an argument without a column name', () => {
    expect(() => parseKeyValues(['=x'])).toThrow(/Expected <column>=<value>/);
    expect(() => parseKeyValues(['novalue'])).toThrow(/Expected <column>=<value>/);
  });
});

describe('parseMeasures', () => {
  it('parses count, fn(column) and aliases', () => {
    expect(parseMeasures('count, sum(amount) as total, avg(amount)')).toEqual([
      {fn: 'count'},
      {fn: 'sum', column: 'amount', as: 'total'},
      {fn: 'avg', column: 'amount'},
    ]);
  });

  it('accepts expanded fk paths and rejects garbage', () => {
    expect(parseMeasures('max(project_id.number) AS last')).toEqual([{fn: 'max', column: 'project_id.number', as: 'last'}]);
    expect(() => parseMeasures('sum amount')).toThrow(/Cannot parse measure/);
  });
});

describe('querySpec', () => {
  it('only carries what was asked for', () => {
    expect(querySpec({}, '', undefined, 0)).toEqual({});
    expect(querySpec({sort: '!created_on', columns: 'a,b', expand: ['x', 'y,z']}, 'status = "open"', 20, 40))
      .toEqual({filter: 'status = "open"', sort: '!created_on', columns: ['a', 'b'], expand: ['x', 'y', 'z'], limit: 20, offset: 40});
  });
});

describe('applyBody', () => {
  it('keeps only the keys the server accepts and folds the flags in', () => {
    const manifest = {name: 'grit', version: '1.0.0', description: 'x', tables: {issue: {}}, dropTables: ['old']};
    expect(applyBody(manifest, {'confirm-destructive': true, 'if-version': 3}))
      .toEqual({tables: {issue: {}}, dropTables: ['old'], confirmDestructive: true, ifVersion: '3'});
    expect(applyBody({tables: {}}, {})).toEqual({tables: {}});
  });
});

describe('table rows', () => {
  it('fills the defaults the entity serializer omits', () => {
    expect(schemaRow({name: 'grit', version: '2.1.0', tables: [{}, {}], id: 'S'}))
      .toEqual({name: 'grit', friendlyName: '', managedBy: 'package', version: '2.1.0', tables: 2, id: 'S'});
    expect(tableRow({name: 'issue', businessKey: ['project_id', 'number'], id: 'T'}))
      .toEqual({name: 'issue', securityMode: 'table', businessKey: 'project_id,number', nameColumn: '', origin: 'package', readOnly: false, description: '', id: 'T'});
  });

  it('lists the manifest columns of a table', () => {
    const rows = columnRows({columns: {
      sku: {type: 'string', required: true, unique: true},
      status: {type: 'string', choices: ['open', 'closed'], default: 'open'},
      project_id: {type: 'ref', ref: 'project'},
    }});
    expect(rows.map((r) => r.column)).toEqual(['sku', 'status', 'project_id']);
    expect(rows[0]).toMatchObject({required: true, unique: true});
    expect(rows[1]).toMatchObject({choices: 'open|closed', default: 'open'});
    expect(rows[2]).toMatchObject({type: 'ref', ref: 'project'});
  });
});

describe('rowsToCsv', () => {
  it('uses the union of keys and serializes nested values', () => {
    const csv = rowsToCsv([{id: '1', name: 'a,b', tags: [{id: 'x'}]}, {id: '2', extra: null}]);
    expect(csv.split('\n')).toEqual(['id,name,tags,extra', '1,"a,b","[{""id"":""x""}]",', '2,,,', '']);
    expect(rowsToCsv([])).toBe('');
  });
});

describe('NodeDomainsDataSource', () => {
  const SCHEMAS = [{id: 'S1', name: 'grit', tables: [{id: 'T1', name: 'issue'}]}];

  it('resolves grant targets from the schema list', async () => {
    const {client} = makeMock((_m, path) => path.startsWith('/domains/schemas') ? SCHEMAS : null);
    const ds = new NodeDomainsDataSource(client);
    expect(await ds.entityId({schema: 'grit'})).toBe('S1');
    expect(await ds.entityId({schema: 'grit', table: 'issue'})).toBe('T1');
    await expect(ds.entityId({schema: 'grit', table: 'nope'})).rejects.toThrow(/not found/);
    await expect(ds.entityId({schema: 'nope'})).rejects.toThrow(/not found/);
  });

  it('hits the row endpoints with the right shapes', async () => {
    const {client, calls} = makeMock(() => ({}));
    const ds = new NodeDomainsDataSource(client);
    await ds.query('grit', 'issue', {filter: 'x', limit: 5});
    await ds.update('grit', 'issue', 'R1', {title: 't'}, 3);
    await ds.insert('grit', 'issue', {title: 't'}, true);
    await ds.deleteWhere('grit', 'issue', 'status = "closed"', 10);
    await ds.revoke('S1', 'G1', 'Edit');
    await ds.revoke('S1', 'G1');
    expect(calls.map((c) => [c.method, c.path])).toEqual([
      ['POST', '/domains/grit/issue/query'],
      ['PATCH', '/domains/grit/issue/R1'],
      ['POST', '/domains/grit/issue?errorOnDuplicate=true'],
      ['POST', '/domains/grit/issue/delete'],
      ['DELETE', '/domains/grants/S1?group=G1&permission=Edit'],
      ['DELETE', '/domains/grants/S1?group=G1'],
    ]);
    expect(calls[1].body).toEqual({values: {title: 't'}, version: 3});
    expect(calls[3].body).toEqual({filter: 'status = "closed"', limit: 10});
  });

  it('counts through the aggregate endpoint', async () => {
    const {client, calls} = makeMock(() => [{count: 42}]);
    const ds = new NodeDomainsDataSource(client);
    expect(await ds.count('grit', 'issue', 'status = "open"')).toBe(42);
    expect(calls[0].body).toEqual({measures: [{fn: 'count'}], filter: 'status = "open"'});
  });

  it('answers null for a 404 row and rethrows anything else', async () => {
    const notFound = Object.assign(new Error('not-found'), {apiError: {error: 'not-found', errorCode: 404}});
    const {client} = makeMock(() => { throw notFound; });
    expect(await new NodeDomainsDataSource(client).getRow('grit', 'issue', 'R1')).toBeNull();
    const {client: c2} = makeMock(() => { throw Object.assign(new Error('boom'), {apiError: {errorCode: 500}}); });
    await expect(new NodeDomainsDataSource(c2).getRow('grit', 'issue', 'R1')).rejects.toThrow('boom');
  });

  it('streams a batch with the mode flags in the query string', async () => {
    const {client, calls} = makeMock(() => ({inserted: 1}));
    const ds = new NodeDomainsDataSource(client);
    await ds.batch('grit', 'issue', Buffer.from('a,b\n1,2'), 'text/csv', {mode: 'upsert', allOrNothing: false});
    expect(calls[0].path).toBe('/domains/grit/issue/batch?mode=upsert&allOrNothing=false');
    expect(calls[0].contentType).toBe('text/csv');
    await ds.queryD42('grit', 'issue', {limit: 1});
    expect(calls[1].body).toEqual({limit: 1, format: 'd42'});
  });
});
