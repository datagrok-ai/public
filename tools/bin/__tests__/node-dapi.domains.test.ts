import {describe, it, expect} from 'vitest';
import {NodeDomainsDataSource} from '../utils/node-dapi';

interface Call {method: string; path: string; body?: any}

function makeMock() {
  const calls: Call[] = [];
  const client: any = {
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      return {};
    },
    get(path: string) { return this.request('GET', path); },
    post(path: string, body?: any) { return this.request('POST', path, body); },
    del(path: string) { return this.request('DELETE', path); },
  };
  return {client, calls};
}

/** A literal copy of the Dart `DomainRowKey.reserved` list: the test checks that node-dapi.ts escapes exactly
 * these (the Datlas suite domain_row_key_test pins that list to the router's route table). */
const RESERVED = ['access', 'aggregate', 'audit', 'batch', 'delete', 'facets', 'filters', 'query', 'update',
  'version', 'watch', 'new', '.', '..'];

describe('NodeDomainsDataSource row paths', () => {
  it('escapes a key equal to a route word, a dot segment or led by ~ with a ~, then percent-encodes once', async () => {
    const {client, calls} = makeMock();
    const domains = new NodeDomainsDataSource(client);
    for (const key of RESERVED) {
      calls.length = 0;
      await domains.getRow('s', 't', key);
      expect(calls[0].path).toBe(`/domains/s/t/${encodeURIComponent(`~${key}`)}`);
    }
    await domains.getRow('s', 't', '~x');
    expect(calls.at(-1)!.path).toBe('/domains/s/t/~~x');
    await domains.getRow('s', 't', 'a%2Fb');
    expect(calls.at(-1)!.path).toBe('/domains/s/t/a%252Fb');
    await domains.getRow('s', 't', '1,a%2Cb');
    expect(calls.at(-1)!.path).toBe('/domains/s/t/1%2Ca%252Cb');
    await domains.getRow('s', 't', '0f1e2d3c-4b5a-6978-8a9b-0c1d2e3f4a5b');
    expect(calls.at(-1)!.path).toBe('/domains/s/t/0f1e2d3c-4b5a-6978-8a9b-0c1d2e3f4a5b');
  });

  it('every row route takes the same segment; the table routes stay unescaped', async () => {
    const {client, calls} = makeMock();
    const domains = new NodeDomainsDataSource(client);
    await domains.update('s', 't', 'version', {n: 1});
    await domains.deleteRow('s', 't', 'version');
    await domains.rowAudit('s', 't', 'version');
    await domains.access('s', 't');
    await domains.tableAudit('s', 't');
    expect(calls.map((c) => `${c.method} ${c.path}`)).toEqual([
      'PATCH /domains/s/t/~version',
      'DELETE /domains/s/t/~version',
      'GET /domains/s/t/~version/audit',
      'GET /domains/s/t/access',
      'GET /domains/s/t/audit',
    ]);
    expect(calls[0].body).toEqual({values: {n: 1}});
  });
});
