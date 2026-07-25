"use strict";

var _vitest = require("vitest");
var _nodeDapi = require("../utils/node-dapi");
function makeMock(responder) {
  const calls = [];
  const client = {
    async request(method, path, body) {
      calls.push({
        method,
        path,
        body
      });
      return responder(method, path, body);
    },
    get(path) {
      return this.request('GET', path);
    },
    post(path, body) {
      return this.request('POST', path, body);
    },
    del(path) {
      return this.request('DELETE', path);
    }
  };
  return {
    client,
    calls
  };
}
const CONN = {
  name: 'TestPg',
  dataSource: 'PostgreSQL',
  server: 'localhost',
  port: 5432,
  db: 'datagrok'
};
(0, _vitest.describe)('NodeConnectionsDataSource.save', () => {
  (0, _vitest.it)('POSTs to /public/v1/connections without saveCredentials by default', async () => {
    const {
      client,
      calls
    } = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/connections') return {
        ...body,
        id: 'new-id'
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    const saved = await ds.save(CONN);
    (0, _vitest.expect)(saved.id).toBe('new-id');
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/connections');
    (0, _vitest.expect)(calls[0].body).toEqual(CONN);
  });
  (0, _vitest.it)('appends saveCredentials=true when requested', async () => {
    const {
      client,
      calls
    } = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/connections?saveCredentials=true') return {
        ...body,
        id: 'new-id'
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await ds.save(CONN, true);
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/connections?saveCredentials=true');
  });
  (0, _vitest.it)('omits saveCredentials query param when false', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'POST' && path === '/public/v1/connections') return {
        id: 'new-id'
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await ds.save(CONN, false);
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/connections');
  });
});
(0, _vitest.describe)('NodeConnectionsDataSource.test', () => {
  (0, _vitest.it)('POSTs the connection body to /public/v1/connections/test and resolves on "ok"', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'POST' && path === '/public/v1/connections/test') return 'ok';
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await ds.test(CONN);
    (0, _vitest.expect)(calls[0].body).toEqual(CONN);
  });
  (0, _vitest.it)('strips surrounding quotes on a JSON-encoded "ok"', async () => {
    const {
      client
    } = makeMock(() => '"ok"');
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await (0, _vitest.expect)(ds.test(CONN)).resolves.toBeUndefined();
  });
  (0, _vitest.it)('throws the server error text when the response is not ok', async () => {
    const {
      client
    } = makeMock(() => 'connection refused: tcp://localhost:5432');
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await (0, _vitest.expect)(ds.test(CONN)).rejects.toThrow(/connection refused/);
  });
  (0, _vitest.it)('throws a generic message on empty response', async () => {
    const {
      client
    } = makeMock(() => '');
    const ds = new _nodeDapi.NodeConnectionsDataSource(client);
    await (0, _vitest.expect)(ds.test(CONN)).rejects.toThrow(/failed/i);
  });
});