const test = require('node:test');
const assert = require('node:assert');
const http = require('node:http');
const {userDirs, UserDirectory} = require('../dist/user/user-dir.js');

const GOOD_ID = 'a1b2c3d4-0000-1111-2222-333344445555';
const users = new Map([
  ['good-token', {id: GOOD_ID}],
  ['evil-token', {id: '../workspace'}],
  ['pin-token', {id: GOOD_ID}],
  ['dedup-token', {id: GOOD_ID}],
]);
let requests = 0;
const server = http.createServer((req, res) => {
  requests++;
  const user = req.url === '/users/current' ? users.get(req.headers['authorization']) : undefined;
  if (!user) {
    res.statusCode = 401;
    return res.end('Unauthorized');
  }
  res.setHeader('content-type', 'application/json');
  res.end(JSON.stringify(user));
});

let apiUrl;
test.before(() => new Promise((r) => server.listen(0, '127.0.0.1', () => {
  server.unref();
  apiUrl = `http://127.0.0.1:${server.address().port}`;
  r();
})));

test('id comes from the server, is cached, and is re-verified after the TTL', async () => {
  assert.equal(await userDirs.resolveUserId('good-token', apiUrl), GOOD_ID);
  const before = requests;
  await userDirs.resolveUserId('good-token');
  assert.equal(requests, before);
  userDirs.verified.get('good-token').at = 0;
  assert.equal(await userDirs.resolveUserId('good-token'), GOOD_ID);
  assert.equal(requests, before + 1);
});

test('a rejected token (forged JWT payloads included) resolves no identity and is not retried within the rejection window', async () => {
  const forged = 'x.' + Buffer.from(JSON.stringify({sub: GOOD_ID})).toString('base64url') + '.y';
  await assert.rejects(() => userDirs.resolveUserId(forged), /HTTP 401/);
  const before = requests;
  await assert.rejects(() => userDirs.resolveUserId(forged), /HTTP 401/);
  assert.equal(requests, before);
});

test('a path-unsafe id from the server is rejected, not joined into /users', async () => {
  await assert.rejects(() => userDirs.resolveUserId('evil-token'), /invalid id/);
});

test('later client input cannot repoint verification at another server', async () => {
  assert.equal(await userDirs.resolveUserId('pin-token', 'http://127.0.0.1:1'), GOOD_ID);
});

test('concurrent calls with an unseen token share one verification request', async () => {
  const before = requests;
  const ids = await Promise.all([userDirs.resolveUserId('dedup-token'), userDirs.resolveUserId('dedup-token')]);
  assert.deepEqual(ids, [GOOD_ID, GOOD_ID]);
  assert.equal(requests, before + 1);
});

test('DATAGROK_API_URL wins over client-supplied URLs; with neither, no identity resolves', async () => {
  process.env.DATAGROK_API_URL = apiUrl;
  const dirs = new UserDirectory();
  delete process.env.DATAGROK_API_URL;
  assert.equal(await dirs.resolveUserId('good-token', 'http://127.0.0.1:1'), GOOD_ID);
  await assert.rejects(() => new UserDirectory().resolveUserId('good-token'), /no Datagrok API URL/);
});
