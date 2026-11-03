const test = require('node:test');
const assert = require('node:assert');
const http = require('node:http');

const TOKEN = 'eyJhbGciOiJSUzI1NiJ9.CONTAINER-TOKEN.sig';
let reply;
let seen;
const datlas = http.createServer((req, res) => {
  seen = {url: req.url, authorization: req.headers.authorization};
  res.writeHead(reply.status, {'content-type': 'application/json', ...reply.headers});
  res.end(JSON.stringify(reply.body ?? {}));
});
let credentials;

test.before(async () => {
  await new Promise((resolve) => datlas.listen(0, '127.0.0.1', resolve));
  Object.assign(process.env, {DATAGROK_API_URL: `http://127.0.0.1:${datlas.address().port}`, apiKey: 'container-key',
    model: 'container-model'});
  credentials = require('../dist/credentials.js');
});

test.after(() => datlas.close());

test('the container token buys the caller\'s provider credentials, and only those', async () => {
  reply = {status: 200, body: {provider: 'Azure', apiKey: 'group-key', model: 'claude-x', region: '', dbPassword: 'p'}};
  const c = await credentials.userCredentials(TOKEN);
  assert.deepEqual(seen, {url: '/docker/credentials', authorization: `Bearer ${TOKEN}`});
  assert.deepEqual(c.values, {provider: 'Azure', apiKey: 'group-key', model: 'claude-x'});
});

test('no token, or none of the caller\'s credentials, leaves the container\'s own', async () => {
  assert.equal(await credentials.userCredentials(undefined), credentials.containerCredentials);
  for (reply of [{status: 404}, {status: 200, headers: {'api-error': 'NotFound', 'api-error-code': '404'}},
    {status: 200, body: {dbPassword: 'p'}}])
    assert.equal(await credentials.userCredentials(TOKEN), credentials.containerCredentials);
});

test('the container\'s key never goes out with a group\'s provider settings', async () => {
  reply = {status: 200, body: {apiKey: 'container-key', model: 'container-model', provider: 'Custom', baseUrl: 'https://elsewhere'}};
  assert.equal(await credentials.userCredentials(TOKEN), credentials.containerCredentials);
});

test('a refused token fails with its status, never the token', async () => {
  reply = {status: 200, headers: {'api-error': 'UserException', 'api-error-code': '401'}, body: {message: TOKEN}};
  await assert.rejects(credentials.userCredentials(TOKEN),
    {message: 'The assistant could not read its credentials: HTTP 401 from /docker/credentials'});
});
