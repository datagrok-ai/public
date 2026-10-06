// request() must send string bodies (file uploads) over the wire untouched while still
// JSON-encoding object bodies — a stringified upload corrupts every stored file.

import test from 'node:test';
import assert from 'node:assert';
import http from 'node:http';
import {runWithContext, uploadFile, downloadFile, request} from '../dist/shared-api-client.js';

const files = new Map();
const server = http.createServer((req, res) => {
  let body = '';
  req.on('data', (c) => body += c);
  req.on('end', () => {
    if (req.method === 'POST')
      files.set(req.url, body);
    res.setHeader('content-type', 'application/octet-stream');
    res.end(files.get(req.url) ?? '');
  });
});
await new Promise((r) => server.listen(0, '127.0.0.1', r));
server.unref();
const ctx = {apiKey: 'key', apiUrl: `http://127.0.0.1:${server.address().port}`};

test('file upload → download round-trips byte-identically', () => runWithContext(ctx, async () => {
  const content = 'line1\nline2 "quoted" back\\slash\ttab юникод ✓\n{"json": [1, 2, 3]}\n';
  await uploadFile('System:DemoFiles', 'roundtrip.txt', content);
  assert.equal(files.get('/public/v1/files/System.DemoFiles/roundtrip.txt'), content,
    'uploaded body must reach the server as-is, not JSON-encoded');
  assert.equal(await downloadFile('System:DemoFiles', 'roundtrip.txt'), content);
}));

test('object bodies are still JSON-encoded', () => runWithContext(ctx, async () => {
  await request('POST', '/public/v1/functions/Sin/call', {x: 0});
  assert.equal(files.get('/public/v1/functions/Sin/call'), '{"x":0}');
}));
