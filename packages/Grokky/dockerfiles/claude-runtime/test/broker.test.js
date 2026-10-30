const test = require('node:test');
const assert = require('node:assert');
const http = require('node:http');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {createBroker, resolveConfig} = require('../dist/broker/broker.js');

function listen(server) {
  return new Promise((resolve) => server.listen(0, '127.0.0.1', () => resolve(server.address().port)));
}
const close = (s) => new Promise((r) => s.close(r));

function stubUpstream() {
  const seen = {};
  const server = http.createServer((req, res) => {
    seen.headers = req.headers;
    seen.url = req.url;
    let body = '';
    req.on('data', (c) => body += c);
    req.on('end', () => {
      seen.body = body;
      res.end('{"ok":true}');
    });
  });
  return {server, seen};
}

const baseCfg = (over) => ({
  mode: 'anthropic', credentialsPath: '/nonexistent', models: {},
  upstreams: {anthropic: 'http://unused'}, ...over,
});

test('/status carries no secret', async () => {
  const broker = createBroker(baseCfg({anthropicKey: 'sk-ant-REALSECRET'}));
  const port = await listen(broker);
  try {
    const status = await (await fetch(`http://127.0.0.1:${port}/status`)).json();
    assert.equal(status.mode, 'anthropic');
    assert.ok(!JSON.stringify(status).includes('REALSECRET'));
  } finally {
    await close(broker);
  }
});

test('anthropic mode injects the real key and strips the inbound placeholder', async () => {
  const up = stubUpstream();
  const upPort = await listen(up.server);
  try {
    const broker = createBroker(baseCfg({
      anthropicKey: 'sk-ant-REALSECRET', upstreams: {anthropic: `http://127.0.0.1:${upPort}`},
    }));
    const port = await listen(broker);
    try {
      await fetch(`http://127.0.0.1:${port}/v1/messages`, {
        method: 'POST', headers: {'x-api-key': 'sk-ant-broker-placeholder', 'authorization': 'Bearer stale'}, body: '{}',
      });
      assert.equal(up.seen.headers['x-api-key'], 'sk-ant-REALSECRET');
      assert.ok(!up.seen.headers['authorization']);
    } finally {
      await close(broker);
    }
  } finally {
    await close(up.server);
  }
});

test('subscription mode injects the OAuth Bearer and strips the placeholder', async () => {
  const credPath = path.join(os.tmpdir(), `creds-${Date.now()}.json`);
  fs.writeFileSync(credPath, JSON.stringify({claudeAiOauth: {accessToken: 'OAUTHTOKEN123'}}));
  const up = stubUpstream();
  const upPort = await listen(up.server);
  try {
    const broker = createBroker(baseCfg({
      mode: 'subscription', credentialsPath: credPath, upstreams: {anthropic: `http://127.0.0.1:${upPort}`},
    }));
    const port = await listen(broker);
    try {
      await fetch(`http://127.0.0.1:${port}/v1/messages`, {
        method: 'POST', headers: {'x-api-key': 'sk-ant-broker-placeholder'}, body: '{}',
      });
      assert.equal(up.seen.headers['authorization'], 'Bearer OAUTHTOKEN123');
      assert.ok(!up.seen.headers['x-api-key']);
    } finally {
      await close(broker);
    }
  } finally {
    await close(up.server);
    fs.unlinkSync(credPath);
  }
});

test('MCP tokenization injects the datagrok key server-side, so it never lands on the CLI argv', async () => {
  const up = stubUpstream();
  const upPort = await listen(up.server);
  try {
    const broker = createBroker(baseCfg({}));
    const port = await listen(broker);
    try {
      const reg = await (await fetch(`http://127.0.0.1:${port}/mcp-session`, {
        method: 'POST',
        headers: {'content-type': 'application/json'},
        body: JSON.stringify({targetUrl: `http://127.0.0.1:${upPort}/mcp`, apiKey: 'DGKEY', apiUrl: 'http://dg'}),
      })).json();
      await fetch(`http://127.0.0.1:${port}/mcp/${reg.token}`, {method: 'POST', body: '{}'});
      assert.equal(up.seen.headers['x-user-api-key'], 'DGKEY');
      assert.equal((await fetch(`http://127.0.0.1:${port}/mcp/nope`, {method: 'POST', body: '{}'})).status, 404);
    } finally {
      await close(broker);
    }
  } finally {
    await close(up.server);
  }
});

test('DATABRICKS_URL selects databricks mode, its Anthropic route, and the endpoint name as every model tier', () => {
  const cfg = resolveConfig({
    DATABRICKS_URL: 'https://dbc-123.cloud.databricks.com/serving-endpoints/claude-sonnet/invocations',
    DATABRICKS_TOKEN: 'dapiSECRET', provider: 'Bedrock', apiKey: 'sk-ant-other',
  });
  assert.equal(cfg.mode, 'databricks');
  assert.equal(cfg.upstreams.databricks, 'https://dbc-123.cloud.databricks.com/serving-endpoints/anthropic');
  assert.deepEqual(cfg.models, {opus: 'claude-sonnet', sonnet: 'claude-sonnet', haiku: 'claude-sonnet'});

  const tierEnvIgnored = resolveConfig({DATABRICKS_URL: 'https://dbc-123.cloud.databricks.com/serving-endpoints/x', opusModel: 'big'});
  assert.deepEqual(tierEnvIgnored.models, {opus: 'x', sonnet: 'x', haiku: 'x'});

  assert.notEqual(resolveConfig({apiKey: 'sk-ant-x'}).mode, 'databricks');
});

test('an AI Gateway DATABRICKS_URL selects the gateway route and DATABRICKS_MODEL as every model tier', () => {
  for (const path of ['anthropic', 'anthropic/v1/messages', 'mlflow/v1/responses']) {
    const cfg = resolveConfig({
      DATABRICKS_URL: `https://dbc-123.cloud.databricks.com/ai-gateway/${path}`,
      DATABRICKS_TOKEN: 'dapiSECRET', DATABRICKS_MODEL: 'system.ai.claude-sonnet-5-5',
    });
    assert.equal(cfg.mode, 'databricks');
    assert.equal(cfg.upstreams.databricks, 'https://dbc-123.cloud.databricks.com/ai-gateway/anthropic');
    assert.deepEqual(cfg.models, {
      opus: 'system.ai.claude-sonnet-5-5', sonnet: 'system.ai.claude-sonnet-5-5', haiku: 'system.ai.claude-sonnet-5-5',
    });
  }
});

test('the AI Gateway gets the token and the coding-agent header, and refuses a call without a model', async () => {
  const up = stubUpstream();
  const upPort = await listen(up.server);
  const upstreams = {anthropic: 'http://unused', databricks: `http://127.0.0.1:${upPort}/ai-gateway/anthropic`};
  const call = (port) => fetch(`http://127.0.0.1:${port}/v1/messages?beta=true`, {
    method: 'POST', headers: {'x-api-key': 'sk-ant-broker-placeholder'}, body: '{}',
  });
  try {
    const broker = createBroker(baseCfg({
      mode: 'databricks', databricksToken: 'dapiSECRET', upstreams, models: {sonnet: 'system.ai.claude-sonnet-5-5'},
    }));
    const noModel = createBroker(baseCfg({mode: 'databricks', databricksToken: 'dapiSECRET', upstreams, models: {}}));
    const port = await listen(broker);
    const noModelPort = await listen(noModel);
    try {
      await call(port);
      assert.equal(up.seen.url, '/ai-gateway/anthropic/v1/messages?beta=true');
      assert.equal(up.seen.headers['authorization'], 'Bearer dapiSECRET');
      assert.equal(up.seen.headers['x-databricks-use-coding-agent-mode'], 'true');
      assert.ok(!up.seen.headers['x-api-key']);

      const refused = await call(noModelPort);
      assert.equal(refused.status, 500);
      assert.match((await refused.json()).error, /DATABRICKS_MODEL/);
    } finally {
      await close(broker);
      await close(noModel);
    }
  } finally {
    await close(up.server);
  }
});

test('databricks mode sends the token as Bearer to the anthropic route and strips the placeholder', async () => {
  const up = stubUpstream();
  const upPort = await listen(up.server);
  try {
    const broker = createBroker(baseCfg({
      mode: 'databricks', databricksToken: 'dapiSECRET',
      upstreams: {anthropic: 'http://unused', databricks: `http://127.0.0.1:${upPort}/serving-endpoints/anthropic`},
    }));
    const port = await listen(broker);
    try {
      await fetch(`http://127.0.0.1:${port}/v1/messages?beta=true`, {
        method: 'POST', headers: {'x-api-key': 'sk-ant-broker-placeholder'}, body: '{}',
      });
      assert.equal(up.seen.url, '/serving-endpoints/anthropic/v1/messages?beta=true');
      assert.equal(up.seen.headers['authorization'], 'Bearer dapiSECRET');
      assert.ok(!up.seen.headers['x-api-key']);
    } finally {
      await close(broker);
    }
  } finally {
    await close(up.server);
  }
});
