const test = require('node:test');
const assert = require('node:assert');
const {spawnSync} = require('node:child_process');
const {buildCliEnv, PLACEHOLDER_KEY, BROKER_BASE} = require('../dist/broker/provider-env.js');
const {buildOptions, createBrowserExecServer} = require('../dist/query-options.js');

const SENTINELS = {
  apiKey: 'LEAK_SENTINEL_ANTHROPIC_KEY',
  awsSecretAccessKey: 'LEAK_SENTINEL_AWS_SECRET',
  awsBearerToken: 'LEAK_SENTINEL_AWS_BEARER',
  foundryApiKey: 'LEAK_SENTINEL_FOUNDRY_KEY',
};
const ALL_SENTINELS = Object.values(SENTINELS);
const found = (text) => ALL_SENTINELS.filter((s) => text.includes(s));

function runBash(cmd, env) {
  const r = spawnSync('bash', ['-c', cmd], {env, encoding: 'utf8'});
  return (r.stdout || '') + (r.stderr || '');
}

const legacyEnv = () => ({...process.env, ...SENTINELS, ANTHROPIC_API_KEY: SENTINELS.apiKey});

const MODES = ['anthropic', 'subscription', 'bedrock', 'foundry'];
const infoFor = (mode) => ({mode, region: 'us-east-1', foundryResource: 'my-res'});

test('LEAK (pre-fix): `env` exposes credentials to bash — proves the harness detects a leak', () => {
  assert.ok(found(runBash('env', legacyEnv())).length > 0);
});

test('FIXED: the scrubbed CLI env hides credentials from bash', () => {
  for (const mode of MODES)
    assert.deepEqual(found(runBash('env; printenv', buildCliEnv(infoFor(mode)))), [], `mode ${mode} leaked`);
});

test('buildCliEnv exposes only the allowlist — no secret values, no raw provider fields', () => {
  for (const mode of MODES) {
    const env = buildCliEnv(infoFor(mode));
    for (const v of Object.values(env))
      assert.ok(!ALL_SENTINELS.some((s) => String(v).includes(s)));
    for (const forbidden of ['apiKey', 'awsSecretAccessKey', 'awsBearerToken', 'foundryApiKey', 'provider'])
      assert.ok(!(forbidden in env), `${mode}: raw field ${forbidden} must not be in the CLI env`);
  }
  const anth = buildCliEnv(infoFor('anthropic'));
  assert.equal(anth.ANTHROPIC_API_KEY, PLACEHOLDER_KEY);
  assert.equal(anth.ANTHROPIC_BASE_URL, BROKER_BASE);
});

test('buildOptions puts no auth headers on the datagrok MCP config (CLI-argv leak)', () => {
  const brokerUrl = `${BROKER_BASE}/mcp/opaque-token`;
  const opts = buildOptions(createBrowserExecServer(async () => ({})), {mode: 'anthropic'}, undefined, brokerUrl);
  assert.equal(opts.mcpServers.datagrok.url, brokerUrl);
  assert.equal(opts.mcpServers.datagrok.headers, undefined);
});
