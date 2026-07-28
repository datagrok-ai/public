// GroundingGate — blocks a Stop once if the turn answered platform questions from memory
// rather than from a source opened this turn.
// Note: the hook callbacks are async; every invocation must be awaited.

const test = require('node:test');
const assert = require('node:assert');
const {GroundingGate} = require('../dist/grounding.js');

const post = (toolName, input) => ({hook_event_name: 'PostToolUse', tool_name: toolName, tool_input: input});
const stop = () => ({hook_event_name: 'Stop'});
const blocked = (r) => r.decision === 'block';

test('reading a help page grounds the turn', async () => {
  const g = new GroundingGate();
  await g.postToolUse(post('Read', {file_path: 'workspace/help/visualize/viewers/scatter-plot.md'}));
  assert.ok(!blocked(await g.stop(stop())));
});

test('reading js-api source grounds the turn', async () => {
  const g = new GroundingGate();
  await g.postToolUse(post('Read', {file_path: 'workspace/js-api/src/dataframe.ts'}));
  assert.ok(!blocked(await g.stop(stop())));
});

test('web lookups count as grounding', async () => {
  for (const tool of ['WebFetch', 'WebSearch']) {
    const g = new GroundingGate();
    await g.postToolUse(post(tool, {url: 'https://datagrok.ai/help'}));
    assert.ok(!blocked(await g.stop(stop())), `${tool} should ground`);
  }
});

test('taking an action grounds the turn — it is not a from-memory answer', async () => {
  const g = new GroundingGate();
  await g.postToolUse(post('mcp__datagrok-browser__datagrok_exec', {code: 'return 1'}));
  assert.ok(!blocked(await g.stop(stop())));
});

test('reading an unrelated file does NOT ground the turn', async () => {
  const g = new GroundingGate();
  await g.postToolUse(post('Read', {file_path: 'workspace/packages/Chem/package.json'}));
  assert.ok(blocked(await g.stop(stop())));
});

test('an answer with no tool use at all is blocked once', async () => {
  const g = new GroundingGate();
  const r = await g.stop(stop());
  assert.ok(blocked(r));
  assert.match(r.reason, /INDEX\.md/, 'the block must point at the index, not at grepping');
  assert.match(r.reason, /NO_REVISION/, 'and must offer the keep-original escape hatch');
});

test('the gate blocks at most once per turn', async () => {
  const g = new GroundingGate();
  assert.ok(blocked(await g.stop(stop())));
  assert.ok(!blocked(await g.stop(stop())), 'a second Stop must pass or the turn cannot end');
});

test('block reason is framed as internal feedback', async () => {
  const g = new GroundingGate();
  const r = await g.stop(stop());
  assert.match(r.reason, /never mention, quote, or allude to it/);
});

test('onBlock fires so the runtime can start a hidden revision', async () => {
  let calls = 0;
  const g = new GroundingGate(() => calls++);
  await g.stop(stop());
  await g.stop(stop());
  assert.equal(calls, 1, 'only the single real block should trigger a revision');
});

test('summary reports both flags', async () => {
  const g = new GroundingGate();
  assert.equal(g.summary(), 'grounded=false blocked=false');
  await g.postToolUse(post('Read', {file_path: 'workspace/help/index.md'}));
  assert.equal(g.summary(), 'grounded=true blocked=false');
});

test('grep and glob over help count, over other paths do not', async () => {
  const yes = new GroundingGate();
  await yes.postToolUse(post('Grep', {pattern: 'scatter', path: 'workspace/help/'}));
  assert.ok(!blocked(await yes.stop(stop())));

  const no = new GroundingGate();
  await no.postToolUse(post('Glob', {pattern: '**/*.ts', path: 'workspace/packages/'}));
  assert.ok(blocked(await no.stop(stop())));
});
