// GroundingGate — blocks a Stop once if the turn answered platform questions from memory
// rather than from a source opened this turn.

const test = require('node:test');
const assert = require('node:assert');
const {GroundingGate, isSmallTalk, makesPlatformClaims} = require('../dist/grounding.js');

const blocked = (feedback) => feedback !== undefined;

test('reading a help page grounds the turn', () => {
  const g = new GroundingGate();
  g.afterTool('read_file', {path: 'workspace/help/visualize/viewers/scatter-plot.md'}, '');
  assert.ok(!blocked(g.atStop()));
});

test('reading js-api source grounds the turn', () => {
  const g = new GroundingGate();
  g.afterTool('read_file', {path: 'workspace/js-api/src/dataframe.ts'}, '');
  assert.ok(!blocked(g.atStop()));
});

test('web lookups count as grounding', () => {
  for (const tool of ['web_extract', 'web_search']) {
    const g = new GroundingGate();
    g.afterTool(tool, {url: 'https://datagrok.ai/help'}, '');
    assert.ok(!blocked(g.atStop()), `${tool} should ground`);
  }
});

test('taking an action grounds the turn — it is not a from-memory answer', () => {
  const g = new GroundingGate();
  g.afterTool('datagrok_exec', {code: 'return 1'}, '');
  assert.ok(!blocked(g.atStop()));
});

test('reading an unrelated file does NOT ground the turn', () => {
  const g = new GroundingGate();
  g.afterTool('read_file', {path: 'workspace/packages/Chem/package.json'}, '');
  assert.ok(blocked(g.atStop()));
});

test('an answer with no tool use at all is blocked once', () => {
  const g = new GroundingGate();
  const r = g.atStop();
  assert.ok(blocked(r));
  assert.match(r, /INDEX\.md/, 'the block must point at the index, not at grepping');
  assert.match(r, /NO_REVISION/, 'and must offer the keep-original escape hatch');
});

test('the gate blocks at most once per turn', () => {
  const g = new GroundingGate();
  assert.ok(blocked(g.atStop()));
  assert.ok(!blocked(g.atStop()), 'a second Stop must pass or the turn cannot end');
});

test('block reason is framed as internal feedback', () => {
  const g = new GroundingGate();
  const r = g.atStop();
  assert.match(r, /never mention, quote, or allude to it/);
});

test('onBlock fires so the runtime can start a hidden revision', () => {
  let calls = 0;
  const g = new GroundingGate(() => calls++);
  g.atStop();
  g.atStop();
  assert.equal(calls, 1, 'only the single real block should trigger a revision');
});

test('summary reports both flags', () => {
  const g = new GroundingGate();
  assert.equal(g.summary(), 'grounded=false blocked=false');
  g.afterTool('read_file', {path: 'workspace/help/index.md'}, '');
  assert.equal(g.summary(), 'grounded=true blocked=false');
});

test('search_files over help counts, over other paths does not', () => {
  const yes = new GroundingGate();
  yes.afterTool('search_files', {pattern: 'scatter', path: 'workspace/help/'}, '');
  assert.ok(!blocked(yes.atStop()));

  const no = new GroundingGate();
  no.afterTool('search_files', {pattern: '*.ts', target: 'files', path: 'workspace/packages/'}, '');
  assert.ok(blocked(no.atStop()));
});

// isSmallTalk decides whether the gate is armed AT ALL (session.ts), so both failure directions
// are real: a greeting that matches nothing costs a hidden revision call and a "Revising…" flash;
// a question that matches by accident answers ungrounded with no safety net.
test('isSmallTalk accepts pure greetings and acknowledgements', () => {
  for (const m of ['hello', 'Hi!', 'hey', 'HELLO', 'thanks', 'Thank you!', 'ok', 'cool',
    'good morning', 'got it', 'hello there', 'hi grokky', 'thanks again', 'bye', 'how are you?'])
    assert.ok(isSmallTalk(m), `"${m}" should be small talk`);
});

test('isSmallTalk rejects anything with actual content', () => {
  for (const m of ['hello, how do I add a scatter plot?', 'hi — can you filter this table?',
    'what is a project?', 'help', 'select all rows', 'thanks, now add a histogram',
    'ok now delete the AGE column', 'hello world program in python', ''])
    assert.ok(!isSmallTalk(m), `"${m}" must NOT be small talk`);
});

test('isSmallTalk sees through the prepended workspace context', () => {
  // The panel sends '<workspace context>\n---\n\nhello' when a view is open (panel.ts
  // prependViewContext) — the dashboard-open case that made the regression visible. Stripping is
  // anchored on the context block's fixed header, not on '---' alone.
  const ctx = 'Workspace state (live, changes as the user navigates):\n' +
    'Current view: "Dashboard" (TableView)\nTable: demog, 5850 rows';
  assert.ok(isSmallTalk(ctx + '\n---\n\nhello'), 'context + greeting is still small talk');
  assert.ok(!isSmallTalk(ctx + '\n---\n\nhow many rows are selected?'),
    'context + real question keeps the gate');
  // A user message that merely contains a --- separator is never split — fail-closed.
  assert.ok(!isSmallTalk('explain this:\n---\n\nok'), 'user-typed --- must not discard the head');
});

test('loading a skill counts as grounding', () => {
  // The system prompt names skills as the source of truth for code/API answers; blocking a
  // skill-grounded turn costs a revision call that always ends NO_REVISION.
  const g = new GroundingGate();
  g.afterTool('skill_view', {name: 'datagrok-viewers'}, '');
  assert.ok(!blocked(g.atStop()));
});

// Content-aware gating (from the 36-turn A/B: every Stop-block ended NO_REVISION, so the block
// is only worth its revision call when the answer actually makes platform how-to claims).

test('makesPlatformClaims: UI-instruction answers are claims', () => {
  for (const a of [
    'Click Tools | Data Science | Cluster, then pick the columns.',
    'Right-click the column header and choose Rename.',
    'Open the Add viewer dialog from the toolbar.',
    'Go to File > Open and select your CSV.',
    'Use the context panel to edit properties.',
    'Press Ctrl+F to open the search.',
  ])
    assert.ok(makesPlatformClaims(a), `should be a claim: "${a}"`);
});

test('makesPlatformClaims: data answers are not claims', () => {
  for (const a of [
    'The **demog** table has **5,850 rows** and **11 columns**: USUBJID, AGE, SEX, RACE.',
    'The average WEIGHT across all patients is 78.4 kg.',
    'There are 878 patients older than 65.',
    'Correlation between HEIGHT and WEIGHT is 0.71 — moderately strong.',
  ])
    assert.ok(!makesPlatformClaims(a), `should NOT be a claim: "${a}"`);
});

test('an ungrounded answer WITHOUT platform claims is not blocked', () => {
  const g = new GroundingGate(undefined, () => 'The table has 5,850 rows and 11 columns.');
  assert.ok(!blocked(g.atStop()));
});

test('an ungrounded answer WITH platform claims is blocked once', () => {
  const g = new GroundingGate(undefined, () => 'Click Tools | X | Y to do that.');
  assert.ok(blocked(g.atStop()));
  assert.ok(!blocked(g.atStop()), 'still at most one block per turn');
});

test('a grounded claimy answer is not blocked', () => {
  const g = new GroundingGate(undefined, () => 'Click Add | Scatter plot (from the docs).');
  g.afterTool('read_file', {path: 'workspace/help/visualize/viewers/scatter-plot.md'}, '');
  assert.ok(!blocked(g.atStop()));
});

test('without an answerText getter the gate keeps the old always-block behavior', () => {
  const g = new GroundingGate();
  assert.ok(blocked(g.atStop()), 'fail-closed when no text is available');
});

test('a terminal search of the docs or the js-api source grounds the turn, any other command does not', () => {
  for (const [args, grounds] of [[{command: 'grep -rn "scatter" workspace/help/visualize'}, true],
    [{command: 'rg addViewer workspace/js-api/src'}, true], [{command: 'ls workspace/packages'}, false]]) {
    const g = new GroundingGate(undefined, () => 'Click Add | Scatter plot.');
    g.afterTool('terminal', args, '');
    assert.equal(!blocked(g.atStop()), grounds, args.command);
  }
});

test('state-changing browser, view and domain calls ground the turn; reads and questions do not', () => {
  for (const [name, args, grounds] of [['datagrok_exec', {code: 'x'}, true],
    ['call_view_function', {name: 'setColumn'}, true],
    ['mcp__datagrok__datagrok_functions', {op: 'call', args: {name: 'X:y'}}, true],
    ['list_view_functions', {query: 'x'}, false], ['get_view_function_result', {name: 'getRange'}, false],
    ['mcp__datagrok__datagrok_functions', {op: 'list'}, false], ['datagrok_show_entities', {entities: []}, false],
    ['AskUserQuestion', {questions: []}, false]]) {
    const g = new GroundingGate(undefined, () => 'Click the toolbar icon.');
    g.afterTool(name, args, '{}');
    assert.equal(!blocked(g.atStop()), grounds, name);
  }
});
