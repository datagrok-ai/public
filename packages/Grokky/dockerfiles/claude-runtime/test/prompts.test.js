// System prompt assembly. The prefix these produce is re-read on every API call of every turn, so anything that
// silently bloats it is a per-turn tax.

const path = require('node:path');
const test = require('node:test');
const assert = require('node:assert');

process.env.GROKKY_PLUGIN_DIR = path.join(__dirname, '..', 'plugin');
const {buildSystemPrompt} = require('../dist/prompts.js');

test('bash and none modes stay minimal', () => {
  assert.match(buildSystemPrompt('bash'), /^Execute the given shell command/);
  assert.equal(buildSystemPrompt('none'), '');
});

test('the full prompt is the default for any other mode', () => {
  const full = buildSystemPrompt();
  assert.ok(full.length > 1000);
  assert.equal(buildSystemPrompt('datagrok'), full);
  assert.equal(buildSystemPrompt(undefined), full);
});

test('the full prompt documents the domain MCP tools', () => {
  const full = buildSystemPrompt();
  for (const tool of ['datagrok_functions', 'datagrok_files', 'datagrok_projects',
    'datagrok_spaces', 'datagrok_platform'])
    assert.ok(full.includes(tool), `system prompt never mentions ${tool}`);
});

test('the full prompt prefers widget quick functions for viewer-targeted requests', () => {
  const full = buildSystemPrompt();
  assert.ok(full.includes('list_view_widgets'), 'widget tool never mentioned');
  assert.match(full, /prefer\s+its own quick functions/i, 'no widget-function preference rule');
  for (const fn of ['zoomIn', 'resetView', 'colorBy', 'splitBy'])
    assert.ok(full.includes(fn), `viewer quick command ${fn} not named in the prompt`);
  // Function names are camelCase across the board; snake_case in the widget list is a regression.
  for (const stale of ['get_dialog_info', 'list_panes', 'set_column', 'list_items'])
    assert.ok(!full.includes(stale), `stale snake_case name ${stale} in the prompt`);
  assert.match(full, /preference, not a hard rule/i, 'the rule must stay soft');
});

test('inlined skills carry no YAML frontmatter', () => {
  const full = buildSystemPrompt();
  assert.ok(full.includes('## Inlined Skills'), 'the inlined skills are missing');
  const inlined = full.slice(full.indexOf('## Inlined Skills'));
  assert.ok(!/^name:\s/m.test(inlined), 'frontmatter leaked into the system prompt');
  assert.ok(!/^description:\s/m.test(inlined), 'frontmatter leaked into the system prompt');
  assert.ok(!/^---\s*$/m.test(inlined.replace(/^---$/gm, (m, i) =>
    // the section separator between skills is a legitimate '---'
    inlined.slice(0, i).includes('### ') ? '' : m)), 'stray frontmatter delimiters');
});

test('the full prompt warns that search_files lists hidden directories and matches names only', () => {
  const full = buildSystemPrompt('datagrok');
  assert.match(full, /`search_files` with `target: "files"` matches file names only/);
  assert.match(full, /hidden directories such as `\.git`/);
});
