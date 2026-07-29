// System prompt assembly. The prefix these produce is re-read on every API call of every turn
// (measured at ~28k tokens/call, see docs/BENCHMARK.md), so anything that silently bloats it is a
// per-turn tax — which is exactly what the frontmatter stripper did: it was written LF-only while
// the skills are checked in CRLF, so it matched nothing and inlined each skill's YAML header.

const test = require('node:test');
const assert = require('node:assert');
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

// Only meaningful when the skills are present (they are inside the image and in a dev checkout);
// skipped rather than failed elsewhere so the suite stays runnable without the plugin.
test('inlined skills carry no YAML frontmatter', (t) => {
  const full = buildSystemPrompt();
  if (!full.includes('## Inlined Skills'))
    return t.skip('plugin skills not mounted at /app/plugin');
  const inlined = full.slice(full.indexOf('## Inlined Skills'));
  assert.ok(!/^name:\s/m.test(inlined), 'frontmatter leaked into the system prompt');
  assert.ok(!/^description:\s/m.test(inlined), 'frontmatter leaked into the system prompt');
  assert.ok(!/^---\s*$/m.test(inlined.replace(/^---$/gm, (m, i) =>
    // the section separator between skills is a legitimate '---'
    inlined.slice(0, i).includes('### ') ? '' : m)), 'stray frontmatter delimiters');
});
