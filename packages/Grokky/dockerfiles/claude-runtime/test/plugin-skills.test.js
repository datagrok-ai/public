// Skills are prose, so nothing type-checks them — a skill can keep instructing the model to call
// a tool that no longer exists and the only symptom is a turn that fails in an odd way. The
// 34-tool -> 5-domain-tool consolidation left three such references behind; these tests are the
// guard against the next rename.

const test = require('node:test');
const assert = require('node:assert');
const fs = require('node:fs');
const path = require('node:path');

const SKILLS = path.join(__dirname, '..', 'plugin', 'skills');

// Retired when the datagrok MCP server moved to one tool per domain (mcp-server/src/ops.ts).
// Matched as whole words so prose like "create a project" is not flagged.
const RETIRED = [
  'list_functions', 'get_function', 'call_function', 'list_scripts', 'list_queries',
  'create_script', 'create_query', 'list_files', 'download_file', 'upload_file',
  'list_projects', 'get_project', 'create_project', 'delete_project', 'search_project',
  'list_recent_projects', 'attach_entity_to_project', 'share_project', 'list_project_shares',
  'list_spaces', 'get_space', 'create_space', 'delete_space', 'create_subspace',
  'list_space_children', 'add_entity_to_space', 'remove_entity_from_space',
  'read_space_file', 'write_space_file', 'delete_space_file',
  'whoami', 'list_connections', 'list_groups', 'list_users',
];

const skillFiles = () => fs.readdirSync(SKILLS)
  .map((d) => path.join(SKILLS, d, 'SKILL.md'))
  .filter((p) => fs.existsSync(p));

test('skills exist and are loadable', () => {
  assert.ok(skillFiles().length >= 10, 'expected the full skill set');
});

test('no skill references a retired MCP tool name', () => {
  const offences = [];
  for (const file of skillFiles()) {
    const text = fs.readFileSync(file, 'utf8');
    for (const name of RETIRED) {
      const re = new RegExp(`\\b${name}\\b`);
      if (re.test(text))
        offences.push(`${path.basename(path.dirname(file))}: ${name}`);
    }
  }
  assert.deepEqual(offences, [],
    'skills name tools that no longer exist — use datagrok_<domain>(op:"…", args:{…})');
});

test('every skill has name and description frontmatter', () => {
  for (const file of skillFiles()) {
    const text = fs.readFileSync(file, 'utf8');
    const fm = text.match(/^---\r?\n([\s\S]*?)\r?\n---/);
    assert.ok(fm, `${file} has no frontmatter`);
    assert.match(fm[1], /^name:\s*\S+/m, `${file} has no name`);
    // The description is what the model routes on; an empty one makes the skill unreachable.
    assert.match(fm[1], /^description:\s*\S.{40,}/m, `${file} needs a routing description`);
  }
});
