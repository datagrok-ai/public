import {test} from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, readdirSync, readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {U2_RULES} from '../src/core/rules.js';
import {U2_KINDS} from '../src/core/kinds.js';
import {agentsRules, kindNotes, readmeRules, ruleBlock} from '../scripts/gen-manifest.mjs';

const require = createRequire(import.meta.url);
const ts = require('typescript');
const root = fileURLToPath(new URL('..', import.meta.url));
const publicRoot = join(root, '..', '..');
const read = (path) => readFileSync(path, 'utf8').replace(/\r\n/g, '\n');

function walk(dir, keep, out = []) {
  for (const e of readdirSync(dir, {withFileTypes: true})) {
    const path = join(dir, e.name);
    if (e.isDirectory())
      walk(path, keep, out);
    else if (keep(path))
      out.push(path);
  }
  return out;
}

/** The rules a guard may name, read off the things that hold them. */
function known() {
  const plugin = require(join(publicRoot, 'libraries', 'eslint-plugin-u2', 'index.js'));
  const usage = new Set();
  const sources = [...walk(join(publicRoot, 'js-api', 'src', 'u2core'), (p) => /(?<!\.d)\.ts$/.test(p)),
    ...walk(join(root, 'src'), (p) => /(?<!\.d)\.ts$/.test(p))];
  for (const file of sources)
    for (const m of read(file).matchAll(/Usage\.(?:report|fail)\('([\w-]+)'/g))
      usage.add(m[1]);
  const operability = /export type OperabilityRule =([^;]+);/
    .exec(read(join(publicRoot, 'js-api', 'src', 'u2core', 'inspect.ts')))[1];
  const audit = read(join(root, 'scripts', 'audit-consumers.mjs'));
  const neutral = /const NEUTRAL = \[([^\]]*)\]/.exec(audit)[1];
  const counted = Object.values(JSON.parse(read(join(root, 'consumers-baseline.json'))).trees)[0];
  const types = new Map();
  for (const file of walk(join(root, '.test-build', 'src'), (p) => p.endsWith('.d.ts'))) {
    const sf = ts.createSourceFile(file, read(file), ts.ScriptTarget.Latest, false);
    for (const st of sf.statements) {
      if ((ts.isInterfaceDeclaration(st) || ts.isClassDeclaration(st)) && st.name !== undefined)
        types.set(st.name.text, new Set([...types.get(st.name.text) ?? [], ...st.members.map((m) => m.name?.text)]));
    }
  }
  const textOf = (file) => read(join(root, 'tests', file));
  return {
    lint: new Set(Object.keys(plugin.rules).map((r) => `u2/${r}`)),
    usage,
    probe: new Set([...operability.matchAll(/'([\w-]+)'/g)].map((m) => m[1])),
    type: types,
    ratchet: new Set(Object.keys(counted).filter((k) => typeof counted[k] === 'number' && !neutral.includes(`'${k}'`))),
    api: new Set([.../const FIX = \{([\s\S]*?)\n\};/.exec(textOf('api-conventions.test.js'))[1]
      .matchAll(/^ {2}'([\w-]+)':/gm)].map((m) => m[1])),
    choke: new Set([...textOf('choke-points.test.js').matchAll(/^\s+id: '([^']+)'/gm)].map((m) => m[1])),
    debug: new Set(readdirSync(join(root, 'tests')).filter((f) => f.endsWith('.test.js')).map((f) => f.slice(0, -8))),
  };
}

const KNOWN = known();

/** What is wrong with each rule's guards: a kind no check has, or a name its check does not know. */
function guardProblems(rules, k = KNOWN) {
  const problems = [];
  for (const rule of rules) {
    for (const guard of rule.guards) {
      const at = guard.indexOf(':');
      const [kind, id] = [guard.slice(0, at), guard.slice(at + 1)];
      const [type, member] = id.split('.');
      const ok = kind === 'type' ? k.type.get(type)?.has(member) === true : k[kind]?.has(id) === true;
      if (!ok)
        problems.push(`${rule.id}: ${guard} — no such ${kind in k ? kind : 'guard kind'}`);
    }
  }
  return problems;
}

/** What is wrong with the rules' shape. */
function shapeProblems(rules) {
  const problems = [];
  const ids = new Set();
  for (const r of rules) {
    if (!/^[a-z]+(?:-[a-z]+)*$/.test(r.id) || ids.has(r.id))
      problems.push(`${r.id}: an id is unique and kebab-case`);
    ids.add(r.id);
    for (const field of ['rule', 'why', 'detail']) {
      if (typeof r[field] !== 'string' || !/[.!?]$/.test(r[field]) || r[field].includes('\n'))
        problems.push(`${r.id}: ${field} is a sentence on one line`);
    }
    for (const field of ['wrong', 'right']) {
      if (typeof r[field] !== 'string' || r[field] === '' || /[`\n]/.test(r[field]))
        problems.push(`${r.id}: ${field} is one line of code with no backtick`);
    }
    if (r.wrong === r.right)
      problems.push(`${r.id}: the wrong and the right line differ`);
  }
  return problems;
}

/** Rules nothing catches yet; adding a guard takes a rule off this list, and a new rule needs a guard. */
const UNGUARDED = ['on-input', 'reset-follow', 'two-verdicts', 'tabs-mark', 'label-column', 'docked-editor', 'dock',
  'find-from-root', 'dialog-cancel', 'settings-scope', 'dialog-outlives-row', 'balloon-content', 'app-header',
  'confirm', 'records-checklist', 'catalog-table'];

test('rules: an id, one sentence for the rule, its why and its detail, one line each for the wrong and the right code',
  () => {
    assert.ok(U2_RULES.length >= 31, `${U2_RULES.length} rules`);
    assert.deepEqual(shapeProblems(U2_RULES), []);
  });

test('rules: every guard names a lint rule, usage code, probe rule, type member, ratchet count, API convention, ' +
  'choke point or debug test that exists', () => {
  assert.deepEqual(guardProblems(U2_RULES), []);
});

test('rules: a rule with no guard is listed as unguarded, and a listed one has none', () => {
  assert.deepEqual(U2_RULES.filter((r) => r.guards.length === 0).map((r) => r.id), UNGUARDED);
});

test('rules: every u2 lint rule and every usage code the engine and u2 report holds some rule', () => {
  const named = new Set(U2_RULES.flatMap((r) => r.guards));
  assert.deepEqual([...KNOWN.lint].filter((r) => !named.has(`lint:${r}`)), [], 'lint rules no rule names');
  assert.deepEqual([...KNOWN.usage].filter((c) => !named.has(`usage:${c}`)), [], 'usage codes no rule names');
});

test('rules: the kinds a rule speaks of exist, each gets one whenToUse at most, and kinds.ts has none of its own',
  () => {
    const kinds = new Set(U2_KINDS.map((k) => k.kind));
    assert.deepEqual(U2_RULES.flatMap((r) => Object.keys(r.kinds ?? {}).filter((k) => !kinds.has(k))
      .map((k) => `${r.id}: ${k}`)), []);
    assert.doesNotThrow(() => kindNotes(U2_RULES));
    assert.deepEqual(U2_KINDS.filter((k) => k.whenToUse !== undefined || k.notFor !== undefined).map((k) => k.kind),
      [], 'a kind\'s whenToUse and notFor come from src/core/rules.ts');
  });

test('rules: README.md and AGENTS.md carry the generated lists, manifest.json the count (`npm run manifest`)', () => {
  assert.equal(ruleBlock(read(join(root, 'README.md'))), readmeRules(U2_RULES));
  assert.equal(ruleBlock(read(join(root, 'AGENTS.md'))), agentsRules(U2_RULES));
  assert.equal(JSON.parse(read(join(root, 'manifest.json'))).rules.count, U2_RULES.length);
});

test('AGENTS.md: two pages at most, the eight sections, and no word of how the library was made', () => {
  const text = read(join(root, 'AGENTS.md'));
  assert.ok(text.split('\n').length <= 200, `${text.split('\n').length} lines`);
  assert.deepEqual([...text.matchAll(/^## (.+)$/gm)].map((m) => m[1]), ['Entry points', 'Ten ideas', 'Rules',
    'Choosing a component', 'Runnable samples', 'Testing', 'Checks', 'When u2 lacks something']);
  assert.equal([...text.matchAll(/^\d+\. /gm)].length, 10, 'ten ideas');
  assert.deepEqual(text.match(/\b(?:WO-\d+\w*|lanes?|phases?|plans?|rounds?|waves?|coordinator|tonight|the lead)\b/gi),
    null);
});

test('mutation: a guard naming nothing, a guard kind nobody checks, a two-line wrong and a duplicate id fail', () => {
  const [rule] = U2_RULES;
  assert.deepEqual(guardProblems([{...rule, guards: ['lint:u2/no-such-rule', 'usage:no-such-code',
    'type:Action.nothing', 'ratchet:files', 'review:by-eye']}]), [
    `${rule.id}: lint:u2/no-such-rule — no such lint`,
    `${rule.id}: usage:no-such-code — no such usage`,
    `${rule.id}: type:Action.nothing — no such type`,
    `${rule.id}: ratchet:files — no such ratchet`,
    `${rule.id}: review:by-eye — no such guard kind`,
  ]);
  assert.deepEqual(shapeProblems([rule, {...rule, wrong: 'a\nb'}]), [`${rule.id}: an id is unique and kebab-case`,
    `${rule.id}: wrong is one line of code with no backtick`]);
});

test('mutation: a README rules list edited by hand no longer matches', () => {
  const block = readmeRules(U2_RULES);
  assert.notEqual(block.replace('Import from the entry points only.', 'Import from anywhere.'), block);
  assert.equal(ruleBlock(`a\n${block}\nb`), block);
  assert.equal(ruleBlock('no markers here'), null);
});

test('the guards are read off real files: a lint rule, a usage code, a probe rule, a type, a count each', () => {
  for (const [kind, sample] of [['lint', 'u2/import-entry'], ['usage', 'owner-missing'],
    ['probe', 'disabled-without-reason'], ['ratchet', 'deepImports'], ['api', 'no-enabled'], ['choke', 'handlers'],
    ['debug', 'rebuilt-scope-warning']])
    assert.ok(KNOWN[kind].has(sample), `${kind}: ${sample}`);
  assert.ok(KNOWN.type.get('Action')?.has('run'));
  assert.ok(existsSync(join(root, 'tests', 'rebuilt-scope-warning.test.js')));
});
