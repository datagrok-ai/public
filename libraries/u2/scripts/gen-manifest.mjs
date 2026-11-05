/* Writes manifest.json from the registry — the machine-readable component list an LLM reads
   before emitting a dg-ui/1 spec — and the kind catalog (src/core/kinds.ts) the bdd library generates
   its vocabulary from; and, from the usage rules (src/core/rules.ts), the kinds' `whenToUse` / `notFor` and
   the rules lists of README.md and AGENTS.md. Run `npm run build` first: this loads the compiled sources
   (`--from .test-build` loads what `npm test` compiled instead).
   With --check it compares instead of writing, so a stale checked-in file fails the build. */
import {readFileSync, realpathSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '../tests/dom-shim.js';

/** The catalog states a tag reaches through a prop it declares: state → [prop, the props that set it]. */
const STATES = {
  disabled: ['disabled', {disabled: 'Fill in the fields above first'}],
  readOnly: ['readOnly', {readOnly: true}],
  invalid: ['required', {required: true}],
  busy: ['busy', {busy: true}],
  empty: ['items', {items: []}],
};

/** States a tag's fresh value keeps it out of: a checkbox holds `false`, a color input black and a
 * font input the default font, so `required` alone leaves them valid. */
const UNREACHED = {invalid: ['u2-bool-input', 'u2-color-input', 'u2-font-input']};

export const RULES_BEGIN = '<!-- rules: generated from src/core/rules.ts by `npm run manifest`; edit them there -->';
export const RULES_END = '<!-- /rules -->';
const WIDTH = 117;

const GUARD_NAMES = {lint: 'lint', usage: 'usage error', probe: 'operability probe', type: 'type', ratchet: 'ratchet',
  api: 'API convention', choke: 'choke point', debug: 'debug warning'};

/** Per kind, what the rules say about it: one `whenToUse`, any number of `notFor`. */
export function kindNotes(rules) {
  const notes = {};
  for (const rule of rules) {
    for (const [kind, note] of Object.entries(rule.kinds ?? {})) {
      const into = notes[kind] ??= {};
      if (note.whenToUse !== undefined) {
        if (into.whenToUse !== undefined)
          throw new Error(`rules: ${kind} gets a second whenToUse from ${rule.id}`);
        into.whenToUse = note.whenToUse;
      }
      if (note.notFor !== undefined)
        (into.notFor ??= []).push(note.notFor);
    }
  }
  return notes;
}

/** "lint `u2/x`; usage error `a`, `b`": the guards by kind, in the order the rule names them. */
export function guardText(guards) {
  const byKind = new Map();
  for (const guard of guards) {
    const at = guard.indexOf(':');
    byKind.set(guard.slice(0, at), [...byKind.get(guard.slice(0, at)) ?? [], guard.slice(at + 1)]);
  }
  return [...byKind].map(([kind, ids]) => kind === 'debug' ?
    `${GUARD_NAMES[kind]} (${ids.map((id) => `\`tests/${id}.test.js\``).join(', ')})` :
    `${GUARD_NAMES[kind]} ${ids.map((id) => `\`${id}\``).join(', ')}`).join('; ');
}

/** Words and whole code spans, greedily to {@link WIDTH}: a code span is never split, so the docs checks read it. */
function wrap(text, first, next) {
  const lines = [];
  let line = first;
  let empty = true;
  for (const word of text.match(/(?:`[^`]*`|[^\s`])+/g) ?? []) {
    if (!empty && line.length + 1 + word.length > WIDTH) {
      lines.push(line);
      line = next + word;
    }
    else
      line += (empty ? '' : ' ') + word;
    empty = false;
  }
  lines.push(line);
  return lines;
}

const held = (rule) => rule.guards.length > 0 ? ` Held by: ${guardText(rule.guards)}.` : '';

/** README.md's "Rules that matter": each rule, how, and what holds it. */
export function readmeRules(rules) {
  return [RULES_BEGIN, ...rules.flatMap((r) => wrap(`**${r.rule}** ${r.detail}${held(r)}`, '- ', '  ')), RULES_END]
    .join('\n');
}

/** AGENTS.md's rules: each rule, why, and the right and the wrong line. */
export function agentsRules(rules) {
  return [RULES_BEGIN, ...rules.flatMap((r) =>
    wrap(`**${r.rule}** ${r.why} Write \`${r.right}\`, not \`${r.wrong}\`.`, '- ', '  ')), RULES_END].join('\n');
}

/** The text between the rules markers, markers included; null when a marker is missing. */
export function ruleBlock(text) {
  const lf = text.replace(/\r\n/g, '\n');
  const from = lf.indexOf(RULES_BEGIN);
  const to = lf.indexOf(RULES_END, from);
  return from < 0 || to < 0 ? null : lf.slice(from, to + RULES_END.length);
}

/** `text` with its rules block replaced by `block`, in the line endings it has. */
export function withRuleBlock(text, block) {
  const current = ruleBlock(text);
  if (current === null)
    throw new Error(`no rules markers (${RULES_BEGIN} … ${RULES_END})`);
  const replaced = text.replace(/\r\n/g, '\n').replace(current, () => block);
  return text.includes('\r\n') ? replaced.replace(/\n/g, '\r\n') : replaced;
}

/** The manifest the checked-in file must equal. `base` is where the compiled sources are: `dist/src/`
 * for the script, the in-place emit (`src/`) for the drift test. */
export async function buildManifest(base = '../dist/src/') {
  const load = (path) => import(new URL(base + path, import.meta.url).href);
  const [{registry}, {registerAll}, {SpecContext, renderSpec}, {U2_KINDS, kindCatalog}, {Scope}, {U2_RULES}] =
    await Promise.all([load('spec/registry.js'), load('spec/registrations.js'), load('spec/spec.js'),
      load('core/kinds.js'), load('core/scope.js'), load('core/rules.js')]);
  registerAll();
  const manifest = registry.manifest();

  // a tag's kind is the data-u2 its built example root carries
  const kindOfTag = new Map();
  const examples = Scope.root('manifest examples');
  for (const meta of registry.metas()) {
    if (meta.visual === false)
      continue;
    const instance = examples.build(() => renderSpec({$schema: 'dg-ui/1', root: meta.example}, new SpecContext(),
      registry));
    const kind = instance.root.children[0]?.dataset.u2;
    if (kind !== undefined)
      kindOfTag.set(meta.tag, kind);
    instance.dispose();
  }
  examples.dispose();
  const tagOf = (kind) => [...kindOfTag].find(([, k]) => k === kind)?.[0];

  const components = manifest.components.map((c) => {
    const names = c.props.map((p) => p.name);
    const states = Object.fromEntries(Object.entries(STATES)
      .filter(([state, [prop]]) => c.visual !== false && names.includes(prop) &&
        !(UNREACHED[state] ?? []).includes(c.tag))
      .map(([state, [, props]]) => [state, props]));
    return {...c, ...kindOfTag.has(c.tag) ? {kind: kindOfTag.get(c.tag)} : {},
      ...Object.keys(states).length > 0 ? {states} : {}};
  });
  const propsOf = (tag) => components.find((c) => c.tag === tag).props
    .filter((p) => p.structural !== true).map((p) => p.name);
  const {version} = JSON.parse(readFileSync(new URL('../package.json', import.meta.url), 'utf8'));
  return {version, rules: {count: U2_RULES.length, source: 'src/core/rules.ts', docs: ['README.md', 'AGENTS.md']},
    ...manifest, components, kinds: kindCatalog(U2_KINDS, tagOf, propsOf, kindNotes(U2_RULES)), kindsVersion: 1};
}

const isMain = process.argv[1] !== undefined &&
  realpathSync(process.argv[1]) === realpathSync(fileURLToPath(import.meta.url));

/** The checked-in copy is CRLF in a Windows checkout, the generated string always LF. */
function compare(path, text, what) {
  const lf = (s) => s.replace(/\r\n/g, '\n');
  let actual = '';
  try {
    actual = lf(readFileSync(path, 'utf8'));
  }
  catch {
    // a missing file is as stale as an old one
  }
  if (actual === lf(text))
    return true;
  const a = actual.split('\n');
  const b = lf(text).split('\n');
  let i = 0;
  while (i < a.length && i < b.length && a[i] === b[i]) i++;
  console.error(`${what} is stale — run \`npm run manifest\` and commit the result.`);
  console.error(`  first difference at line ${i + 1}:`);
  console.error(`    checked in: ${i < a.length ? a[i] : '<end of file>'}`);
  console.error(`    generated:  ${i < b.length ? b[i] : '<end of file>'}`);
  return false;
}

if (isMain) {
  const from = process.argv.indexOf('--from');
  const base = from < 0 ? '../dist/src/' : `../${process.argv[from + 1]}/src/`;
  const manifest = await buildManifest(base);
  const {U2_RULES: rules} = await import(new URL(base + 'core/rules.js', import.meta.url).href);
  const file = (name) => fileURLToPath(new URL(`../${name}`, import.meta.url));
  const outputs = [
    ['manifest.json', JSON.stringify(manifest, null, 2) + '\n'],
    ['README.md', withRuleBlock(readFileSync(file('README.md'), 'utf8'), readmeRules(rules))],
    ['AGENTS.md', withRuleBlock(readFileSync(file('AGENTS.md'), 'utf8'), agentsRules(rules))],
  ];
  const summary = `${manifest.components.length} components, ${manifest.kinds.length} kinds, ${rules.length} rules`;

  if (!process.argv.includes('--check')) {
    for (const [name, text] of outputs)
      writeFileSync(file(name), text);
    console.log(`${summary} → ${outputs.map(([name]) => name).join(', ')}`);
  } else {
    if (!outputs.every(([name, text]) => compare(file(name), text, name)))
      process.exit(1);
    console.log(`manifest.json, README.md and AGENTS.md are up to date (${summary})`);
  }
}
