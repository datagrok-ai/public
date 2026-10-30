import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {createRequire} from 'node:module';

const ts = createRequire(import.meta.url)('typescript');
const root = fileURLToPath(new URL('..', import.meta.url)).replace(/\\/g, '/');
const BUILT = `${root}.test-build/src/`;
const ALLOW_FILE = 'tests/api-conventions.allow.json';
const ALLOW = JSON.parse(readFileSync(join(root, ALLOW_FILE), 'utf8'));

/** The `.d.ts` of every entry point package.json names, as `npm test` emits it into .test-build. */
function entries() {
  const {exports} = JSON.parse(readFileSync(join(root, 'package.json'), 'utf8'));
  return Object.entries(exports).filter(([key, v]) => !key.includes('*') && typeof v === 'object' && v.types)
    .map(([, v]) => v.types.replace(/^\.\/dist\/src\//, BUILT));
}

const OPTIONS = {noLib: true, types: [], noEmit: true, skipLibCheck: true,
  module: ts.ModuleKind.ESNext, moduleResolution: ts.ModuleResolutionKind.Bundler};
const parsed = new Map();

function program(overrides = new Map()) {
  const host = ts.createCompilerHost(OPTIONS);
  const original = host.getSourceFile.bind(host);
  host.getSourceFile = (name, language) => {
    if (overrides.has(name))
      return ts.createSourceFile(name, overrides.get(name), language, true);
    if (!parsed.has(name))
      parsed.set(name, original(name, language));
    return parsed.get(name);
  };
  return ts.createProgram(entries(), OPTIONS, host);
}

const FIX = {
  'owner-positional': 'take the owner as the receiver (`owner.subscribe(…)`) or as `{owner}`, never positionally',
  'no-enabled': 'say why it is off: `disabled: Reason` (the text the user sees, null while on, or a signal of it), ' +
    'never `enabled`',
  'disabled-reason': 'a `disabled` carries the reason the user sees (`Reason`), never a boolean',
  'caption-role': 'captions by role: `label` for a field, `text` for an action, button or menu item, `title` for a ' +
    'container; `name` is the automation id',
  'value-signal': 'a member named `value` is the plain value; the signal is `.signal`',
  'sort-order': 'a table\'s order is the SortOrder string (\'name\', \'!name\', \'name,!date\'; ' +
    'src/components/collections/sort-order.ts), held as a `string` or a signal of one: one contract for every table',
  'deprecated': 'u2 keeps no aliases: move every caller to what its `@deprecated` names and delete it; until then it ' +
    'is allowlisted with what removes it',
};

const ROLES = [
  {caption: 'text', type: /(?:Action|ActionGroup|Button|ButtonGroupItem|ButtonOptions|MenuItem|MenuItemOptions|MenuBarItem|Link)$/},
  {caption: 'label', type: /(?:Input|Combobox|TypeAhead|Slider)\w*Options$/},
  {caption: 'title', type: /^(?:Dialog|ConfirmDialog|Section|Card|Accordion|AccordionPane|Wizard|WizardStep|Tab|TabStrip|Panel|Page|AppHeader)(?:Options)?$/},
];
const CAPTIONS = ['label', 'text', 'title', 'caption'];
const SIGNALS = new Set(['Signal', 'ReadonlySignal', 'WritableSignal', 'Computed']);
const OWNER_TYPES = new Set(['Owner', 'Scope']);

function violations(prog) {
  const checker = prog.getTypeChecker();
  const out = new Map();
  const add = (rule, where) => out.set(`${rule}|${where}`, {rule, where});
  const inU2 = (node) => node.getSourceFile().fileName.startsWith(BUILT);
  const module = (node) => node.getSourceFile().fileName.slice(BUILT.length).replace(/\.d\.ts$/, '');
  const hidden = (node) => ts.getCombinedModifierFlags(node) & (ts.ModifierFlags.Private | ts.ModifierFlags.Protected);
  const typeText = (node) => checker.typeToString(checker.getTypeAtLocation(node.name), undefined,
    ts.TypeFormatFlags.NoTruncation);
  const ownerTyped = (type) => type !== undefined && (ts.isParenthesizedTypeNode(type) ? ownerTyped(type.type) :
    ts.isUnionTypeNode(type) ? type.types.some(ownerTyped) :
      ts.isTypeReferenceNode(type) && OWNER_TYPES.has(type.typeName.getText()));
  const deprecated = (node) => ts.getJSDocTags(node).some((t) => t.tagName.text === 'deprecated');
  const isSignal = (t) => SIGNALS.has((t.aliasSymbol ?? t.getSymbol())?.name);
  const signalTyped = (node) => {
    const t = checker.getTypeAtLocation(node.name);
    const parts = t.isUnion() ? t.types.filter((p) => !(p.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null))) : [t];
    return parts.length > 0 && parts.every(isSignal);
  };

  const params = (where, fn) => {
    for (const p of fn.parameters) {
      const name = p.name.getText();
      if (name === 'owner' || name === 'scope' || ownerTyped(p.type))
        add('owner-positional', `${where}(${name})`);
      if (name === 'enabled' && !out.has(`no-enabled|${where}`))
        add('no-enabled', `${where}(${name})`);
      if (name === 'disabled' && /\b(boolean|true|false)\b/.test(typeText(p)))
        add('disabled-reason', `${where}(${name})`);
    }
  };

  const member = (md, holder, names) => {
    const name = md.name?.getText();
    if (name === undefined || hidden(md))
      return;
    const where = `${module(md)}: ${holder}.${name}`;
    if (deprecated(md))
      add('deprecated', where);
    if (/^(?:enabled|\w+Enabled)$/.test(name))
      add('no-enabled', where);
    if (ts.isMethodDeclaration(md) || ts.isMethodSignature(md))
      params(where, md);
    const valued = ts.isPropertySignature(md) || ts.isPropertyDeclaration(md) || ts.isGetAccessor(md) || ts.isSetAccessor(md);
    if (!valued)
      return;
    if (name === 'disabled' && /\b(boolean|true|false)\b/.test(typeText(md)))
      add('disabled-reason', where);
    if (name === 'value' && signalTyped(md))
      add('value-signal', where);
    if (name === 'sort' && !/^(?:(?:Readonly)?Signal<string>|string)(?: \| undefined)?$/.test(typeText(md)))
      add('sort-order', where);
    const role = ROLES.find((r) => r.type.test(holder));
    if (role && CAPTIONS.includes(name) && name !== role.caption)
      add('caption-role', where);
    if (role?.caption === 'text' && name === 'name' && !CAPTIONS.some((c) => names.has(c)))
      add('caption-role', where);
  };

  const holderName = (md, fallback) => {
    const p = md.parent;
    if ((ts.isClassDeclaration(p) || ts.isInterfaceDeclaration(p)) && p.name)
      return p.name.text;
    return ts.isTypeLiteralNode(p) && ts.isTypeAliasDeclaration(p.parent) ? p.parent.name.text : fallback;
  };

  const seen = new Set();
  for (const file of entries()) {
    const sf = prog.getSourceFile(file);
    assert.ok(sf, `${file} is missing: run \`npm test\`, which compiles .test-build first`);
    for (let s of checker.getExportsOfModule(checker.getSymbolAtLocation(sf))) {
      if (s.flags & ts.SymbolFlags.Alias)
        s = checker.getAliasedSymbol(s);
      for (const d of s.declarations ?? []) {
        if (!inU2(d) || seen.has(d))
          continue;
        seen.add(d);
        if (deprecated(d))
          add('deprecated', `${module(d)}: ${s.name}`);
        if (ts.isFunctionDeclaration(d))
          params(`${module(d)}: ${s.name}`, d);
        if (!(ts.isClassDeclaration(d) || ts.isInterfaceDeclaration(d) || ts.isTypeAliasDeclaration(d)))
          continue;
        const props = [...checker.getPropertiesOfType(checker.getDeclaredTypeOfSymbol(s))];
        if (ts.isClassDeclaration(d)) {
          props.push(...checker.getPropertiesOfType(checker.getTypeOfSymbolAtLocation(s, d))
            .filter((p) => p.name !== 'prototype'));
          for (const c of d.members.filter(ts.isConstructorDeclaration))
            params(`${module(d)}: new ${s.name}`, c);
        }
        const names = new Set(props.map((p) => p.name));
        for (const p of props) {
          for (const md of p.declarations ?? []) {
            if (inU2(md) && !seen.has(md)) {
              seen.add(md);
              member(md, holderName(md, s.name), names);
            }
          }
        }
      }
    }
  }
  return [...out.values()];
}

function verdict(found, allow = ALLOW) {
  const problems = [];
  const keys = new Set(found.map((v) => `${v.rule}|${v.where}`));
  for (const v of found)
    if (allow[v.rule]?.[v.where] === undefined)
      problems.push(`${v.rule}: ${v.where} — ${FIX[v.rule]}`);
  for (const [rule, sites] of Object.entries(allow))
    for (const where of Object.keys(sites))
      if (!keys.has(`${rule}|${where}`))
        problems.push(`${rule}: ${where} is allowlisted but gone — delete it from ${ALLOW_FILE} (the allowlist only shrinks)`);
  return problems;
}

const real = program();
const realFound = violations(real);

test('the emitted .d.ts keeps the API conventions, but for the allowlisted sites', () => {
  assert.deepEqual(verdict(realFound), []);
});

test('every allowlisted site names a known rule and what removes it', () => {
  for (const [rule, sites] of Object.entries(ALLOW)) {
    assert.ok(FIX[rule], `unknown rule ${rule}`);
    for (const [where, by] of Object.entries(sites))
      assert.match(by, /\S/, `${rule}: ${where}`);
  }
});

const ELEMENTS = `${BUILT}core/elements.d.ts`;
const mutated = (file, change) => violations(program(new Map([[file, change(readFileSync(file, 'utf8'))]])));
const appended = (text) => mutated(ELEMENTS, (t) => `${t}\n${text}\n`);

const MUTANTS = [
  ['no-enabled', 'export interface ProbeOptions { enabled?: boolean; }', 'core/elements: ProbeOptions.enabled'],
  ['owner-positional', 'export declare function watch(source: unknown, handler: () => void, owner?: Owner): () => void;',
    'core/elements: watch(owner)'],
  ['disabled-reason', 'export interface ProbeItem { disabled?: boolean; }', 'core/elements: ProbeItem.disabled'],
  ['caption-role', 'export interface ProbeAction { label: string; run: () => void; }', 'core/elements: ProbeAction.label'],
  ['caption-role', 'export interface ProbeMenuItem { name: string; run: () => void; }', 'core/elements: ProbeMenuItem.name'],
  ['value-signal', 'export declare class ProbeInput { readonly value: ReadonlySignal<string>; }',
    'core/elements: ProbeInput.value'],
  ['sort-order', 'export interface ProbeTableOptions { sort?: {column: string, descending: boolean}; }',
    'core/elements: ProbeTableOptions.sort'],
  ['deprecated', '/** @deprecated Use probeNew. */\nexport declare function probeOld(): void;', 'core/elements: probeOld'],
];

for (const [rule, text, where] of MUTANTS) {
  test(`mutation: ${where} fails "${rule}"`, () => {
    assert.deepEqual(verdict(appended(text)), [`${rule}: ${where} — ${FIX[rule]}`]);
  });
}

test('mutation: what is not a positional owner or a caption passes', () => {
  assert.deepEqual(verdict(appended([
    'export declare function probe(build: (scope: Owner) => void, options?: {owner?: Owner}): void;',
    'export interface ProbeStatOptions { value?: string | ReadonlySignal<string>; }',
    'export interface ProbeMenuItem { text: string; name: string; disabled?: string | null; }',
    'export declare class ProbeHidden { private enabled; protected setEnabled(on: boolean): void; }',
  ].join('\n'))), []);
});

test('mutation: a violation fixed must leave the allowlist', () => {
  const [rule, where] = ['no-enabled', 'core/elements: ProbeOptions.enabled'];
  const allow = {...ALLOW, [rule]: {...ALLOW[rule], [where]: 'the next major version drops it'}};
  assert.deepEqual(verdict(appended('export interface ProbeOptions { enabled?: boolean; }'), allow), []);
  assert.deepEqual(verdict(realFound, allow), [`${rule}: ${where} is allowlisted but gone — ` +
    `delete it from ${ALLOW_FILE} (the allowlist only shrinks)`]);
});
