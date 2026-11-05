import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {createRequire} from 'node:module';

const ts = createRequire(import.meta.url)('typescript');
const root = fileURLToPath(new URL('..', import.meta.url)).replace(/\\/g, '/');
const BUILT = `${root}.test-build/src/`;

function entries() {
  const {exports} = JSON.parse(readFileSync(join(root, 'package.json'), 'utf8'));
  return Object.entries(exports).filter(([key, v]) => !key.includes('*') && typeof v === 'object' && v.types)
    .map(([key, v]) => ({key, file: v.types.replace(/^\.\/dist\/src\//, BUILT)}));
}

/** Every value a page script can name without importing it: the browser's window members and the language's globals. */
function globals() {
  const name = `${root}globals.ts`;
  const options = {lib: ['lib.es2023.d.ts', 'lib.dom.d.ts', 'lib.dom.iterable.d.ts'], types: [], noEmit: true};
  const host = ts.createCompilerHost(options);
  const original = host.getSourceFile.bind(host);
  host.getSourceFile = (file, language) => file === name ? ts.createSourceFile(file, '', language) : original(file, language);
  const prog = ts.createProgram([name], options, host);
  const values = prog.getTypeChecker().getSymbolsInScope(prog.getSourceFile(name), ts.SymbolFlags.Value);
  return new Set(values.map((s) => s.name));
}

const GLOBALS = globals();
const OPTIONS = {noLib: true, types: [], noEmit: true, skipLibCheck: true,
  module: ts.ModuleKind.ESNext, moduleResolution: ts.ModuleResolutionKind.Bundler};

function shadowing(overrides = new Map()) {
  const host = ts.createCompilerHost(OPTIONS);
  const original = host.getSourceFile.bind(host);
  host.getSourceFile = (file, language) =>
    overrides.has(file) ? ts.createSourceFile(file, overrides.get(file), language, true) : original(file, language);
  const all = entries().filter((e, i, list) => list.findIndex((o) => o.file === e.file) === i);
  const prog = ts.createProgram(all.map((e) => e.file), OPTIONS, host);
  const checker = prog.getTypeChecker();
  const typeOnly = (s) => (s.declarations ?? []).some((d) => ts.isExportSpecifier(d) && (d.isTypeOnly || d.parent.parent.isTypeOnly));
  const out = [];
  for (const {key, file} of all) {
    const sf = prog.getSourceFile(file);
    assert.ok(sf, `${file} is missing: run \`npm test\`, which compiles .test-build first`);
    for (const s of checker.getExportsOfModule(checker.getSymbolAtLocation(sf))) {
      const target = s.flags & ts.SymbolFlags.Alias ? checker.getAliasedSymbol(s) : s;
      if (target.flags & ts.SymbolFlags.Value && !typeOnly(s) && GLOBALS.has(s.name))
        out.push(`${key}: ${s.name} shadows the global \`${s.name}\` in every module that imports it — rename it`);
    }
  }
  return out;
}

test('the globals list holds the window members a careless name would hide', () => {
  for (const name of ['confirm', 'alert', 'prompt', 'close', 'open', 'print', 'name', 'status', 'event', 'focus',
    'Text', 'Option', 'Image', 'File', 'Range', 'Selection', 'Notification', 'Map', 'Set'])
    assert.ok(GLOBALS.has(name), name);
});

test('no u2 export shadows a DOM or language global', () => {
  assert.deepEqual(shadowing(), []);
});

test('mutation: exporting `confirm` fails (the CR4 shadowing case)', () => {
  const confirm = `${BUILT}components/containers/confirm.d.ts`;
  const index = `${BUILT}index.d.ts`;
  const found = shadowing(new Map([
    [confirm, `${readFileSync(confirm, 'utf8')}\nexport declare function confirm(message: string): Promise<boolean>;\n`],
    [index, `${readFileSync(index, 'utf8')}\nexport { confirm } from './components/containers/confirm.js';\n`],
  ]));
  assert.deepEqual(found, ['.: confirm shadows the global `confirm` in every module that imports it — rename it']);
});

test('every input\'s options type is named from the root: a bag a consumer cannot name gets retyped by hand', () => {
  const host = ts.createCompilerHost(OPTIONS);
  const index = `${BUILT}index.d.ts`;
  const prog = ts.createProgram([index], OPTIONS, host);
  const checker = prog.getTypeChecker();
  const exported = new Set(checker.getExportsOfModule(checker.getSymbolAtLocation(prog.getSourceFile(index)))
    .map((s) => s.name));
  const dir = `${BUILT}components/inputs/`;
  const missing = [];
  for (const file of readdirSync(dir).filter((f) => f.endsWith('.d.ts'))) {
    for (const [, name] of readFileSync(join(dir, file), 'utf8').matchAll(/^export interface (\w+Options)\b/gm)) {
      if (!exported.has(name))
        missing.push(`${file}: ${name}`);
    }
  }
  assert.deepEqual(missing, []);
});

test('mutation: a type that shares a DOM name is not a shadowed value', () => {
  const elements = `${BUILT}core/elements.d.ts`;
  const text = readFileSync(elements, 'utf8');
  assert.deepEqual(shadowing(new Map([[elements, `${text}\nexport type Selection = string[];\n` +
    'export interface Range { from: number; to: number; }\nexport type { Control as Option } from \'./component.js\';\n']])), []);
  assert.deepEqual(shadowing(new Map([[elements, `${text}\nexport { Control as Option } from './component.js';\n`]])),
    ['.: Option shadows the global `Option` in every module that imports it — rename it']);
});
