// Migrates code to the u2 v1 API shape, type-checked against the u2 it compiles with today.
// node scripts/codemod/u2-v1.mjs [--project <tsconfig.json>] [--u2-src <u2 dir>] [--write] [--report <file>]
//   [--rules V1,V2,…] [--loose] [--json] [--diff] [--root <dir>] [files or dirs…]
// A consumer runs it with its own tsconfig, before it moves to the new u2: the old declarations say what to rewrite.
import {createRequire} from 'node:module';
import {existsSync, readdirSync, readFileSync, statSync, writeFileSync} from 'node:fs';
import {dirname, extname, join, relative, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const ts = createRequire(import.meta.url)('typescript');

export const RULES = {
  V1: '`x.value.value` → `x.value`',
  V2: '`x.value.peek()` → `x.peek()`',
  V3: 'bare `x.value` used as a signal → `x.signal`',
  D1: '`enabled` → `disabled: Reason`, with `U2_REASON_NEEDED` where a reason must be written; boolean `disabled` → a reason',
  O1: '`x.runInScope(fn)` → `x.build(fn)`',
  O2: '`Scope.runWith(s, fn)` → `s.build(fn)`; `Scope.runWith(undefined, …)` reported',
  O3: 'a positional owner → the receiver (`o.subscribe(src, h)`) or `{owner}`',
  A1: 'captions by role: an action\'s, a button group item\'s, a menu bar item\'s `name`/`label` → `text`; a tab\'s `label` → `title`',
  W1: 'the view wrappers → `ViewHost`: `tableViewChrome(v, o)` → `ViewHost.of(v).ribbon(o)`, `showInView(d, v)` → ' +
    '`ViewHost.of(v).show(d)`, `viewSignals(v)` → `ViewHost.of(v).signals`; `ownView`, `appView` reported',
  I1: 'a path inside u2 (`@datagrok-libraries/u2/src/…`) → the entry point that exports every name the import takes',
  R1: 'an action\'s or a ribbon\'s `run`, an input\'s `onEnter`/`onEscape` handing back a value → `() => void x()`, ' +
    'a promise of one → `async () => { await x(); }`',
  S1: 'spec documents (JSON files, ```json blocks): a u2 node\'s `enabled` prop → `disabled`, with "U2_REASON_NEEDED" ' +
    'where a reason must be written; a bound one renamed and reported',
};

/** The identifier the codemod writes where only a person can say why something is off: undeclared, so the build
 * fails until it is replaced by the reason the user reads. */
export const REASON_NEEDED = 'U2_REASON_NEEDED';

const VALUE_CLASSES = new Set(['Input', 'AsyncValue', 'ProgressBar', 'Combobox', 'StateSource']);
const U2_FILE = /[\\/](?:libraries[\\/]u2|@datagrok-libraries[\\/]u2)[\\/](?:src|dist|\.test-build)[\\/]/;
const CODE = new Set(['.ts', '.tsx', '.js', '.mjs', '.cjs']);
const SKIP_DIRS = new Set(['node_modules', 'dist', '.test-build', '.git']);

const DEFAULT_OPTIONS = {
  target: ts.ScriptTarget.ES2022, module: ts.ModuleKind.ESNext, moduleResolution: ts.ModuleResolutionKind.Bundler,
  allowJs: true, checkJs: false, noEmit: true, skipLibCheck: true, esModuleInterop: true, experimentalDecorators: true,
  useDefineForClassFields: false, lib: ['lib.esnext.d.ts', 'lib.dom.d.ts', 'lib.dom.iterable.d.ts'], types: [],
};

function listFiles(path, out = []) {
  const st = statSync(path);
  if (st.isFile()) {
    out.push(resolve(path));
    return out;
  }
  for (const e of readdirSync(path, {withFileTypes: true})) {
    if (e.isDirectory()) {
      if (!SKIP_DIRS.has(e.name))
        listFiles(join(path, e.name), out);
    }
    else if ((CODE.has(extname(e.name)) || e.name.endsWith('.md') || e.name.endsWith('.json')) &&
      !/\.d\.[mc]?ts$/.test(e.name))
      out.push(resolve(path, e.name));
  }
  return out;
}

/** A JSON file worth reading for specs: it names a u2 tag. */
const holdsSpec = (text) => text.includes('"tag"') && text.includes('"u2-');

/** The ```ts (or ```json) blocks of a Markdown file as virtual sources next to it, each with its offset in the file. */
function markdownBlocks(file, text, langs = 'ts|typescript', ext = 'ts') {
  const blocks = [];
  const fence = new RegExp(`^([ \\t]*)\`\`\`(${langs})\\b[^\\n]*\\r?\\n([\\s\\S]*?)^\\1\`\`\``, 'gm');
  for (const m of text.matchAll(fence)) {
    const start = m.index + m[0].indexOf('\n') + 1;
    blocks.push({name: `${file}.${blocks.length}.${ext}`, start, text: m[3]});
  }
  return blocks;
}

/** u2's entry points, in the order an import that fits several takes them. */
const ENTRIES = ['@datagrok-libraries/u2', '@datagrok-libraries/u2/dg', '@datagrok-libraries/u2/filter',
  '@datagrok-libraries/u2/dg/filter', '@datagrok-libraries/u2/testing'];
const U2_SPECIFIER = /^@datagrok-libraries\/u2(\/.*)?$/;
const ENTRY_PATHS = new Set(['', '/dg', '/filter', '/dg/filter', '/testing', '/package.json']);

/** The subpaths of u2's explicit `exports` keys outside `./src/` (its test helpers, its stylesheets): what an import
 * may name beside the entry points. `file` is a file of the u2 package. */
function exportedPaths(file) {
  for (let dir = dirname(file); dir !== dirname(dir); dir = dirname(dir)) {
    const manifest = join(dir, 'package.json');
    const pkg = existsSync(manifest) ? JSON.parse(readFileSync(manifest, 'utf8')) : null;
    if (pkg?.name === '@datagrok-libraries/u2') {
      return Object.keys(pkg.exports ?? {}).filter((k) => !k.includes('*') && !/^\.\/src(\/|$)/.test(k))
        .map((k) => k.slice(1));
    }
  }
  return [];
}

/** The entry a path inside u2 belongs to, as the `u2/import-entry` lint rule names it: the first one tried. */
function entryFor(subpath) {
  if (/^\/src\/(core\/)?filter(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/filter';
  if (/^\/src\/dg\/filter(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/dg/filter';
  if (/^\/src\/dg(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/dg';
  if (/^\/(src\/)?test(s|ing)(\/|$)/.test(subpath))
    return '@datagrok-libraries/u2/testing';
  return '@datagrok-libraries/u2';
}

function u2SrcPaths(u2) {
  return {baseUrl: u2, paths: {
    '@datagrok-libraries/u2': ['src/index.ts'],
    '@datagrok-libraries/u2/dg': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/filter': ['src/core/filter/index.ts'],
    '@datagrok-libraries/u2/dg/filter': ['src/dg/filter/index.ts'],
    '@datagrok-libraries/u2/testing': ['src/testing/index.ts'],
    '@datagrok-libraries/u2/src/dg/index.js': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/src/dg': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/*': ['*'],
    'datagrok-api/dg': ['node_modules/datagrok-api/dist/dg.d.ts'],
    'datagrok-api/grok': ['node_modules/datagrok-api/dist/grok.d.ts'],
    'datagrok-api/ui': ['node_modules/datagrok-api/dist/ui.d.ts'],
    'datagrok-api/u2core': ['node_modules/datagrok-api/dist/src/u2core/index.d.ts'],
  }};
}

function createProgram({project, u2Src, files}) {
  let options = {...DEFAULT_OPTIONS};
  let roots = [];
  if (project) {
    const parsed = ts.getParsedCommandLineOfConfigFile(resolve(project), {}, {...ts.sys, onUnRecoverableConfigFileDiagnostic: (d) => {
      throw new Error(ts.flattenDiagnosticMessageText(d.messageText, '\n'));
    }});
    options = {...parsed.options, noEmit: true, allowJs: true, checkJs: false, skipLibCheck: true};
    roots = parsed.fileNames;
  }
  if (u2Src)
    options = {...options, ...u2SrcPaths(resolve(u2Src))};
  const virtual = new Map();
  const blocksOf = new Map();
  const code = [];
  for (const f of files) {
    if (f.endsWith('.md')) {
      const blocks = markdownBlocks(f, readFileSync(f, 'utf8'));
      blocksOf.set(f, blocks);
      for (const b of blocks)
        virtual.set(b.name, b.text);
    }
    else if (!f.endsWith('.json'))
      code.push(f);
  }
  // the entries resolved as the code being migrated resolves them, so I1 compares the very symbols it imports
  const from = join(project ? dirname(resolve(project)) : dirname(code[0] ?? files[0] ?? resolve('x')), '__u2__.ts');
  const entries = new Map();
  for (const spec of ENTRIES) {
    const resolved = ts.resolveModuleName(spec, from, options, ts.sys).resolvedModule;
    if (resolved !== undefined)
      entries.set(spec, resolved.resolvedFileName);
  }
  const host = ts.createCompilerHost(options, true);
  const getSourceFile = host.getSourceFile.bind(host);
  host.getSourceFile = (name, version, onError, create) => virtual.has(name) ?
    ts.createSourceFile(name, virtual.get(name), version, true, ts.ScriptKind.TS) :
    getSourceFile(name, version, onError, create);
  const fileExists = host.fileExists.bind(host);
  host.fileExists = (name) => virtual.has(name) || fileExists(name);
  const readFile = host.readFile.bind(host);
  host.readFile = (name) => virtual.get(name) ?? readFile(name);
  const program = ts.createProgram([...new Set([...roots, ...code, ...virtual.keys(), ...entries.values()])], options,
    host);
  return {program, blocksOf, entries};
}

/** What a rule needs about the program it runs in. */
class Context {
  constructor(program, options, entries = new Map()) {
    this.program = program;
    this.checker = program.getTypeChecker();
    this.loose = options.loose;
    this.rules = options.rules;
    this.entries = entries;
    const main = entries.get('@datagrok-libraries/u2');
    this.entryPaths = new Set([...ENTRY_PATHS, ...main === undefined ? [] : exportedPaths(main)]);
    this.exports = new Map();
    this.diagnostics = new Map();
  }

  /** A bare `x.value` of a u2 v0 value class: 'signal' where the code asks for a signal there (a signal-typed slot,
   * `.subscribe`), 'v1' where the v0 type check rejects the expression it sits in — code in the v1 shape run in the
   * old world — and 'unclear' where nothing says which. */
  signalUse(node) {
    if (this.rejected(node))
      return 'v1';
    if (isAccess(node.parent, 'subscribe') && node.parent.expression === node)
      return 'signal';
    const wanted = this.checker.getContextualType(node);
    const parts = wanted === undefined ? [] : wanted.isUnion() ? wanted.types : [wanted];
    return parts.some((t) => SIGNAL_TYPES.has((t.aliasSymbol ?? t.getSymbol())?.name) ||
      t.getProperty('subscribe') !== undefined) ? 'signal' : 'unclear';
  }

  /** The v0 type check reports an error in the expression `node` is part of, up to an argument, a statement or a
   * function. */
  rejected(node) {
    // a Markdown snippet rarely compiles whole: its errors say nothing about its shape
    if (/\.md\.\d+\.ts$/.test(node.getSourceFile().fileName))
      return false;
    let unit = node;
    while (ts.isExpression(unit.parent) && !ts.isFunctionLike(unit.parent) &&
      !((ts.isCallExpression(unit.parent) || ts.isNewExpression(unit.parent)) && unit.parent.arguments?.includes(unit)))
      unit = unit.parent;
    const sf = node.getSourceFile();
    if (!this.diagnostics.has(sf))
      this.diagnostics.set(sf, this.program.getSemanticDiagnostics(sf));
    const start = unit.getStart(sf);
    return this.diagnostics.get(sf).some((d) => d.start !== undefined && d.start < unit.end &&
      d.start + (d.length ?? 0) > start);
  }

  /** What an entry point exports — name → the symbol it stands for — or null where it did not resolve. */
  entryExports(spec) {
    if (!this.exports.has(spec)) {
      const sf = this.entries.has(spec) ? this.program.getSourceFile(this.entries.get(spec)) : undefined;
      const module = sf && this.checker.getSymbolAtLocation(sf);
      this.exports.set(spec, module ?
        new Map(this.checker.getExportsOfModule(module).map((s) => [s.name, this.target(s)])) : null);
    }
    return this.exports.get(spec);
  }

  target(symbol) {
    return symbol.flags & ts.SymbolFlags.Alias ? this.checker.getAliasedSymbol(symbol) : symbol;
  }

  holder(decl) {
    const p = decl.parent;
    return (ts.isClassDeclaration(p) || ts.isInterfaceDeclaration(p)) && p.name ? p.name.text : undefined;
  }

  inU2(decl) {
    return U2_FILE.test(decl.getSourceFile().fileName.replace(/\\/g, '/'));
  }

  /** A member u2 declares on a class, an interface or a named type — not a property of an object literal, nor a
   * platform shape u2 describes. */
  u2Member(d) {
    return this.inU2(d) && !PLATFORM_SHAPES.has(this.holder(d)) && (ts.isPropertySignature(d) ||
      ts.isPropertyDeclaration(d) || ts.isGetAccessor(d) || ts.isSetAccessor(d) || ts.isMethodDeclaration(d) ||
      ts.isMethodSignature(d));
  }

  declarations(name) {
    return this.checker.getSymbolAtLocation(name)?.declarations ?? [];
  }

  isSignalField(d) {
    return (ts.isPropertyDeclaration(d) || ts.isPropertySignature(d)) && VALUE_CLASSES.has(this.holder(d)) && this.inU2(d);
  }

  /** `x.value` is the signal member u2 v0 declares on its value classes (a field; v1 has an accessor there):
   * 'typed', or 'helper' where the untyped `x` was built as `const x = mount(new TextInput(…))` through a helper
   * of the same file that answers its argument. */
  signalValue(access) {
    if (access.name.text !== 'value')
      return null;
    if (this.declarations(access.name).some((d) => this.isSignalField(d)))
      return 'typed';
    if (!this.isAny(access.expression))
      return null;
    const built = this.builtThroughHelper(access.expression);
    const prop = built && this.checker.getTypeAtLocation(built).getProperty('value');
    return prop?.declarations?.some((d) => this.isSignalField(d)) ? 'helper' : null;
  }

  builtThroughHelper(expr) {
    if (!ts.isIdentifier(expr))
      return null;
    const decl = this.checker.getSymbolAtLocation(expr)?.valueDeclaration;
    let init = decl && ts.isVariableDeclaration(decl) && ts.isIdentifier(decl.name) ? decl.initializer : undefined;
    for (let depth = 0; init && depth < 4; depth++) {
      while (ts.isAwaitExpression(init) || ts.isParenthesizedExpression(init))
        init = init.expression;
      if (ts.isNewExpression(init))
        return init;
      if (!ts.isCallExpression(init) || init.arguments.length === 0 || !this.answersArgument(init.expression))
        return null;
      init = init.arguments[0];
    }
    return null;
  }

  answersArgument(callee) {
    const fn = ts.isIdentifier(callee) ? this.checker.getSymbolAtLocation(callee)?.valueDeclaration : undefined;
    const body = fn && (ts.isFunctionDeclaration(fn) ? fn : ts.isVariableDeclaration(fn) &&
      fn.initializer && (ts.isArrowFunction(fn.initializer) || ts.isFunctionExpression(fn.initializer)) ? fn.initializer : null);
    const param = body?.parameters[0]?.name;
    if (!param || !ts.isIdentifier(param) || body.getSourceFile() !== callee.getSourceFile())
      return false;
    if (!ts.isBlock(body.body))
      return ts.isIdentifier(body.body) && body.body.text === param.text;
    const returns = [];
    const visit = (n) => {
      if (ts.isReturnStatement(n))
        returns.push(n.expression);
      else if (!ts.isFunctionLike(n))
        ts.forEachChild(n, visit);
    };
    ts.forEachChild(body.body, visit);
    return returns.length > 0 && returns.every((r) => r && ts.isIdentifier(r) && r.text === param.text);
  }

  isAny(node) {
    return (this.checker.getTypeAtLocation(node).flags & ts.TypeFlags.Any) !== 0;
  }

  /** `expr` is typed by an anonymous `{value: {value: …}}` the file wrote itself — the shape a cast over a u2 control
   * (`Control.forElement(el) as unknown as {value: …}`) gives it. */
  localValueShape(expr) {
    const type = this.checker.getNonNullableType(this.checker.getTypeAtLocation(expr));
    const decls = type.getProperty('value')?.declarations ?? [];
    return decls.length > 0 && decls.every((d) => ts.isPropertySignature(d) && ts.isTypeLiteralNode(d.parent) &&
      !this.inU2(d) && d.type !== undefined && ts.isTypeLiteralNode(d.type) &&
      d.type.members.some((m) => m.name?.getText() === 'value'));
  }

  /** `type`, or a member of a union, has a member `name` that u2 declares. */
  hasU2Member(type, ...names) {
    if (type === undefined || type.flags & ts.TypeFlags.Any)
      return false;
    const parts = type.isUnion() ? type.types : [type];
    return parts.some((t) => names.some((name) => this.checker.getApparentType(t).getProperty(name)?.declarations
      ?.some((d) => this.u2Member(d))));
  }

  isSignalType(type) {
    const parts = type.isUnion() ? type.types.filter((t) => !(t.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null))) :
      [type];
    return parts.length > 0 && parts.every((t) => SIGNAL_TYPES.has((t.aliasSymbol ?? t.getSymbol())?.name));
  }

  /** A live option (`T | ReadonlySignal<T>`): some member of `type` is a signal, some is not. */
  maybeSignal(type) {
    const parts = type.isUnion() ? type.types : [type];
    return !this.isSignalType(type) && parts.some((t) => SIGNAL_TYPES.has((t.aliasSymbol ?? t.getSymbol())?.name));
  }

  /** `node` may be undefined: its type says so or, where no strict null checks keep that in the type, it names an
   * optional property or parameter. */
  optional(node) {
    if (nullable(this, node))
      return true;
    const symbol = ts.isShorthandPropertyAssignment(node.parent) && node.parent.name === node ?
      this.checker.getShorthandAssignmentValueSymbol(node.parent) :
      this.checker.getSymbolAtLocation(ts.isPropertyAccessExpression(node) ? node.name : node);
    const decl = symbol?.valueDeclaration;
    return symbol !== undefined && ((symbol.flags & ts.SymbolFlags.Optional) !== 0 ||
      decl !== undefined && ts.isParameter(decl) && decl.questionToken !== undefined);
  }

  /** A u2 options shape whose `disabled` is (or becomes) a reason: it takes a string — or was a tree node's flag —
   * and is no state block's boolean. */
  takesReason(type) {
    if (type === undefined || type.flags & ts.TypeFlags.Any)
      return false;
    const parts = type.isUnion() ? type.types : [type];
    return parts.some((t) => {
      const prop = this.checker.getApparentType(t).getProperty('disabled');
      const decls = prop?.declarations?.filter((d) => this.u2Member(d)) ?? [];
      if (decls.length === 0 || decls.some((d) => STATE_SHAPES.has(this.holder(d))))
        return false;
      const takes = this.checker.getTypeOfSymbol(prop);
      const strings = (takes.isUnion() ? takes.types : [takes]).some((p) => (p.flags & ts.TypeFlags.StringLike) !== 0);
      return strings || decls.some((d) => FLAG_TO_REASON.has(this.holder(d)));
    });
  }

  /** As {@link hasU2Member}, declared by a u2 class (an accessor that always holds a value), not an options shape. */
  classMember(type, ...names) {
    if (type === undefined || type.flags & ts.TypeFlags.Any)
      return false;
    const parts = type.isUnion() ? type.types : [type];
    return parts.some((t) => names.some((name) => this.checker.getApparentType(t).getProperty(name)?.declarations
      ?.some((d) => this.u2Member(d) && ts.isClassDeclaration(d.parent))));
  }

  signalValueType(type) {
    const part = type.isUnion() ? type.types.find((t) => !(t.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null))) : type;
    const value = part?.getProperty('value');
    return value === undefined ? this.checker.getAnyType() : this.checker.getTypeOfSymbol(value);
  }

  /** Booleans beside strings: an off that is sometimes a reason and sometimes a flag. */
  mixesFlags(type) {
    const parts = type.isUnion() ? type.types : [type];
    return parts.some((t) => (t.flags & ts.TypeFlags.BooleanLike) !== 0) &&
      parts.some((t) => (t.flags & ts.TypeFlags.StringLike) !== 0);
  }

  isBooleanType(type) {
    const parts = type.isUnion() ? type.types.filter((t) => !(t.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null))) :
      [type];
    return parts.length > 0 && parts.every((t) => (t.flags & ts.TypeFlags.BooleanLike) !== 0);
  }
}

const SIGNAL_TYPES = new Set(['Signal', 'ReadonlySignal', 'Computed', 'WritableSignal']);
/** The state protocol's own `disabled` flags: booleans by contract. */
const STATE_SHAPES = new Set(['IItemState', 'IControlState', 'IItemsState']);
/** Platform objects u2 types structurally: their `enabled` is the platform's. */
const PLATFORM_SHAPES = new Set(['DartInputLike']);
/** u2 v0 shapes whose boolean `disabled` becomes a reason in v1. */
const FLAG_TO_REASON = new Set(['TreeNode']);

const isAccess = (n, name) => n !== undefined && ts.isPropertyAccessExpression(n) && n.name.text === name;

const importsU2 = (sf) => sf.statements.some((s) => (ts.isImportDeclaration(s) || ts.isExportDeclaration(s)) &&
  s.moduleSpecifier !== undefined && ts.isStringLiteral(s.moduleSpecifier) &&
  U2_SPECIFIER.test(s.moduleSpecifier.text));

/** Passing the signal where its text was meant reads the same: a template, `String(x)`, a `+` with a string. */
function coerced(node) {
  const p = node.parent;
  if (ts.isTemplateSpan(p) && p.expression === node)
    return true;
  if (ts.isCallExpression(p) && ts.isIdentifier(p.expression) && ['String', 'Number', 'Boolean'].includes(p.expression.text))
    return true;
  return ts.isBinaryExpression(p) && p.operatorToken.kind === ts.SyntaxKind.PlusToken &&
    [p.left, p.right].some((x) => ts.isStringLiteralLike(x) || ts.isTemplateExpression(x));
}

function valueRules(ctx, sf, out) {
  const signals = [];
  const unclear = [];
  let v1Shaped = false;
  const visit = (node) => {
    if (ts.isPropertyAccessExpression(node) && node.name.text === 'value') {
      const how = ctx.signalValue(node);
      const typed = how !== null;
      const loose = !typed && ctx.loose && /\.[mc]?js$/.test(sf.fileName) &&
        (ctx.isAny(node.expression) || ctx.declarations(node.name).length === 0);
      const parent = node.parent;
      if (typed || loose) {
        if (isAccess(parent, 'value') && parent.expression === node) {
          if (ctx.rules.has('V1'))
            out.edit('V1', node.end, parent.end, '', loose);
        }
        else if (isAccess(parent, 'peek') && parent.expression === node && ts.isCallExpression(parent.parent) &&
          parent.parent.expression === parent) {
          if (ctx.rules.has('V2'))
            out.edit('V2', node.expression.end, node.end, node.questionDotToken ? '?' : '', loose);
        }
        else if (typed && ctx.rules.has('V3') && !coerced(node)) {
          const use = ctx.signalUse(node);
          if (use === 'v1')
            v1Shaped = true;
          else
            (use === 'signal' ? signals : unclear).push(node);
        }
      }
      else if (isAccess(parent, 'value') && parent.expression === node && ctx.rules.has('V1') && importsU2(sf) &&
        ctx.localValueShape(node.expression))
        out.hand('V1', node.getStart(sf), 'a `.value.value` through a local type (`{value: {value: …}}`): if it is a u2 ' +
          'control, read `.value` — v1 holds the value itself, and the cast hides it from the type check');
      else if (ctx.isAny(node.expression) && (importsU2(sf) || /\bU2\b/.test(node.expression.getText(sf)))) {
        if (isAccess(parent, 'value') && parent.expression === node && ctx.rules.has('V1'))
          out.hand('V1', node.getStart(sf), 'a `.value.value` on an untyped (`any`) receiver: if it is a u2 control, ' +
            'read `.value` — v1 holds the value itself, and the type check cannot see it');
        else if (isAccess(parent, 'peek') && parent.expression === node && ts.isCallExpression(parent.parent) &&
          parent.parent.expression === parent && ctx.rules.has('V2'))
          out.hand('V2', node.getStart(sf), 'a `.value.peek()` on an untyped (`any`) receiver: if it is a u2 control, ' +
            'read `.peek()` — v1 holds the value itself, and the type check cannot see it');
      }
    }
    if (ts.isBindingElement(node) && ts.isObjectBindingPattern(node.parent) && ctx.rules.has('V3')) {
      const key = node.propertyName ?? node.name;
      const decl = node.parent.parent;
      if (ts.isIdentifier(key) && key.text === 'value' && (ts.isVariableDeclaration(decl) || ts.isParameter(decl))) {
        const prop = ctx.checker.getTypeAtLocation(node.parent).getProperty('value');
        if (prop?.declarations?.some((d) => VALUE_CLASSES.has(ctx.holder(d)) && ctx.inU2(d)))
          out.hand('V3', node.getStart(sf), '`{value}` destructured from a u2 value class: take `signal` (live) or read `.value`');
      }
    }
    ts.forEachChild(node, visit);
  };
  visit(sf);
  // the v0 type check rejecting some `.value` reads says the file is v1 already: its signal-shaped ones may be values
  const left = [...signals, ...unclear];
  if (v1Shaped && left.length > 0) {
    out.hand('V3', Math.min(...left.map((n) => n.getStart(sf))), `the v0 type check reads this file as v1 already: ` +
      `${left.length} bare \`.value\` left as they are — check a migrated tree against u2's v1 sources (\`--u2-src\`)`);
    return;
  }
  for (const node of signals)
    out.edit('V3', node.name.getStart(sf), node.name.end, 'signal');
  for (const node of unclear)
    out.hand('V3', node.getStart(sf), 'a bare `.value` of a u2 value class where nothing asks for a signal: ' +
      '`.signal` passes the signal on; in v1 `.value` reads the value');
}

const SIMPLE =(e) => ts.isIdentifier(e) || ts.isPropertyAccessExpression(e) || ts.isElementAccessExpression(e) ||
  ts.isCallExpression(e) || ts.isParenthesizedExpression(e) || e.kind === ts.SyntaxKind.ThisKeyword;

/** Rewrites what `on` (an `enabled` value) says into the reason `disabled` takes, as edits around `on`, so rewrites
 * inside it still apply: `true` → `null`, `false` → the placeholder, a boolean `e` → `e ? null : …` (`e === false ?
 * … : null` where `e` may be undefined, which was on), a signal `s` → `computed(() => s.value ? null : …)`; one that
 * may be either is reported. `inverse`: `on` is a boolean `disabled`, true meaning off. */
function asReason(ctx, sf, on, out, {inverse = false} = {}) {
  const off = REASON_NEEDED;
  const [yes, no] = inverse ? [off, 'null'] : ['null', off];
  if (on.kind === ts.SyntaxKind.TrueKeyword || on.kind === ts.SyntaxKind.FalseKeyword) {
    const value = on.kind === ts.SyntaxKind.TrueKeyword ? yes : no;
    out.edit('D1', on.getStart(sf), on.end, value);
    if (value === off)
      out.hand('D1', on.getStart(sf), `write the reason the user reads in place of ${off}`);
    return;
  }
  const type = ctx.checker.getTypeAtLocation(on);
  if (ctx.maybeSignal(type)) {
    out.hand('D1', on.getStart(sf), 'a value that may be a signal (a live option): a signal is always truthy, so no ' +
      'test of it says whether it is on — write the reason by hand, a `computed(() => …)` over the signal\'s value');
    return;
  }
  const signal = ctx.isSignalType(type);
  const body = ts.isCallExpression(on) && ts.isIdentifier(on.expression) && on.expression.text === 'computed' &&
    on.arguments.length === 1 && ts.isArrowFunction(on.arguments[0]) && !ts.isBlock(on.arguments[0].body) ?
    on.arguments[0].body : undefined;
  if (signal && body !== undefined) {
    asReason(ctx, sf, body, out, {inverse});
    return;
  }
  const wrap = !SIMPLE(on) && !(ts.isPrefixUnaryExpression(on) && on.operator === ts.SyntaxKind.ExclamationToken);
  if (signal) {
    out.edit('D1', on.getStart(sf), on.getStart(sf), 'computed(() => ' + (wrap ? '(' : ''));
    out.edit('D1', on.end, on.end, `${wrap ? ')' : ''}.value ? ${yes} : ${no})`);
    out.needs('computed', on.getStart(sf));
  }
  else {
    const undefinedOn = !inverse && ctx.optional(on);
    if (wrap)
      out.edit('D1', on.getStart(sf), on.getStart(sf), '(');
    out.edit('D1', on.end, on.end,
      `${wrap ? ')' : ''}${undefinedOn ? ` === false ? ${no} : ${yes}` : ` ? ${yes} : ${no}`}`);
  }
  out.hand('D1', on.getStart(sf), `write the reason the user reads in place of ${off}`);
}

/** Deletes an object-literal property with the comma that separates it. */
function dropProperty(sf, prop, out) {
  dropItem(sf, prop.parent.properties, prop, out, 'D1');
}

/** Deletes `item` from a comma-separated `list` (properties, import specifiers) with the comma that separates it. */
function dropItem(sf, list, item, out, rule) {
  const i = list.indexOf(item);
  const text = sf.text;
  const line = (pos) => sf.getLineAndCharacterOfPosition(pos).line;
  // the last item of a line after others goes with the comma before it, so the next line is not pulled up
  const endsLine = i > 0 && i < list.length - 1 && line(list[i + 1].getStart(sf)) > line(item.end) &&
    line(list[i - 1].end) === line(item.getStart(sf));
  if (i < list.length - 1 && !endsLine)
    out.edit(rule, item.getStart(sf), list[i + 1].getStart(sf), '');
  else if (i > 0)
    out.edit(rule, list[i - 1].end, item.end, '');
  else {
    let end = item.end;
    while (end < text.length && /[\s,]/.test(text[end]) && text[end] !== '\n')
      end++;
    out.edit(rule, item.getStart(sf), end, '');
  }
}

const isBoolLiteral = (e) => e.kind === ts.SyntaxKind.TrueKeyword || e.kind === ts.SyntaxKind.FalseKeyword;

/** A spec node literal: `{tag: 'u2-…', props: {…}, bind: {…}}`. */
function specNodeOf(obj) {
  const parent = obj.parent;
  if (!ts.isPropertyAssignment(parent) || !ts.isObjectLiteralExpression(parent.parent))
    return null;
  const key = parent.name.getText();
  if (key !== 'props' && key !== 'bind')
    return null;
  const tag = parent.parent.properties.find((p) => ts.isPropertyAssignment(p) && p.name.getText() === 'tag');
  return tag && ts.isStringLiteralLike(tag.initializer) && tag.initializer.text.startsWith('u2-') ? key : null;
}

function reasonRules(ctx, sf, out) {
  if (!ctx.rules.has('D1'))
    return;
  const visit = (node) => {
    if ((ts.isPropertyAssignment(node) || ts.isShorthandPropertyAssignment(node)) &&
      ts.isObjectLiteralExpression(node.parent)) {
      const key = node.name.getText(sf);
      const spec = specNodeOf(node.parent);
      const options = ctx.checker.getContextualType(node.parent);
      if (key === 'deleteEnabled' && ctx.hasU2Member(options, 'deleteEnabled', 'deletable'))
        out.edit('D1', node.name.getStart(sf), node.name.end, ts.isShorthandPropertyAssignment(node) ?
          'deletable: deleteEnabled' : 'deletable');
      else if (key === 'enabled' && spec === 'bind') {
        out.edit('D1', node.name.getStart(sf), node.name.end, 'disabled');
        out.hand('D1', node.getStart(sf), 'a bound `disabled` reads a reason (string | null), not a boolean: bind it ' +
          'to one');
      }
      else if (key === 'enabled' && (spec === 'props' || ctx.hasU2Member(options, 'enabled', 'disabled'))) {
        if (ts.isPropertyAssignment(node) && node.initializer.kind === ts.SyntaxKind.TrueKeyword)
          dropProperty(sf, node, out);
        else if (ts.isShorthandPropertyAssignment(node)) {
          out.edit('D1', node.name.getStart(sf), node.name.getStart(sf), 'disabled: ');
          asReason(ctx, sf, node.name, out);
        }
        else {
          out.edit('D1', node.name.getStart(sf), node.name.end, 'disabled');
          asReason(ctx, sf, node.initializer, out);
        }
      }
      else if (key === 'disabled' && ctx.takesReason(options)) {
        const value = ts.isPropertyAssignment(node) ? node.initializer : node.name;
        const type = ctx.checker.getTypeAtLocation(value);
        const said = ctx.isSignalType(type) ? ctx.signalValueType(type) : type;
        if (value.kind === ts.SyntaxKind.FalseKeyword)
          dropProperty(sf, node, out);
        else if (ctx.isBooleanType(said)) {
          if (ts.isShorthandPropertyAssignment(node))
            out.edit('D1', node.name.getStart(sf), node.name.getStart(sf), 'disabled: ');
          asReason(ctx, sf, value, out, {inverse: true});
        }
        else if (ctx.mixesFlags(said))
          out.hand('D1', node.getStart(sf), '`disabled` mixes booleans and reasons: answer the reason while off and ' +
            '`null` while on — never `true` or `false`');
      }
    }
    const receiver = ts.isPropertyAccessExpression(node) && node.name.text === 'enabled' &&
      !(ts.isCallExpression(node.parent) && node.parent.expression === node) ?
      ctx.checker.getTypeAtLocation(node.expression) : undefined;
    const untyped = receiver !== undefined && ctx.loose && /\.[mc]?js$/.test(sf.fileName) &&
      (receiver.flags & ts.TypeFlags.Any || ctx.declarations(node.name).length === 0);
    if (receiver !== undefined && ctx.hasU2Member(receiver, 'enabled', 'disabled') &&
      !ctx.classMember(receiver, 'enabled', 'disabled'))
      out.hand('D1', node.getStart(sf), 'an options object\'s `enabled`: read its `disabled` reason instead');
    else if (receiver !== undefined && (untyped || ctx.classMember(receiver, 'enabled', 'disabled'))) {
      const p = node.parent;
      out.edit('D1', node.name.getStart(sf), node.name.end, 'disabled');
      if (ts.isBinaryExpression(p) && p.left === node && p.operatorToken.kind === ts.SyntaxKind.EqualsToken)
        asReason(ctx, sf, p.right, out);
      else if (ts.isPrefixUnaryExpression(p) && p.operator === ts.SyntaxKind.ExclamationToken) {
        out.edit('D1', p.getStart(sf), node.getStart(sf), '');
        out.edit('D1', node.end, node.end, ' !== null');
      }
      else if (ts.isBinaryExpression(p) && [ts.SyntaxKind.EqualsEqualsEqualsToken, ts.SyntaxKind.ExclamationEqualsEqualsToken]
        .includes(p.operatorToken.kind) && (isBoolLiteral(p.left) || isBoolLiteral(p.right))) {
        const literal = p.left === node ? p.right : p.left;
        const on = (literal.kind === ts.SyntaxKind.TrueKeyword) === (p.operatorToken.kind === ts.SyntaxKind.EqualsEqualsEqualsToken);
        out.edit('D1', p.getStart(sf), node.getStart(sf), '');
        out.edit('D1', node.end, p.end, on ? ' === null' : ' !== null');
      }
      else {
        const bare = ts.isCallExpression(p) || ts.isVariableDeclaration(p) || ts.isReturnStatement(p) ||
          ts.isIfStatement(p) || ts.isPropertyAssignment(p) || ts.isArrowFunction(p) || ts.isConditionalExpression(p) &&
          p.condition === node || ts.isArrayLiteralExpression(p);
        if (!bare)
          out.edit('D1', node.getStart(sf), node.getStart(sf), '(');
        out.edit('D1', node.end, node.end, bare ? ' === null' : ' === null)');
      }
    }
    if (ts.isCallExpression(node) && ts.isPropertyAccessExpression(node.expression)) {
      const callee = node.expression;
      const receiver = ctx.checker.getTypeAtLocation(callee.expression);
      const name = callee.name.text;
      if (name === 'setEnabled' && node.arguments.length === 2 && ctx.hasU2Member(receiver, 'setEnabled', 'setDisabled')) {
        out.edit('D1', callee.name.getStart(sf), callee.name.end, 'setDisabled');
        asReason(ctx, sf, node.arguments[1], out);
      }
      else if (name === 'okEnabled' && node.arguments.length === 1 && ctx.hasU2Member(receiver, 'okEnabled', 'okDisabled')) {
        out.edit('D1', callee.name.getStart(sf), callee.name.end, 'okDisabled');
        asReason(ctx, sf, node.arguments[0], out);
      }
    }
    ts.forEachChild(node, visit);
  };
  visit(sf);
}

/** `computed` in scope for what D1 wrote: imported already, added to the file's import from the u2 entry, or a
 * hand-finish item. */
function importsFor(ctx, sf, out) {
  for (const [name, pos] of out.needed) {
    const imports = sf.statements.filter(ts.isImportDeclaration);
    const has = imports.some((d) => d.importClause?.namedBindings && ts.isNamedImports(d.importClause.namedBindings) &&
      d.importClause.namedBindings.elements.some((e) => e.name.text === name));
    if (has)
      continue;
    const u2 = imports.find((d) => ts.isStringLiteral(d.moduleSpecifier) && d.moduleSpecifier.text === '@datagrok-libraries/u2' &&
      d.importClause?.namedBindings && ts.isNamedImports(d.importClause.namedBindings) && !d.importClause.isTypeOnly);
    if (u2 !== undefined) {
      const elements = u2.importClause.namedBindings.elements;
      out.edit('D1', elements[elements.length - 1].end, elements[elements.length - 1].end, `, ${name}`);
    }
    else
      out.hand('D1', pos, `import \`${name}\` (from the u2 entry) for the reason this rewrite wrote`);
  }
  out.needed.clear();
}

const ENGINE_FILE = /[\\/](?:js-api[\\/](?:dist[\\/])?src[\\/]u2core|datagrok-api[\\/](?:dist[\\/])?src[\\/]u2core)[\\/]/;

/** The positional owners v1 takes as the receiver or as `{owner}`: [function or Holder.method, owner's index,
 * form, the fewest arguments a v0 call with the owner has]. `receiver`: `f(…, o)` → `o.f(…)`; `option`: the owner
 * becomes a last `{owner}` argument; `merge`: the owner joins the options argument after it (`also` names the plain
 * arguments that join it too). */
const POSITIONAL_OWNERS = [
  {name: 'subscribe', owner: 2, form: 'receiver', arity: 3},
  {name: 'toSignal', owner: 2, form: 'option', arity: 3},
  {name: 'bindText', owner: 0, form: 'option', arity: 3},
  {name: 'collectionSelection', owner: 0, form: 'option', arity: 3},
  {name: 'Menu.bindContext', owner: 1, form: 'option', arity: 2},
  {name: 'Tooltip.bind', owner: 2, form: 'option', arity: 3},
  {name: 'Input.bindOwned', owner: 0, form: 'option', arity: 2},
  {name: 'new ViewSync', owner: 1, form: 'option', arity: 2},
  {name: 'Overlay.show', owner: 2, form: 'merge', arity: 3},
  {name: 'new FrameRows', owner: 1, form: 'merge', arity: 2},
  {name: 'frameRows', owner: 1, form: 'merge', arity: 2},
  {name: 'bindValue', owner: 0, form: 'merge', also: ['commitOn'], arity: 3},
];

/** A call already in the v1 shape — an `{owner}` option among its arguments — or one too short to carry a v0 owner
 * (`bindText(el, source)`): what an untyped (`any`) argument would otherwise be taken for. */
const ownerless = (sf, rule, args) => args.length < rule.arity || args.some((a) => ts.isObjectLiteralExpression(a) &&
  a.properties.some((p) => p.name?.getText(sf) === 'owner'));

/** Methods v1 renames: Holder.v0 → v1. */
const METHOD_RENAMES = {'TabStrip.setTabLabel': 'setTabTitle'};

/** Captions by role (M2): per v0 shape, the key that becomes the role's caption. */
const CAPTION_RENAMES = {
  Action: {name: 'text'}, DesignerAction: {name: 'text'}, ActionGroup: {label: 'text'},
  ButtonGroupItem: {label: 'text'}, MenuBarItem: {label: 'text'}, RibbonAction: {label: 'text'},
  TabOptions: {label: 'title'},
};

function calleeName(ctx, call) {
  const target = ts.isNewExpression(call) ? call.expression : call.expression;
  const name = ts.isPropertyAccessExpression(target) ? target.name : ts.isIdentifier(target) ? target : undefined;
  if (name === undefined)
    return null;
  let symbol = ctx.checker.getSymbolAtLocation(name);
  if (symbol && symbol.flags & ts.SymbolFlags.Alias)
    symbol = ctx.checker.getAliasedSymbol(symbol);
  const own = (s) => s?.declarations?.find((d) => ctx.inU2(d) || ENGINE_FILE.test(d.getSourceFile().fileName));
  // a function or class taken off a module another way (`const {f} = await import('…')`) is known by its type
  const typed = ts.isIdentifier(target) && own(symbol) === undefined ? ctx.checker.getTypeAtLocation(name).getSymbol() :
    undefined;
  if (typed !== undefined && typed.flags & (ts.SymbolFlags.Function | ts.SymbolFlags.Class))
    symbol = typed;
  const decl = own(symbol);
  if (decl === undefined)
    return null;
  if (ts.isNewExpression(call))
    return ts.isClassDeclaration(decl) ? `new ${symbol.name}` : null;
  const holder = ctx.holder(decl);
  return holder === undefined ? symbol.name : `${holder}.${symbol.name}`;
}

/** An owner: a scope or a component (something that owns and disposes); `any` counts, a missing argument not. */
function ownerTyped(ctx, arg) {
  if (arg === undefined || ts.isObjectLiteralExpression(arg))
    return false;
  const type = ctx.checker.getTypeAtLocation(arg);
  if (type.flags & ts.TypeFlags.Any)
    return true;
  const parts = type.isUnion() ? type.types.filter((t) => !(t.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null))) : [type];
  return parts.length > 0 && parts.every((t) => t.getProperty('own') !== undefined && t.getProperty('dispose') !== undefined);
}

const nullable = (ctx, arg) => {
  const type = ctx.checker.getTypeAtLocation(arg);
  return type.isUnion() && type.types.some((t) => t.flags & (ts.TypeFlags.Undefined | ts.TypeFlags.Null));
};

function ownerRules(ctx, sf, out) {
  const moved = new Map();
  const visit = (node) => {
    if (ts.isCallExpression(node) && ts.isPropertyAccessExpression(node.expression)) {
      const name = calleeName(ctx, node);
      const args = node.arguments;
      const method = node.expression.name.text;
      const rename = Object.entries(METHOD_RENAMES).find(([from]) => from.endsWith(`.${method}`));
      if (ctx.rules.has('A1') && rename !== undefined &&
        ctx.hasU2Member(ctx.checker.getTypeAtLocation(node.expression.expression), method, rename[1]))
        out.edit('A1', node.expression.name.getStart(sf), node.expression.name.end, rename[1]);
      if (ctx.rules.has('O1') && node.expression.name.text === 'runInScope' && name !== null &&
        name.endsWith('.runInScope'))
        out.edit('O1', node.expression.name.getStart(sf), node.expression.name.end, 'build');
      if (ctx.rules.has('O2') && name === 'Scope.runWith' && args.length === 2) {
        const [owner] = args;
        if (owner.kind === ts.SyntaxKind.UndefinedKeyword || ts.isIdentifier(owner) && owner.text === 'undefined')
          out.hand('O2', node.getStart(sf), '`Scope.runWith(undefined, …)` builds for nobody: give it an owner ' +
            '(`owner.build(…)`) or start a root on purpose (`Scope.root(\'<what>\').build(…)`)');
        else if (nullable(ctx, owner))
          out.hand('O2', node.getStart(sf), '`Scope.runWith(s, …)` with an `s` that may be undefined: build under an ' +
            'owner that is there (`s.build(…)`) or say what runs unowned');
        else {
          const text = SIMPLE(owner) ? owner.getText(sf) : `(${owner.getText(sf)})`;
          out.edit('O2', node.getStart(sf), args[1].getStart(sf), `${text}.build(`);
        }
      }
    }
    if ((ts.isCallExpression(node) || ts.isNewExpression(node)) && ctx.rules.has('O3')) {
      const name = calleeName(ctx, node);
      const rule = name && POSITIONAL_OWNERS.find((r) => r.name === name);
      const args = node.arguments ?? [];
      if (rule && !ownerless(sf, rule, args) && ownerTyped(ctx, args[rule.owner]) &&
        positionalOwner(ctx, sf, node, rule, args, out) && ts.isIdentifier(node.expression)) {
        const alias = ctx.checker.getSymbolAtLocation(node.expression);
        moved.set(alias, (moved.get(alias) ?? 0) + 1);
      }
    }
    if (ctx.rules.has('A1') && (ts.isPropertyAssignment(node) || ts.isShorthandPropertyAssignment(node)) &&
      ts.isObjectLiteralExpression(node.parent))
      caption(ctx, sf, node, out);
    if (ctx.rules.has('A1') && ts.isPropertyAccessExpression(node) &&
      !(ts.isCallExpression(node.parent) && node.parent.expression === node))
      captionRead(ctx, sf, node, out);
    ts.forEachChild(node, visit);
  };
  visit(sf);
  for (const [alias, calls] of moved)
    dropMovedImport(ctx, sf, alias, calls, out);
}

/** The import of a function whose calls O3 all moved onto the owner (`o.subscribe(…)`) goes: that name only, and
 * only when nothing else in the file still uses it. */
function dropMovedImport(ctx, sf, alias, calls, out) {
  const spec = alias?.declarations?.find((d) => ts.isImportSpecifier(d) && d.getSourceFile() === sf);
  if (spec === undefined)
    return;
  let uses = 0;
  const visit = (n) => {
    if (ts.isImportDeclaration(n))
      return;
    if (ts.isIdentifier(n) && n.text === spec.name.text) {
      const p = n.parent;
      const symbol = ts.isShorthandPropertyAssignment(p) ? ctx.checker.getShorthandAssignmentValueSymbol(p) :
        ts.isExportSpecifier(p) ? ctx.checker.getExportSpecifierLocalTargetSymbol(p) : ctx.checker.getSymbolAtLocation(n);
      if (symbol === alias)
        uses++;
    }
    ts.forEachChild(n, visit);
  };
  visit(sf);
  if (uses > calls)
    return;
  const named = spec.parent;
  const clause = named.parent;
  if (named.elements.length > 1)
    dropItem(sf, named.elements, spec, out, 'O3');
  else if (clause.name !== undefined)
    out.edit('O3', clause.name.end, named.end, '');
  else {
    const st = clause.parent;
    out.edit('O3', st.getStart(sf), st.end + (/^[ \t]*\r?\n/.exec(sf.text.slice(st.end))?.[0].length ?? 0), '');
    out.removed.add(st);
  }
}

/** Answers true where the call no longer names the function: the owner became its receiver. */
function positionalOwner(ctx, sf, call, rule, args, out) {
  const owner = args[rule.owner];
  const o = owner.getText(sf);
  const last = args[args.length - 1];
  if (rule.form === 'receiver') {
    if (rule.owner !== args.length - 1 || nullable(ctx, owner)) {
      out.hand('O3', call.getStart(sf), `\`${rule.name}(…, owner)\`: call it on the owner (\`owner.${rule.name}(…)\`)`);
      return;
    }
    const callee = call.expression;
    out.edit('O3', callee.getStart(sf), callee.end, `${SIMPLE(owner) ? o : `(${o})`}.${rule.name}`);
    out.edit('O3', args[rule.owner - 1].end, owner.end, '');
    return true;
  }
  else if (rule.form === 'option') {
    if (owner === last) {
      out.edit('O3', owner.getStart(sf), owner.getStart(sf), '{owner: ');
      out.edit('O3', owner.end, owner.end, '}');
    }
    else {
      out.edit('O3', owner.getStart(sf), args[rule.owner + 1].getStart(sf), '');
      out.edit('O3', last.end, last.end, `, {owner: ${o}}`);
    }
  }
  else {
    const rest = args.slice(rule.owner + 1);
    const also = rule.also ?? [];
    if (also.length > 0) {
      const named = rest.slice(rest.length - Math.min(also.length, Math.max(0, rest.length - 2)));
      const kept = rest.slice(0, rest.length - named.length);
      const props = [`owner: ${o}`, ...named.map((a, i) => `${also[i]}: ${a.getText(sf)}`)];
      out.edit('O3', owner.getStart(sf), rest[0].getStart(sf), '');
      out.edit('O3', kept[kept.length - 1].end, last.end, `, {${props.join(', ')}}`);
      return;
    }
    const options = rest[0];
    if (options === undefined)
      out.edit('O3', args[rule.owner - 1].end, owner.end, `, {owner: ${o}}`);
    else if (ts.isObjectLiteralExpression(options)) {
      out.edit('O3', owner.getStart(sf), options.getStart(sf), '');
      const open = options.getStart(sf) + 1;
      out.edit('O3', open, open, options.properties.length > 0 ? `owner: ${o}, ` : `owner: ${o}`);
    }
    else {
      out.edit('O3', owner.getStart(sf), options.getStart(sf), `{owner: ${o}, ...`);
      out.edit('O3', options.end, options.end, '}');
    }
  }
}

/** A v0 caption read or written off its shape (`action.name`) → the role's member (`action.text`): only where the
 * shape has no such member yet — in v1 an action's `name` is its automation id, so the old read would still compile. */
function captionRead(ctx, sf, node, out) {
  const key = node.name.text;
  const type = ctx.checker.getTypeAtLocation(node.expression);
  if (type.flags & ts.TypeFlags.Any) {
    if (ctx.loose && /\.[mc]?js$/.test(sf.fileName) && (key === 'name' || key === 'label') &&
      overActions(ctx, node.expression))
      out.hand('A1', node.getStart(sf), `an untyped \`.${key}\` over what reads as a list of actions: in v1 a caption is ` +
        '`text` (an action\'s `name` is its automation id) — read `.text` if these are u2 actions');
    return;
  }
  for (const part of type.isUnion() ? type.types : [type]) {
    const apparent = ctx.checker.getApparentType(part);
    for (const [holder, renames] of Object.entries(CAPTION_RENAMES)) {
      const to = renames[key];
      if (to === undefined || apparent.getProperty(to) !== undefined)
        continue;
      if (apparent.getProperty(key)?.declarations?.some((d) => ctx.u2Member(d) && ctx.holder(d) === holder)) {
        out.edit('A1', node.name.getStart(sf), node.name.end, to);
        return;
      }
    }
  }
}

const ARRAY_VISITS = new Set(['map', 'flatMap', 'filter', 'find', 'findIndex', 'some', 'every', 'forEach']);

/** `id` is the item of a callback over a list named for actions (`row.actions.map((a) => …)`, `actionsFor(row)`). */
function overActions(ctx, id) {
  const decl = ts.isIdentifier(id) ? ctx.checker.getSymbolAtLocation(id)?.valueDeclaration : undefined;
  const fn = decl !== undefined && ts.isParameter(decl) && decl.parent.parameters[0] === decl ? decl.parent : undefined;
  const call = fn?.parent;
  return call !== undefined && ts.isCallExpression(call) && call.arguments[0] === fn &&
    ts.isPropertyAccessExpression(call.expression) && ARRAY_VISITS.has(call.expression.name.text) &&
    /action/i.test(call.expression.expression.getText());
}

/** Untyped JS literals shaped like a v0 action (`name` + `run`), designer verb (`name` + `produce`) or menu bar
 * item (`label` + `build`). */
const LOOSE_CAPTIONS = [{key: 'name', with: 'run', to: 'text'}, {key: 'name', with: 'produce', to: 'text'},
  {key: 'label', with: 'build', to: 'text'}];

function caption(ctx, sf, node, out) {
  const key = node.name.getText(sf);
  const type = ctx.checker.getContextualType(node.parent);
  const keys = node.parent.properties.map((p) => p.name?.getText(sf));
  if ((type === undefined || type.flags & ts.TypeFlags.Any) && ctx.loose && /\.[mc]?js$/.test(sf.fileName)) {
    const shape = LOOSE_CAPTIONS.find((s) => s.key === key && keys.includes(s.with) && !keys.includes(s.to));
    if (shape !== undefined)
      out.edit('A1', node.name.getStart(sf), node.name.end, ts.isShorthandPropertyAssignment(node) ?
        `${shape.to}: ${key}` : shape.to, true);
    return;
  }
  const shapes = [type, mappedItemType(ctx, node.parent)]
    .filter((t) => t !== undefined && !(t.flags & ts.TypeFlags.Any));
  for (const part of shapes.flatMap((t) => t.isUnion() ? t.types : [t])) {
    const apparent = ctx.checker.getApparentType(part);
    for (const [holder, renames] of Object.entries(CAPTION_RENAMES)) {
      const to = renames[key];
      if (to === undefined)
        continue;
      const decl = [key, to].map((n) => apparent.getProperty(n)?.declarations?.find((d) => ctx.u2Member(d) &&
        ctx.holder(d) === holder)).find((d) => d !== undefined);
      if (decl === undefined || node.parent.properties.some((p) => p.name?.getText(sf) === to))
        continue;
      out.edit('A1', node.name.getStart(sf), node.name.end, ts.isShorthandPropertyAssignment(node) ? `${to}: ${key}` : to);
      return;
    }
  }
}

/** The item a literal returned from a `.map` callback must be, read off the array the `.map(…)` call is passed as:
 * an unannotated callback gives its literal only the shape it infers. */
function mappedItemType(ctx, literal) {
  let at = literal;
  while (ts.isParenthesizedExpression(at.parent))
    at = at.parent;
  const fn = ts.isReturnStatement(at.parent) ? ts.findAncestor(at.parent, ts.isFunctionLike) :
    ts.isArrowFunction(at.parent) && at.parent.body === at ? at.parent : undefined;
  const call = fn?.parent;
  if (call === undefined || !ts.isCallExpression(call) || call.arguments[0] !== fn || !isAccess(call.expression, 'map'))
    return undefined;
  const wanted = ctx.checker.getContextualType(call);
  const array = (wanted === undefined ? [] : wanted.isUnion() ? wanted.types : [wanted])
    .find((t) => ctx.checker.isArrayType(t));
  return array && ctx.checker.getTypeArguments(array)[0];
}

/** The wrappers over `ViewHost` that v1 removes. */
const VIEW_WRAPPERS = new Set(['tableViewChrome', 'showInView', 'viewSignals', 'ownView', 'appView']);

/** The import of u2 that `id` stands for, and the name it imports — read off the import itself, so code compiled
 * against a u2 that no longer exports the name is read too. */
function importedFromU2(ctx, id) {
  const symbol = ctx.checker.getSymbolAtLocation(id);
  const spec = symbol?.declarations?.find((d) => ts.isImportSpecifier(d));
  if (spec === undefined)
    return null;
  const from = spec.parent.parent.parent.moduleSpecifier;
  const aliased = symbol.flags & ts.SymbolFlags.Alias ? ctx.checker.getAliasedSymbol(symbol) : undefined;
  const u2 = ts.isStringLiteral(from) && from.text.startsWith('@datagrok-libraries/u2') ||
    aliased?.declarations?.some((d) => ctx.inU2(d));
  return u2 ? {name: (spec.propertyName ?? spec.name).text, spec} : null;
}

/** One call of a view wrapper as `ViewHost` says it; false where a person has to (reported). */
function viewWrapper(sf, call, name, out) {
  const args = call.arguments;
  const at = call.getStart(sf);
  const [first, second] = args;
  if (name === 'tableViewChrome' && (args.length === 1 || args.length === 2)) {
    if (second !== undefined && !(ts.isObjectLiteralExpression(second) &&
      !second.properties.some((p) => p.name?.getText(sf) === 'owner'))) {
      out.hand('W1', at, '`tableViewChrome(v, {owner})`: `ViewHost.of(v).ribbon(options)` is owned by the view — drop ' +
        'the owner, or dispose the chrome with it by hand (`owner.own(() => chrome.dispose())`)');
      return false;
    }
    out.edit('W1', call.expression.getStart(sf), first.getStart(sf), 'ViewHost.of(');
    out.edit('W1', first.end, second === undefined ? first.end : second.getStart(sf), ').ribbon(');
    return true;
  }
  if (name === 'showInView' && (args.length === 2 || args.length === 3)) {
    out.edit('W1', at, second.end, `ViewHost.of(${second.getText(sf)}).show(${first.getText(sf)}`);
    return true;
  }
  if (name === 'viewSignals' && args.length === 1) {
    out.edit('W1', at, call.end, `ViewHost.of(${first.getText(sf)}).signals`);
    return true;
  }
  out.hand('W1', at, {
    viewSignals: '`viewSignals(v, owner)`: `ViewHost.of(v).signals`, owned by the view',
    ownView: '`ownView(view, control)`: build the control under the view\'s host, `ViewHost.of(view).build(() => …)`',
    appView: '`appView({name, content, …})`: `ViewHost.create({name, content: () => ({root, ribbon, status, path})}).view`, ' +
      'the content built under the view\'s host',
  }[name] ?? `\`${name}(…)\`: say it through \`ViewHost\``);
  return false;
}

function viewRules(ctx, sf, out) {
  if (!ctx.rules.has('W1'))
    return;
  const uses = new Map();
  const visit = (node) => {
    if (ts.isIdentifier(node) && !ts.isImportSpecifier(node.parent)) {
      const imported = importedFromU2(ctx, node);
      if (imported !== null && VIEW_WRAPPERS.has(imported.name)) {
        const use = uses.get(imported.spec) ?? {total: 0, rewritten: 0};
        uses.set(imported.spec, use);
        use.total++;
        const call = node.parent;
        if (ts.isCallExpression(call) && call.expression === node && viewWrapper(sf, call, imported.name, out))
          use.rewritten++;
      }
    }
    ts.forEachChild(node, visit);
  };
  visit(sf);
  const rewritten = [...uses].filter(([, use]) => use.rewritten > 0);
  if (rewritten.length === 0)
    return;
  const imported = sf.statements.some((s) => ts.isImportDeclaration(s) && s.importClause?.namedBindings &&
    ts.isNamedImports(s.importClause.namedBindings) && s.importClause.namedBindings.elements.some((e) => e.name.text === 'ViewHost'));
  let added = imported;
  for (const [spec, use] of uses) {
    if (use.rewritten !== use.total)
      continue;
    if (!added) {
      out.edit('W1', spec.getStart(sf), spec.end, 'ViewHost');
      added = true;
    }
    else
      dropItem(sf, spec.parent.elements, spec, out, 'W1');
  }
  if (!added) {
    const elements = rewritten[0][0].parent.elements;
    out.edit('W1', elements[elements.length - 1].end, elements[elements.length - 1].end, ', ViewHost');
  }
}

/** The entry point an import (or a re-export) from a path inside u2 can take instead: the first entry that exports
 * every name it takes, as the very symbols it takes; a string to say why where none does. */
function entryOf(ctx, st, subpath) {
  const clause = ts.isImportDeclaration(st) ? st.importClause : undefined;
  const bindings = ts.isImportDeclaration(st) ? clause?.namedBindings : st.exportClause;
  if (ts.isImportDeclaration(st) && clause === undefined)
    return 'a side-effect import of a file inside u2: the entry points load what they need — drop it';
  if (clause?.name !== undefined || bindings === undefined || !ts.isNamedImports(bindings) && !ts.isNamedExports(bindings))
    return 'a default, `* as` or `export *` import of a file inside u2: name what it takes and import that from an entry';
  const names = bindings.elements.map((e) => {
    const symbol = ts.isExportSpecifier(e) ? ctx.checker.getExportSpecifierLocalTargetSymbol(e) :
      ctx.checker.getSymbolAtLocation(e.name);
    const target = symbol && ctx.target(symbol);
    return {name: (e.propertyName ?? e.name).text, target: target?.declarations?.length ? target : undefined};
  });
  const first = entryFor(subpath);
  for (const entry of [first, ...ENTRIES.filter((e) => e !== first)]) {
    const exported = ctx.entryExports(entry);
    if (exported && names.every((n) => exported.has(n.name) && (n.target === undefined || exported.get(n.name) === n.target)))
      return {entry};
  }
  const missing = names.filter((n) => !ENTRIES.some((e) => ctx.entryExports(e)?.has(n.name))).map((n) => `\`${n.name}\``);
  return missing.length > 0 ? `${missing.join(', ')} ${missing.length > 1 ? 'are' : 'is'} on no entry point: ask u2 ` +
    'to export it, or do without' : 'its names are on different entry points: split the import, one per entry';
}

function importRules(ctx, sf, out) {
  if (!ctx.rules.has('I1'))
    return;
  for (const st of sf.statements) {
    const spec = ts.isImportDeclaration(st) || ts.isExportDeclaration(st) ? st.moduleSpecifier : undefined;
    const m = spec !== undefined && ts.isStringLiteral(spec) ? U2_SPECIFIER.exec(spec.text) : null;
    const subpath = m?.[1] ?? '';
    if (m === null || ctx.entryPaths.has(subpath) || /^\/css\//.test(subpath) || out.removed.has(st))
      continue;
    const to = entryOf(ctx, st, subpath);
    if (typeof to === 'string')
      out.hand('I1', st.getStart(sf), `\`${spec.text}\`: ${to}`);
    else
      out.edit('I1', spec.getStart(sf) + 1, spec.end - 1, to.entry);
  }
  const visit = (node) => {
    if (ts.isCallExpression(node) && node.expression.kind === ts.SyntaxKind.ImportKeyword &&
      node.arguments[0] !== undefined && ts.isStringLiteralLike(node.arguments[0])) {
      const m = U2_SPECIFIER.exec(node.arguments[0].text);
      if (m !== null && !ctx.entryPaths.has(m[1] ?? '') && !/^\/css\//.test(m[1] ?? ''))
        out.hand('I1', node.getStart(sf), `\`import('${node.arguments[0].text}')\`: load an entry point instead`);
    }
    ts.forEachChild(node, visit);
  };
  visit(sf);
}

/** The handlers whose answer v1 types as a promise or nothing (`boolean` too, for the keys): key → the u2 shapes. */
const HANDLERS = {run: new Set(['Action', 'RibbonAction', 'RibbonButton']),
  onEnter: new Set(['InputOptions']), onEscape: new Set(['InputOptions'])};

const typeIs = (type, flags) => (type.isUnion() ? type.types : [type]).every((t) => (t.flags & flags) !== 0);

function handlerRules(ctx, sf, out) {
  // a type contract: JavaScript has none to keep
  if (!ctx.rules.has('R1') || /\.[mc]?js$/.test(sf.fileName))
    return;
  const visit = (node) => {
    if (ts.isPropertyAssignment(node) && ts.isObjectLiteralExpression(node.parent) &&
      Object.hasOwn(HANDLERS, node.name.getText(sf)))
      handler(ctx, sf, node, out);
    ts.forEachChild(node, visit);
  };
  visit(sf);
}

function handler(ctx, sf, node, out) {
  const key = node.name.getText(sf);
  const shape = ctx.checker.getContextualType(node.parent);
  if (shape === undefined || shape.flags & ts.TypeFlags.Any)
    return;
  const owned = (shape.isUnion() ? shape.types : [shape]).some((t) => ctx.checker.getApparentType(t).getProperty(key)
    ?.declarations?.some((d) => ctx.u2Member(d) && HANDLERS[key].has(ctx.holder(d))));
  const fn = node.initializer;
  if (!owned || !ts.isArrowFunction(fn) && !ts.isFunctionExpression(fn))
    return;
  const nothing = ts.TypeFlags.Void | ts.TypeFlags.Undefined | ts.TypeFlags.Never | ts.TypeFlags.Any;
  // a promise answers nothing in the end; a plain answer is nothing — or, for a key, whether it was taken
  const fine = (type, promise) => typeIs(type, promise || key === 'run' ? nothing : nothing | ts.TypeFlags.BooleanLike);
  const async = fn.modifiers?.some((m) => m.kind === ts.SyntaxKind.AsyncKeyword) === true;
  const thenable = (type) => type.getProperty('then') !== undefined && !(type.flags & ts.TypeFlags.Any);
  if (!ts.isArrowFunction(fn) || ts.isBlock(fn.body)) {
    const signature = ctx.checker.getSignatureFromDeclaration(fn);
    const answer = signature && ctx.checker.getReturnTypeOfSignature(signature);
    const promise = answer !== undefined && (async || thenable(answer));
    if (answer !== undefined && !fine(promise ? ctx.checker.getAwaitedType(answer) ?? answer : answer, promise))
      out.hand('R1', node.getStart(sf), `\`${key}\` returns a value: return nothing${key === 'run' ? '' : ' or a boolean'}` +
        ' — or a promise of nothing');
    return;
  }
  const body = fn.body;
  const type = ctx.checker.getTypeAtLocation(body);
  const promise = async || thenable(type);
  if (fine(promise ? ctx.checker.getAwaitedType(type) ?? type : type, promise))
    return;
  if (promise) {
    if (!async)
      out.edit('R1', fn.getStart(sf), fn.getStart(sf), 'async ');
    out.edit('R1', body.getStart(sf), body.getStart(sf), thenable(type) ? '{ await ' : '{ ');
    out.edit('R1', body.end, body.end, '; }');
    return;
  }
  const bare = SIMPLE(body) || ts.isNewExpression(body) || ts.isLiteralExpression(body);
  out.edit('R1', body.getStart(sf), body.getStart(sf), bare ? 'void ' : 'void (');
  if (!bare)
    out.edit('R1', body.end, body.end, ')');
}

/** In a spec document: a u2 node's `enabled` prop → `disabled` (a reason), a bound one renamed and reported. */
function specRules(sf, out) {
  const visit = (node) => {
    if (ts.isObjectLiteralExpression(node)) {
      const tag = node.properties.find((p) => ts.isPropertyAssignment(p) && p.name.getText(sf).replace(/['"]/g, '') === 'tag');
      if (tag && ts.isStringLiteral(tag.initializer) && tag.initializer.text.startsWith('u2-')) {
        for (const p of node.properties) {
          const key = ts.isPropertyAssignment(p) ? p.name.getText(sf).replace(/['"]/g, '') : '';
          if ((key === 'props' || key === 'bind') && ts.isObjectLiteralExpression(p.initializer))
            specProps(sf, p.initializer, key, out);
        }
      }
    }
    ts.forEachChild(node, visit);
  };
  visit(sf);
}

function specProps(sf, props, key, out) {
  const enabled = props.properties.find((p) => ts.isPropertyAssignment(p) && p.name.getText(sf).replace(/['"]/g, '') ===
    'enabled');
  if (enabled === undefined)
    return;
  const name = enabled.name;
  const quote = name.getText(sf)[0] === '\'' ? '\'' : '"';
  const value = enabled.initializer;
  if (key === 'props' && value.kind === ts.SyntaxKind.TrueKeyword) {
    dropItem(sf, props.properties, enabled, out, 'S1');
    return;
  }
  out.edit('S1', name.getStart(sf), name.end, `${quote}disabled${quote}`);
  if (key === 'props' && value.kind === ts.SyntaxKind.FalseKeyword) {
    out.edit('S1', value.getStart(sf), value.end, `${quote}${REASON_NEEDED}${quote}`);
    out.hand('S1', value.getStart(sf), `write the reason the user reads in place of "${REASON_NEEDED}"`);
  }
  else
    out.hand('S1', enabled.getStart(sf), key === 'bind' ? 'a bound `disabled` reads a reason (string | null): bind ' +
      'it to one' : '`disabled` takes the reason the user reads (or null), not a boolean');
}

const PASSES = [valueRules, reasonRules, ownerRules, viewRules, importRules, handlerRules, importsFor];

/** The edits and hand-finish items one file collects; overlapping edits keep the first and report the rest. */
class FileResult {
  constructor(file, offset = 0, text = undefined) {
    this.file = file;
    this.offset = offset;
    this.text = text;
    this.edits = [];
    this.hands = [];
    this.needed = new Map();
    this.removed = new Set();
  }

  edit(rule, start, end, text, loose = false) {
    this.edits.push({rule, start: start + this.offset, end: end + this.offset, text, loose, seq: this.edits.length});
  }

  hand(rule, pos, message) {
    this.hands.push({rule, pos: pos + this.offset, message});
  }

  /** A name the rewrites use that the file has to import. */
  needs(name, pos) {
    if (!this.needed.has(name))
      this.needed.set(name, pos);
  }
}

/** Applies span edits back to front. An insertion (empty span) never blocks an edit starting where it does; edits
 * at one place keep the order they were made in; an edit overlapping a kept one is dropped and reported. */
function applyEdits(text, edits) {
  const sorted = [...edits].sort((a, b) => a.start - b.start || (a.end - a.start) - (b.end - b.start) || a.seq - b.seq);
  const kept = [];
  const dropped = [];
  let last = -1;
  for (const e of sorted) {
    if (e.start < last)
      dropped.push(e);
    else {
      kept.push(e);
      last = Math.max(last, e.end);
    }
  }
  let result = text;
  for (const e of [...kept].sort((a, b) => b.start - a.start || (b.end - b.start) - (a.end - a.start) || b.seq - a.seq))
    result = result.slice(0, e.start) + e.text + result.slice(e.end);
  return {result, kept: kept.length, dropped};
}

/** The changed lines, `-` before and `+` after, for a dry run a person reads. */
function lineDiff(before, after) {
  const a = before.split('\n');
  const b = after.split('\n');
  let head = 0;
  while (head < a.length && head < b.length && a[head] === b[head])
    head++;
  const out = [];
  let i = head;
  let j = head;
  while (i < a.length || j < b.length) {
    if (i < a.length && j < b.length && a[i] === b[j]) {
      i++;
      j++;
      continue;
    }
    const next = b.indexOf(a[i], j);
    if (i < a.length && (next < 0 || next - j > 3)) {
      out.push(`${i + 1}- ${a[i].trim()}`);
      i++;
      if (j < b.length && (i >= a.length || b.indexOf(a[i], j) !== j)) {
        out.push(`${j + 1}+ ${b[j].trim()}`);
        j++;
      }
    }
    else {
      out.push(`${j + 1}+ ${b[j].trim()}`);
      j++;
    }
  }
  return out;
}

function lineOf(text, pos) {
  let line = 1;
  for (let i = 0; i < pos && i < text.length; i++)
    if (text.charCodeAt(i) === 10)
      line++;
  return line;
}

export function run(options) {
  const files = options.paths.flatMap((p) => listFiles(p))
    .filter((f) => !f.endsWith('.json') || holdsSpec(readFileSync(f, 'utf8')));
  const {program, blocksOf, entries} = createProgram({project: options.project, u2Src: options.u2Src, files});
  const ctx = new Context(program, options, entries);
  const targets = files.length > 0 ? files : program.getSourceFiles().map((sf) => sf.fileName)
    .filter((f) => !f.includes('/node_modules/') && !/\.d\.[mc]?ts$/.test(f));
  const report = {files: {}, totals: {}, hand: [], diffs: {}};
  let unresolved = 0;
  for (const file of targets) {
    const text = readFileSync(file, 'utf8');
    const parts = file.endsWith('.md') ?
      (blocksOf.get(file) ?? []).map((b) => ({sf: program.getSourceFile(b.name), offset: b.start})) :
      file.endsWith('.json') ? [] : [{sf: program.getSourceFile(file), offset: 0}];
    const result = new FileResult(file);
    for (const {sf, offset} of parts) {
      if (sf === undefined)
        continue;
      result.offset = offset;
      for (const pass of PASSES)
        pass(ctx, sf, result);
    }
    if (options.rules.has('S1')) {
      const docs = file.endsWith('.json') ? [{name: file, start: 0, text}] :
        file.endsWith('.md') ? markdownBlocks(file, text, 'json', 'json') : [];
      for (const doc of docs.filter((d) => holdsSpec(d.text))) {
        result.offset = doc.start;
        specRules(ts.parseJsonText(doc.name, doc.text), result);
      }
    }
    for (const {sf} of parts) {
      for (const st of sf?.statements ?? []) {
        const spec = (ts.isImportDeclaration(st) || ts.isExportDeclaration(st)) ? st.moduleSpecifier : undefined;
        if (spec && ts.isStringLiteral(spec) && spec.text.startsWith('@datagrok-libraries/u2') &&
          ctx.checker.getSymbolAtLocation(spec) === undefined)
          unresolved++;
      }
    }
    const {result: next, kept, dropped} = applyEdits(text, result.edits);
    for (const e of dropped)
      result.hands.push({rule: e.rule, pos: e.start, message: 'overlaps another rewrite: finish by hand'});
    const rel = relative(options.root ?? process.cwd(), file).replace(/\\/g, '/');
    const counts = {};
    for (const e of result.edits)
      counts[e.rule] = (counts[e.rule] ?? 0) + 1;
    for (const h of result.hands)
      report.hand.push({file: rel, line: lineOf(text, h.pos), rule: h.rule, message: h.message});
    if (kept > 0 || result.hands.length > 0) {
      report.files[rel] = {edits: kept, rules: counts, loose: result.edits.filter((e) => e.loose).length,
        hand: result.hands.length, sites: result.edits.map((e) => `${e.rule}:${lineOf(text, e.start)}${e.loose ? '?' : ''}`)};
    }
    for (const [rule, n] of Object.entries(counts))
      report.totals[rule] = (report.totals[rule] ?? 0) + n;
    if (options.write && next !== text)
      writeFileSync(file, next);
    if (options.texts && next !== text)
      (report.texts ??= {})[file] = next;
    if (options.diff && next !== text)
      report.diffs[rel] = lineDiff(text, next);
  }
  report.edits = Object.values(report.files).reduce((n, f) => n + f.edits, 0);
  report.handCount = report.hand.length;
  report.unresolvedU2Imports = unresolved;
  return report;
}

function parseArgs(argv) {
  const options = {paths: [], write: false, loose: false, json: false, rules: new Set(Object.keys(RULES))};
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--project')
      options.project = argv[++i];
    else if (a === '--u2-src')
      options.u2Src = argv[++i];
    else if (a === '--report')
      options.report = argv[++i];
    else if (a === '--root')
      options.root = argv[++i];
    else if (a === '--rules')
      options.rules = new Set(argv[++i].split(','));
    else if (a === '--write')
      options.write = true;
    else if (a === '--loose')
      options.loose = true;
    else if (a === '--json')
      options.json = true;
    else if (a === '--diff')
      options.diff = true;
    else if (a.startsWith('--'))
      throw new Error(`unknown option ${a}`);
    else
      options.paths.push(a);
  }
  if (!options.project && options.paths.length === 0)
    throw new Error('give --project <tsconfig.json> and/or the files or directories to migrate');
  for (const r of options.rules)
    if (!RULES[r])
      throw new Error(`unknown rule ${r}; rules: ${Object.keys(RULES).join(', ')}`);
  return options;
}

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const options = parseArgs(process.argv.slice(2));
  const report = run(options);
  if (options.report)
    writeFileSync(options.report, JSON.stringify(report, null, 2) + '\n');
  if (options.json)
    console.log(JSON.stringify(report, null, 2));
  else {
    for (const [file, f] of Object.entries(report.files))
      console.log(`${file}: ${f.edits} edit(s) ${JSON.stringify(f.rules)}${f.loose ? `, ${f.loose} untyped` : ''}${f.hand ? `, ${f.hand} by hand` : ''}`);
    for (const h of report.hand)
      console.log(`  hand ${h.file}:${h.line} ${h.rule}: ${h.message}`);
    for (const [file, lines] of Object.entries(report.diffs))
      console.log(`--- ${file}\n${lines.join('\n')}`);
    console.log(`${options.write ? 'wrote' : 'would write'} ${report.edits} edit(s) in ${Object.keys(report.files).length} ` +
      `file(s): ${JSON.stringify(report.totals)}; ${report.handCount} to finish by hand` +
      (report.unresolvedU2Imports ? `; ${report.unresolvedU2Imports} u2 import(s) did not resolve — build u2 first` : ''));
  }
}
