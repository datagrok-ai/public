/// TypeScript source files and their declarations (build-plan.md WO-3b): one compiler-API parse per file
/// (`ts.createSourceFile`, no Program) shared by the ts-declarations, ts-imports and ts-uses extractors;
/// `source-file` and `declaration` nodes, `declares` file->declaration and class->member, `extends` and
/// `implements` resolved by name in the same file, through the imports, then the same package or library.
/// A declaration is a node only when it names a type or belongs to the JS API surface (`isNode`).
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {ts} from 'ts-morph';
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../context';
import {pkgId, libId, fileId, declId, countLines, sourceFileRow} from '../../../ids';
import {listPackages, listLibraries, entryPoints} from './packages';

const SOURCE_IGNORE = ['**/node_modules/**', '**/dist/**', '**/*.d.ts', '**/__tests__/fixtures/**'];
const JS_API = libId('js-api');
/** Where a unit's TypeScript is: `src/` for packages and libraries; the CLI keeps its sources beside their Babel output under `bin/`. */
const SOURCE_GLOBS: Record<string, string> = {[libId('tools')]: 'bin/**/*.ts'};
const DEFAULT_SOURCES = 'src/**/*.{ts,tsx}';
/** The test folders of a library, the JS API and the CLI (u2's node:test suites, the JS API's scripts/unit), fixtures left out. */
const LIB_TEST_SOURCES = '{tests,test,__tests__,scripts/unit}/**/*.{test,spec}.{ts,js,mjs,cjs}';
const TEST_IGNORE = [...SOURCE_IGNORE, '**/fixtures/**'];
const JS_API_ROOTS = ['ui.ts', 'grok.ts', 'dg.ts'];
const SIGNATURE_CAP = 160;
/** Files an import may name that are not TypeScript sources: present on disk, absent from the graph. */
const ASSET_EXTENSIONS = ['.js', '.mjs', '.cjs', '.json', '.vue', '.css', '.wasm', '/index.js'];
/** Declarations a class or interface may extend or implement. */
const TYPE_KINDS = new Set(['class', 'interface']);
/** Declaration kinds that name a type, and so are nodes wherever they are declared. */
const NODE_KINDS = new Set(['class', 'interface', 'enum', 'type', 'mixin']);
const PRIVATE = ts.ModifierFlags.Private | ts.ModifierFlags.Protected;

export interface TsDecl {
  /** `Name`, `Owner.member` or `Namespace.member`. */
  name: string;
  kind: string;
  accessor?: 'get' | 'set';
  /** Name of the lexical container: the class of a member, the namespace of what it holds. */
  container?: string;
  exported: boolean;
  line: number;
  signature: string;
  documented: boolean;
  deprecated: boolean;
  extends: string[];
  implements: string[];
  /** Leftmost identifier of a variable initializer (`new Shell()` -> Shell), which is the class behind a `grok.*` namespace. */
  alias?: string;
}

export interface TsImport {
  specifier: string;
  /** What the statement names: imported names, `default`, `*`. */
  symbols: string[];
  /** Local name -> imported name. */
  locals: Record<string, string>;
  /** The local name of `* as X`. */
  namespace?: string;
  /** The exported name of `export * as X from`. */
  alias?: string;
  /** An `export ... from` statement: the file re-exports what it names. */
  reexport?: boolean;
}

/** A package, a library or the JS API: the owner of the files under its folder. */
export interface TsUnit {
  id: string;
  dir: string;
  npm?: string;
  /** A package's entry points (`entryPoints`), flagged `entry` on their file rows. */
  entries?: Set<string>;
}

export interface TsFile {
  /** Posix path from the monorepo root. */
  path: string;
  unit: TsUnit;
  loc: number;
  generated: boolean;
  decls: TsDecl[];
  imports: TsImport[];
}

export interface ApiEntry {
  file: TsFile;
  decl: TsDecl;
  id: string;
}

export interface ImportTarget {
  to?: string;
  /** A specifier that should have resolved (relative or Datagrok-scoped) but did not; bare npm names and asset files are neither. */
  unresolved?: boolean;
}

const CACHE = new WeakMap<BuildContext, TsSources>();

/** A declaration is a node when it names a type, or when it is part of the JS API surface a plugin can call;
 * a unit's own members, consts and helper functions are represented by their file and their type. */
export function isNode(file: TsFile, decl: TsDecl): boolean {
  return NODE_KINDS.has(decl.kind) || (file.unit.id === JS_API && decl.exported);
}

/** The parsed sources of this build, parsed once for the three extractors; problems are reported by whichever runs first. */
export function tsSources(ctx: BuildContext, emitter: Emitter): TsSources {
  let sources = CACHE.get(ctx);
  if (!sources) CACHE.set(ctx, sources = new TsSources(ctx, emitter));
  return sources;
}

export class TsSources {
  readonly files: TsFile[] = [];
  readonly byPath = new Map<string, TsFile>();
  readonly units = new Map<string, TsUnit>();
  /** Files that could not be read or parsed cleanly. */
  failed = 0;
  private byNpm = new Map<string, TsUnit>();
  private api?: Map<string, ApiEntry[]>;
  private unitTypes = new Map<string, Map<string, string>>();

  constructor(private ctx: BuildContext, emitter: Emitter) {
    const jsApi = listLibraries(ctx.repoRoot).find((l) => l.folder === 'js-api');
    const units: TsUnit[] = jsApi ? [{id: JS_API, dir: jsApi.dir, npm: jsApi.json.name}] : [];
    if (ctx.mode === 'full') {
      for (const p of listPackages(ctx.repoRoot)) units.push({id: pkgId(p.folder), dir: p.dir, npm: p.json.name, entries: new Set(entryPoints(p))});
      for (const l of listLibraries(ctx.repoRoot)) if (l.folder !== 'js-api') units.push({id: libId(l.folder), dir: l.dir, npm: l.json.name});
    }
    for (const unit of units) {
      this.units.set(unit.id, unit);
      if (unit.npm) this.byNpm.set(unit.npm, unit);
      const files = globSync(`${unit.dir}/${SOURCE_GLOBS[unit.id] ?? DEFAULT_SOURCES}`, {cwd: ctx.repoRoot, ignore: SOURCE_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true});
      if (unit.id === JS_API) files.push(...JS_API_ROOTS.map((f) => `${unit.dir}/${f}`).filter((f) => fs.existsSync(path.join(ctx.repoRoot, f))));
      if (unit.id.startsWith('lib:')) files.push(...globSync(`${unit.dir}/${LIB_TEST_SOURCES}`, {cwd: ctx.repoRoot, ignore: TEST_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}));
      for (const file of [...new Set(files)].sort()) {
        const parsed = this.parse(file, unit, emitter);
        this.files.push(parsed);
        this.byPath.set(file, parsed);
      }
    }
  }

  declId(file: TsFile, decl: TsDecl): string {
    return declId(file.path, decl.name, decl.accessor);
  }

  /** The node an import specifier names (build-plan.md WO-3b, ts_imports.py rules). */
  resolveImport(file: TsFile, specifier: string): ImportTarget {
    if (specifier.startsWith('.')) {
      const base = path.posix.normalize(path.posix.join(path.posix.dirname(file.path), specifier));
      const hit = this.resolveFile(base);
      if (hit) return {to: fileId(hit)};
      return ['', ...ASSET_EXTENSIONS].some((ext) => fs.existsSync(path.join(this.ctx.repoRoot, base + ext))) ? {} : {unresolved: true};
    }
    const lib = /^@datagrok-libraries\/([^/]+)(?:\/(.+))?$/.exec(specifier);
    if (lib) {
      const unit = this.units.get(libId(lib[1]));
      if (!unit) return {unresolved: true};
      const hit = lib[2] ? this.resolveFile(`${unit.dir}/${lib[2]}`) : undefined;
      return {to: hit ? fileId(hit) : unit.id};
    }
    const pkg = /^(@datagrok\/[^/]+)(?:\/(.+))?$/.exec(specifier);
    if (pkg) {
      const unit = this.byNpm.get(pkg[1]);
      if (!unit) return {unresolved: true};
      if (!pkg[2]) return {to: unit.id};
      const hit = this.resolveFile(`${unit.dir}/${pkg[2]}`);
      return hit ? {to: fileId(hit)} : {unresolved: true};
    }
    if (specifier === 'datagrok-api' || specifier.startsWith('datagrok-api/')) return {to: JS_API};
    return {};
  }

  /** Exported JS API declarations by dotted name (`DataFrame`, `DataFrame.name`, `input.string`); a getter precedes its setter. */
  apiIndex(): Map<string, ApiEntry[]> {
    if (this.api) return this.api;
    const api = this.api = new Map<string, ApiEntry[]>();
    const add = (name: string, file: TsFile, decl: TsDecl) => {
      let list = api.get(name);
      if (!list) api.set(name, list = []);
      list.push({file, decl, id: this.declId(file, decl)});
    };
    const apiFiles = this.files.filter((f) => f.unit.id === JS_API);
    for (const file of apiFiles)
      for (const decl of file.decls)
        if (decl.exported) add(decl.name, file, decl);
    for (const file of apiFiles)
      for (const imp of file.imports) {
        const target = imp.alias ? this.fileOf(this.resolveImport(file, imp.specifier).to) : undefined;
        for (const decl of target?.decls ?? [])
          if (decl.exported && !decl.container) add(`${imp.alias}.${decl.name}`, target!, decl);
      }
    for (const list of api.values()) list.sort((a, b) => (a.decl.accessor === 'set' ? 1 : 0) - (b.decl.accessor === 'set' ? 1 : 0));
    return api;
  }

  /** The id of the class or interface a heritage name means from inside [scope] (the namespace of the extending declaration), or
   * undefined when nothing in reach declares it. */
  resolveHeritage(file: TsFile, name: string, scope?: string): string | undefined {
    const segments = name.split('.');
    const simple = segments[segments.length - 1];
    if (segments.length === 1) {
      for (let s = scope; s; s = s.includes('.') ? s.slice(0, s.lastIndexOf('.')) : undefined)
        if (this.typeIn(file, `${s}.${simple}`)) return declId(file.path, `${s}.${simple}`);
      if (this.typeIn(file, simple)) return declId(file.path, simple);
      for (const imp of file.imports) {
        const imported = imp.locals[simple];
        if (imported === undefined) continue;
        const hit = this.typeBehind(file, imp.specifier, imported);
        if (hit) return hit;
      }
      return this.unitTypeIndex(file.unit).get(simple);
    }
    const imp = file.imports.find((i) => i.namespace === segments[0]);
    return imp ? this.typeBehind(file, imp.specifier, segments.slice(1).join('.')) : undefined;
  }

  private typeBehind(file: TsFile, specifier: string, name: string): string | undefined {
    const {to} = this.resolveImport(file, specifier);
    if (to === JS_API) return this.apiIndex().get(name)?.find((e) => TYPE_KINDS.has(e.decl.kind))?.id;
    const target = this.fileOf(to);
    return target && this.typeIn(target, name) ? declId(target.path, name) : undefined;
  }

  private fileOf(id: string | undefined): TsFile | undefined {
    return id?.startsWith('file:') ? this.byPath.get(id.slice(5)) : undefined;
  }

  private typeIn(file: TsFile, name: string): boolean {
    return file.decls.some((d) => d.name === name && TYPE_KINDS.has(d.kind));
  }

  private unitTypeIndex(unit: TsUnit): Map<string, string> {
    let index = this.unitTypes.get(unit.id);
    if (index) return index;
    this.unitTypes.set(unit.id, index = new Map());
    for (const file of this.files)
      if (file.unit === unit)
        for (const d of file.decls)
          if (!d.container && TYPE_KINDS.has(d.kind) && !index.has(d.name)) index.set(d.name, declId(file.path, d.name));
    return index;
  }

  private resolveFile(base: string): string | undefined {
    for (const candidate of [base, `${base}.ts`, `${base}.tsx`, base.replace(/\.js$/, '.ts'), `${base}/index.ts`, `${base}/index.tsx`])
      if (this.byPath.has(candidate)) return candidate;
    return undefined;
  }

  private parse(file: string, unit: TsUnit, emitter: Emitter): TsFile {
    const parsed: TsFile = {path: file, unit, loc: 0, generated: file.endsWith('.g.ts'), decls: [], imports: []};
    let text: string;
    try {
      text = fs.readFileSync(path.join(this.ctx.repoRoot, file), 'utf8');
    } catch (e: any) {
      this.failed++;
      emitter.problem('invalid_rows', `${file}: ${e.message}`);
      return parsed;
    }
    parsed.loc = countLines(text);
    const kind = file.endsWith('.tsx') ? ts.ScriptKind.TSX : /\.[cm]?js$/.test(file) ? ts.ScriptKind.JS : ts.ScriptKind.TS;
    const sf = ts.createSourceFile(file, text, ts.ScriptTarget.Latest, true, kind);
    const diagnostics = (sf as any).parseDiagnostics as ts.Diagnostic[] ?? [];
    if (diagnostics.length) {
      this.failed++;
      const d = diagnostics[0];
      emitter.problem('invalid_rows', `${file}:${sf.getLineAndCharacterOfPosition(d.start ?? 0).line + 1}: ${ts.flattenDiagnosticMessageText(d.messageText, ' ')}`);
    }
    new FileWalk(sf, text, parsed).run();
    return parsed;
  }
}

/** One file's statements into declarations and imports; namespaces recurse with their name as the prefix. */
class FileWalk {
  private exportList = new Set<string>();
  private byName = new Map<string, TsDecl>();

  constructor(private sf: ts.SourceFile, private text: string, private file: TsFile) {}

  run(): void {
    for (const st of this.sf.statements)
      if (ts.isExportDeclaration(st) && !st.moduleSpecifier && st.exportClause && ts.isNamedExports(st.exportClause))
        for (const e of st.exportClause.elements) this.exportList.add((e.propertyName ?? e.name).text);
    this.visitStatements(this.sf.statements, '', true);
    ts.forEachChild(this.sf, (node) => this.visitDynamicImports(node));
  }

  private visitStatements(statements: ts.NodeArray<ts.Statement>, prefix: string, exportedScope: boolean): void {
    const container = prefix ? prefix.slice(0, -1) : undefined;
    for (const st of statements) {
      if (ts.isImportDeclaration(st) || ts.isExportDeclaration(st) || ts.isImportEqualsDeclaration(st)) {
        this.importOf(st);
        continue;
      }
      const exported = exportedScope && this.isExported(st);
      if (ts.isClassDeclaration(st) || ts.isInterfaceDeclaration(st)) {
        if (!st.name) continue;
        const name = prefix + st.name.text;
        const decl = this.add(st, name, ts.isClassDeclaration(st) ? 'class' : 'interface', exported, container);
        for (const clause of st.heritageClauses ?? []) {
          const names = clause.types.map((t) => t.expression.getText(this.sf));
          if (clause.token === ts.SyntaxKind.ExtendsKeyword) decl.extends.push(...names);
          else decl.implements.push(...names);
        }
        for (const m of st.members) this.member(m, name, exported);
      }
      else if (ts.isFunctionDeclaration(st)) {
        if (st.name) this.add(st, prefix + st.name.text, 'function', exported, container);
      }
      else if (ts.isEnumDeclaration(st)) this.add(st, prefix + st.name.text, 'enum', exported, container);
      else if (ts.isTypeAliasDeclaration(st)) this.add(st, prefix + st.name.text, 'type', exported, container);
      else if (ts.isVariableStatement(st)) {
        for (const d of st.declarationList.declarations) {
          if (!ts.isIdentifier(d.name)) continue;
          const decl = this.add(st, prefix + d.name.text, 'const', exported, container);
          decl.alias = d.initializer && leftmostIdentifier(d.initializer);
        }
      }
      else if (ts.isModuleDeclaration(st) && ts.isIdentifier(st.name)) {
        const name = prefix + st.name.text;
        this.add(st, name, 'const', exported, container);
        let body = st.body;
        while (body && ts.isModuleDeclaration(body)) body = body.body;
        if (body && ts.isModuleBlock(body)) this.visitStatements(body.statements, `${name}.`, exported);
      }
    }
  }

  private member(m: ts.ClassElement | ts.TypeElement, owner: string, ownerExported: boolean): void {
    if (ts.isConstructorDeclaration(m) || ts.isIndexSignatureDeclaration(m) || !m.name || ts.isPrivateIdentifier(m.name)) return;
    const kind = ts.isMethodDeclaration(m) || ts.isMethodSignature(m) ? 'method' : ts.isGetAccessor(m) ? 'getter' : ts.isSetAccessor(m) ? 'setter' :
      ts.isPropertyDeclaration(m) || ts.isPropertySignature(m) ? 'prop' : undefined;
    if (!kind) return;
    const exported = ownerExported && !(ts.getCombinedModifierFlags(m) & PRIVATE);
    this.add(m, `${owner}.${m.name.getText(this.sf)}`, kind, exported, owner, kind === 'getter' ? 'get' : kind === 'setter' ? 'set' : undefined);
  }

  /** Overloads and duplicate names share one declaration with the signatures joined. */
  private add(node: ts.Node, name: string, kind: string, exported: boolean, container: string | undefined, accessor?: 'get' | 'set'): TsDecl {
    const key = accessor ? `${name}:${accessor}` : name;
    const docs = ts.getJSDocCommentsAndTags(node);
    const deprecated = docs.some((d) => /@deprecated\b/i.test(d.getFullText(this.sf)));
    const signature = this.signature(node);
    const existing = this.byName.get(key);
    if (existing) {
      if (signature && !existing.signature.includes(signature)) existing.signature = `${existing.signature} | ${signature}`.slice(0, SIGNATURE_CAP);
      existing.documented ||= docs.length > 0;
      existing.deprecated ||= deprecated;
      existing.exported ||= exported;
      return existing;
    }
    const decl: TsDecl = {name, kind, accessor, container, exported, line: this.sf.getLineAndCharacterOfPosition(this.start(node)).line + 1,
      signature, documented: docs.length > 0, deprecated, extends: [], implements: []};
    this.byName.set(key, decl);
    this.file.decls.push(decl);
    return decl;
  }

  /** The head of the declaration up to its body, decorators left out, whitespace collapsed, at most 160 characters. */
  private signature(node: ts.Node): string {
    const start = this.start(node);
    const body = (node as any).body as ts.Node | undefined;
    const raw = body && ts.isFunctionLike(node) ? this.text.slice(start, body.getStart(this.sf)) : this.text.slice(start, Math.min(node.getEnd(), start + SIGNATURE_CAP * 2)).split('{')[0];
    return raw.replace(/\s+/g, ' ').replace(/\s*[=;,]?\s*$/, '').trim().slice(0, SIGNATURE_CAP);
  }

  private start(node: ts.Node): number {
    const decorators = ts.canHaveDecorators(node) ? ts.getDecorators(node) : undefined;
    return decorators?.length ? skipSpace(this.text, decorators[decorators.length - 1].getEnd()) : node.getStart(this.sf);
  }

  private isExported(st: ts.Statement): boolean {
    if (ts.getCombinedModifierFlags(st as unknown as ts.Declaration) & ts.ModifierFlags.Export) return true;
    const name = (st as any).name;
    if (name && ts.isIdentifier(name)) return this.exportList.has(name.text);
    return ts.isVariableStatement(st) && st.declarationList.declarations.some((d) => ts.isIdentifier(d.name) && this.exportList.has(d.name.text));
  }

  private importOf(st: ts.ImportDeclaration | ts.ExportDeclaration | ts.ImportEqualsDeclaration): void {
    const imp: TsImport = {specifier: '', symbols: [], locals: {}};
    if (ts.isImportEqualsDeclaration(st)) {
      if (!ts.isExternalModuleReference(st.moduleReference) || !ts.isStringLiteral(st.moduleReference.expression)) return;
      imp.specifier = st.moduleReference.expression.text;
      imp.symbols.push('*');
      imp.namespace = st.name.text;
    }
    else {
      if (!st.moduleSpecifier || !ts.isStringLiteral(st.moduleSpecifier)) return;
      imp.specifier = st.moduleSpecifier.text;
      if (ts.isExportDeclaration(st)) imp.reexport = true;
      const clause = ts.isImportDeclaration(st) ? st.importClause : st.exportClause;
      if (ts.isImportDeclaration(st) && clause && ts.isImportClause(clause) && clause.name) {
        imp.symbols.push('default');
        imp.locals[clause.name.text] = 'default';
      }
      const bindings = clause && ts.isImportClause(clause) ? clause.namedBindings : clause;
      if (!bindings) {
        if (ts.isExportDeclaration(st)) imp.symbols.push('*');
      }
      else if (ts.isNamespaceImport(bindings) || ts.isNamespaceExport(bindings)) {
        imp.symbols.push('*');
        if (ts.isNamespaceImport(bindings)) imp.namespace = bindings.name.text;
        else imp.alias = bindings.name.text;
      }
      else
        for (const e of bindings.elements) {
          const imported = (e.propertyName ?? e.name).text;
          imp.symbols.push(imported);
          imp.locals[e.name.text] = imported;
        }
    }
    this.file.imports.push(imp);
  }

  private visitDynamicImports(node: ts.Node): void {
    if (ts.isCallExpression(node) && node.expression.kind === ts.SyntaxKind.ImportKeyword && node.arguments.length && ts.isStringLiteral(node.arguments[0]))
      this.file.imports.push({specifier: node.arguments[0].text, symbols: [], locals: {}});
    ts.forEachChild(node, (child) => this.visitDynamicImports(child));
  }
}

export const declarationsExtractor: Extractor = {
  name: 'ts-declarations',
  describes: {'ts-declarations': 'source files, declarations and their inheritance'},
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const sources = tsSources(ctx, emitter);
    for (const file of sources.files) {
      const isPackage = file.unit.id.startsWith('pkg:');
      const isApi = file.unit.id === JS_API;
      const admitted = emitter.node(sourceFileRow(file.path, {loc: file.loc, generated: file.generated ? true : undefined, package: isPackage ? file.unit.id : undefined,
        entry: file.unit.entries?.has(file.path) ? true : undefined}));
      if (!admitted.accepted) continue;
      emitter.edge({type: 'declares', from: file.unit.id, to: fileId(file.path), derived_by: 'ast', confidence: 1, evidence: [file.path]});
      const emitted = new Set<string>();
      for (const decl of file.decls) {
        const publicApi = isApi && decl.exported;
        if (file.generated && decl.container && !publicApi) continue;
        if (!isNode(file, decl)) continue;
        const id = sources.declId(file, decl);
        // a type inside a namespace the graph does not hold is declared by the file, not by the namespace
        const container = decl.container !== undefined && emitted.has(decl.container) ? declId(file.path, decl.container) : undefined;
        if (!emitter.node({type: 'declaration', id, name: decl.name.slice(decl.name.lastIndexOf('.') + 1), kind: decl.kind, exported: decl.exported,
          public_api: publicApi ? true : undefined, generated: file.generated ? true : undefined, deprecated: decl.deprecated ? true : undefined, documented: decl.documented,
          signature: decl.signature || undefined, line: decl.line, language: 'ts', path: file.path, provenance: 'ast', source_layer: 'public'}).accepted) continue;
        emitted.add(decl.name);
        emitter.edge({type: 'declares', from: container ?? fileId(file.path), to: id, derived_by: 'ast', confidence: 1, evidence: [file.path]});
        for (const [type, names] of [['extends', decl.extends], ['implements', decl.implements]] as const)
          for (const name of names) {
            const to = sources.resolveHeritage(file, name, decl.container);
            if (to) emitter.edge({type, from: id, to, derived_by: 'ast', confidence: 1, evidence: [file.path]});
            else emitter.problem('unresolved_ids', `${file.path}:${decl.line}: ${decl.name} ${type} ${name}, which no file in reach declares`);
          }
      }
    }
    emitter.source('ts-declarations', sources.failed ? 'partial' : 'ok');
  },
};

function leftmostIdentifier(expression: ts.Expression): string | undefined {
  let node: ts.Node = expression;
  for (;;) {
    if (ts.isNewExpression(node) || ts.isCallExpression(node) || ts.isPropertyAccessExpression(node) || ts.isAsExpression(node) || ts.isParenthesizedExpression(node)) node = node.expression;
    else return ts.isIdentifier(node) ? node.text : undefined;
  }
}

function skipSpace(text: string, at: number): number {
  while (at < text.length && /\s/.test(text[at])) at++;
  return at;
}
