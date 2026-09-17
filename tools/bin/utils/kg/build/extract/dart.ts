/// The Dart sources of `core/` (build-plan.md WO-6): one lexical pass per file — the file itself, its
/// top-level declarations, its import directives, the tests of a test file (`test()` or the client's
/// `regTest()`) and the `~id` markers of conventions.md §6. No analyzer and no AST: no members, no calls,
/// no heritage, and nothing a regex cannot see.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../emitter';
import {BuildContext, Extractor} from '../context';
import {HomeSet} from '../../homes';
import {fileId, declId, testId, suiteId, docId, docKind, HELP_DIR, DART_PACKAGES, REG_TEST_DIR, helpPage, countLines, sourceFileRow} from '../../ids';
import {homesOf, resolveMention, MARKER_LINE} from './markers';
import {blankComments, matchBrace} from './ts/tests';

const SOURCES = 'core/**/*.dart';
const SOURCE_IGNORE = ['**/.dart_tool/**', '**/build/**', '**/packages/**', '**/node_modules/**'];
/** A Dart package of the checkout: a folder under `core/` with a pubspec.yaml, `core/<area>/libs/<pkg>` included. */
const PUBSPECS = ['core/*/*/pubspec.yaml', 'core/*/*/*/pubspec.yaml'];
/** A top-level declaration: at column 0, so a nested class or a string holding the word is not one. */
const TYPE_DECL = /^(?:abstract\s+)?(class|mixin|enum)\s+([A-Za-z_$][\w$]*)/;
/** Dart 1 spells an alias `typedef void Action(x)` and Dart 2 `typedef Action = ...`: the name precedes the parameters. */
const TYPEDEF = /^typedef\s+.*?([A-Za-z_$][\w$]*)\s*(?:<[^(]*>)?\s*[(=]/;
const KINDS: Record<string, string> = {class: 'class', mixin: 'mixin', enum: 'enum'};
/** `/// ~id` opening a doc comment: the ownership marker of §6 in its Dart spelling. */
const OWN_MARKER = /^\/\/\/\s*~((?:[A-Z][A-Za-z]{0,5}:)?[a-z][a-z0-9]*(?:-[a-z0-9]+)*(?:\/[a-z0-9]+(?:-[a-z0-9]+)*)*)(?:#[\w-]+)?\s*$/;
const ANNOTATION = /^@/;
const DEPRECATED = /^@(?:deprecated\b|Deprecated\()/;
const TEST_DIR = /(?:^|\/)test\//;
const TEST_FILE = /(?:^|\/)test\/|_test\.dart$/;
const TEST_CALL = /(?<![\w.$])(group|test)\s*\(\s*(['"])((?:\\.|(?!\2).)*)\2/g;
/** The client's browser tests: `regTest('Area | Group | Name', ...)`, run by DevTools rather than `dart test`. */
const REG_TEST_CALL = /(?<![\w.$])regTest\s*\(\s*(['"])((?:\\.|(?!\1).)*)\1/g;
/** An `import`, `export` or `part` directive: its first specifier, then whatever clauses follow it. */
const DIRECTIVE = /^\s*(import|export|part)\s+(['"])([^'"]+)\2([^;]*)/;
const SHOW = /\bshow\s+([\w$,\s]+)/;
/** A capitalized word of a test file, the lexical `uses` candidate of §8.1; shorter than MIN_TOKEN it is noise. */
const TYPE_TOKEN = /\b[A-Z][A-Za-z0-9_]*\b/g;
const MIN_TOKEN = 4;
/** String literals, blanked before tokenizing: a type name inside a function expression or a message is not a use. */
const STRING_LITERAL = /r?(?:'''[\s\S]*?'''|"""[\s\S]*?"""|'(?:[^'\\\n]|\\.)*'|"(?:[^"\\\n]|\\.)*")/g;
const LEXICAL_CONFIDENCE = 0.7;

/** The table of help urls every other file reaches a page through; one of its constants lacks the leading slash. */
const HELP_TABLE = 'core/shared/grok_shared/lib/src/help_url.dart';
const HELP_CONST = /static\s+const\s+String\s+([A-Za-z_$][\w$]*)\s*=\s*'(\/?help\/[^']+)'/g;
const HELP_REF = /\bHelpUrl\.([A-Za-z_$][\w$]*)/g;
const HELP_LITERAL = /'(\/help\/[^'\s]+)'/g;
/** A help page a doc comment names by its repo path, as the viewer cores do. */
const HELP_DOC_PATH = /^\s*\/\/\/.*?(public\/help\/[\w./-]+\.mdx?)/;

export const dartExtractor: Extractor = {
  name: 'dart',
  describes: {dart: 'Dart source files, their top-level types, imports, tests and markers (a lexical pass)'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const homes = homesOf(ctx);
    const helpUrls = helpConstants(ctx.repoRoot);
    const files = globSync(SOURCES, {cwd: ctx.repoRoot, ignore: SOURCE_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
    const known = new Set(files);
    const perPackage = new Map<string, number>();
    const declared = new Map<string, Map<string, string[]>>();
    const testFiles: {file: string, text: string, packages: string[]}[] = [];
    const libraryOf = new Map<string, string>();
    const packagesOf = new Map<string, string[]>();
    let unresolved = 0;
    for (const file of files) {
      const full = path.join(ctx.repoRoot, file);
      const text = fs.readFileSync(full, 'utf8');
      const lines = text.split(/\r?\n/);
      const generated = file.endsWith('.g.dart');
      const directives = lines.map((l) => DIRECTIVE.exec(l)).filter((m): m is RegExpExecArray => !!m);
      const entry = directives.some((m) => m[1] === 'part');
      if (!emitter.node(sourceFileRow(file, {loc: countLines(text), generated: generated ? true : undefined, entry: entry ? true : undefined})).accepted) continue;
      const pkg = packageOf(file);
      if (pkg) perPackage.set(pkg, (perPackage.get(pkg) ?? 0) + 1);
      const names = declarations(emitter, file, lines, generated);
      if (pkg && names.length) {
        const index = declared.get(pkg) ?? new Map<string, string[]>();
        for (const name of names) index.set(name, [...index.get(name) ?? [], file]);
        declared.set(pkg, index);
      }
      const packages = imports(emitter, file, directives, known);
      packagesOf.set(file, packages);
      for (const m of directives)
        if (m[1] === 'part' && importTarget(file, m[3]) !== undefined) libraryOf.set(importTarget(file, m[3])!, file);
      const testing = file.startsWith(REG_TEST_DIR) ? regTests(emitter, file, text) : TEST_FILE.test(file) && tests(emitter, file, text);
      // a helper of a test folder (a server setup, a fixture) uses types the tests importing it never name
      if (testing || TEST_DIR.test(file)) testFiles.push({file, text, packages: [...new Set([pkg, ...packages])].filter((p): p is string => p !== undefined)});
      unresolved += markers(emitter, homes, file, lines);
      helpRefs(emitter, ctx.repoRoot, file, text, lines, helpUrls);
    }
    for (const t of testFiles)
      lexicalUses(emitter, t.file, t.text, [...new Set([...t.packages, ...packagesOf.get(libraryOf.get(t.file) ?? '') ?? []])], declared);
    emitter.manifest('dart_packages', coverage(ctx.repoRoot, perPackage));
    emitter.manifest('dart_depth', 'lexical');
    emitter.source('dart', unresolved ? 'partial' : 'ok');
  },
};

/** The top-level types of a file, with the doc comment and the annotations above each of them; returns their names. */
function declarations(emitter: Emitter, file: string, lines: string[], generated: boolean): string[] {
  const names: string[] = [];
  for (let i = 0; i < lines.length; i++) {
    const declared = declarationOf(lines[i]);
    if (!declared) continue;
    const {name, kind} = declared;
    const id = declId(file, name);
    const {documented, deprecated} = above(lines, i);
    if (!emitter.node({type: 'declaration', id, name, kind, exported: !name.startsWith('_'),
      generated: generated ? true : undefined, documented, deprecated: deprecated ? true : undefined,
      line: i + 1, language: 'dart', path: file, provenance: 'ast', source_layer: 'core'}).accepted) continue;
    emitter.edge({type: 'declares', from: fileId(file), to: id, derived_by: 'ast', confidence: 1, evidence: [file]});
    names.push(name);
  }
  return names;
}

function declarationOf(line: string): {name: string, kind: string} | undefined {
  const type = TYPE_DECL.exec(line);
  if (type) return {name: type[2], kind: KINDS[type[1]]};
  const alias = TYPEDEF.exec(line);
  return alias ? {name: alias[1], kind: 'type'} : undefined;
}

/** What sits above a declaration: whether its annotations deprecate it, and whether the line before them documents it. */
function above(lines: string[], at: number): {documented: boolean, deprecated: boolean} {
  let deprecated = false;
  for (let i = at - 1; i >= 0; i--) {
    const line = lines[i].trim();
    if (!line) continue;
    if (!ANNOTATION.test(line)) return {documented: line.startsWith('///'), deprecated};
    deprecated = deprecated || DEPRECATED.test(line);
  }
  return {documented: false, deprecated};
}

/** One `imports` edge per file a directive names and the pass walked, with the names of a `show` clause (`part` for a
 * part); a `dart:` library, a pub package or a file outside the walk resolves to nothing and is not counted. Returns
 * the names of the `package:` imports the table knows. */
function imports(emitter: Emitter, file: string, directives: RegExpExecArray[], known: Set<string>): string[] {
  const targets = new Map<string, Set<string>>();
  const packages = new Set<string>();
  for (const m of directives) {
    const pkg = /^package:([^/]+)\//.exec(m[3])?.[1];
    if (pkg !== undefined && DART_PACKAGES[pkg] !== undefined) packages.add(pkg);
    const to = importTarget(file, m[3]);
    if (to === undefined || to === file || !known.has(to)) continue;
    let symbols = targets.get(to);
    if (!symbols) targets.set(to, symbols = new Set());
    if (m[1] === 'part') symbols.add('part');
    for (const s of SHOW.exec(m[4])?.[1].split(',') ?? [])
      if (s.trim()) symbols.add(s.trim());
  }
  for (const [to, symbols] of targets)
    emitter.edge({type: 'imports', from: fileId(file), to: fileId(to), symbols: symbols.size ? [...symbols].sort() : undefined, derived_by: 'ast', confidence: 1, evidence: [file]});
  return [...packages];
}

/** §8.1 uses, lexically: a capitalized word of a test file that exactly one type of its own package or of a package it
 * imports declares; a name several files declare is counted as ambiguous and drawn for nothing. */
function lexicalUses(emitter: Emitter, file: string, text: string, packages: string[], declared: Map<string, Map<string, string[]>>): void {
  const candidates = new Map<string, string[]>();
  for (const pkg of packages)
    for (const [name, files] of declared.get(pkg) ?? [])
      candidates.set(name, [...candidates.get(name) ?? [], ...files.filter((f) => f !== file)]);
  const seen = new Set<string>();
  for (const m of blankComments(text).replace(STRING_LITERAL, (s) => ' '.repeat(s.length)).matchAll(TYPE_TOKEN)) {
    const name = m[0];
    if (name.length < MIN_TOKEN || seen.has(name)) continue;
    seen.add(name);
    const files = candidates.get(name);
    if (!files?.length) continue;
    if (files.length > 1) emitter.problem('ambiguous_uses', `${file}: ${name} is declared in ${files.join(' and ')}`);
    else emitter.edge({type: 'uses', from: fileId(file), to: declId(files[0], name), kind: 'type', derived_by: 'lexical', confidence: LEXICAL_CONFIDENCE, evidence: [file]});
  }
}

/** `package:<pkg>/<path>` through the package table into `<dir>/lib/<path>`; a relative specifier against the file. */
function importTarget(file: string, specifier: string): string | undefined {
  const pkg = /^package:([^/]+)\/(.+)$/.exec(specifier);
  if (pkg) return DART_PACKAGES[pkg[1]] === undefined ? undefined : `${DART_PACKAGES[pkg[1]]}/lib/${pkg[2]}`;
  return specifier.includes(':') ? undefined : path.posix.normalize(path.posix.join(path.posix.dirname(file), specifier));
}

/** One suite per test file, and its tests under the `group` titles that enclose them; whether the file holds any. */
function tests(emitter: Emitter, file: string, text: string): boolean {
  const found = parseDartTests(text);
  if (!found.length) return false;
  const suite = suiteId('dart', file);
  emitter.node({type: 'test-suite', id: suite, name: path.posix.basename(file), framework: 'dart', path: file,
    provenance: 'filesystem', source_layer: 'core'});
  for (const t of found) {
    const id = testId('dart', file, t.category, t.name);
    if (emitter.node({type: 'test', id, name: t.name, path: file, framework: 'dart',
      level: 'unit', category: t.category || undefined, suite, provenance: 'ast', source_layer: 'core'}).accepted)
      emitter.edge({type: 'declares', from: fileId(file), to: id, derived_by: 'ast', confidence: 1, evidence: [file]});
  }
  return true;
}

/** `regTest('A | B | C'` under the client tests folder: category `A | B`, name `C`, one suite per file. Whether a
 * category is a d4 or an xamgle one is the runner's call, made from the text later. `${DartLibraryTestCategoryName.x}`
 * is the constant `x` (tests.dart); any other interpolation makes the test dynamic and keeps the title as written. */
function regTests(emitter: Emitter, file: string, text: string): boolean {
  const found = [...blankComments(text).matchAll(REG_TEST_CALL)]
    .map((m) => m[2].trim().replace(/\$\{DartLibraryTestCategoryName\.(\w+)\}/g, '$1'));
  if (!found.length) return false;
  const suite = suiteId('xamgle', file);
  emitter.node({type: 'test-suite', id: suite, name: path.posix.basename(file), framework: 'xamgle', path: file,
    provenance: 'filesystem', source_layer: 'core'});
  for (const title of found) {
    const segments = title.split('|').map((s) => s.trim());
    const name = segments.pop()!;
    const category = segments.join(' | ') || undefined;
    const id = testId('xamgle', file, category ?? '', name);
    if (emitter.node({type: 'test', id, name, path: file, framework: 'xamgle', level: 'e2e', category, suite,
      dynamic: title.includes('$') ? true : undefined, provenance: 'ast', source_layer: 'core'}).accepted)
      emitter.edge({type: 'declares', from: fileId(file), to: id, derived_by: 'ast', confidence: 1, evidence: [file]});
  }
  return true;
}

/** `test('name'` calls with the `group('title'` chain around each; a commented-out test is not a test. */
export function parseDartTests(source: string): {category: string, name: string}[] {
  const text = blankComments(source);
  const groups: {title: string, start: number, end: number}[] = [];
  const out: {category: string, name: string}[] = [];
  for (const m of text.matchAll(TEST_CALL)) {
    const at = m.index!;
    if (m[1] === 'group') {
      const open = text.indexOf('{', at + m[0].length);
      if (open >= 0) groups.push({title: m[3].trim(), start: open, end: matchBrace(text, open)});
      continue;
    }
    out.push({category: groups.filter((g) => g.start < at && at < g.end).map((g) => g.title).join(' > '), name: m[3].trim()});
  }
  return out;
}

/** §6 in Dart: `/// ~id` opening a doc comment owns the file, `// ~id` alone on a line only participates. Returns
 * how many markers named a feature no home declares, which are counted and drawn for nothing. */
function markers(emitter: Emitter, homes: HomeSet, file: string, lines: string[]): number {
  let unresolved = 0;
  const seen = new Set<string>();
  for (let i = 0; i < lines.length; i++) {
    const owns = OWN_MARKER.exec(lines[i].trim());
    // a `///` line inside a doc comment is prose about the feature, not the marker that opens one
    if (owns && i > 0 && lines[i - 1].trim().startsWith('///')) continue;
    const token = owns ? owns[1] : MARKER_LINE.exec(lines[i])?.[1];
    if (token === undefined || seen.has(`${owns ? 'owns' : 'in'}:${token}`)) continue;
    seen.add(`${owns ? 'owns' : 'in'}:${token}`);
    const target = resolveMention(emitter, homes, token, `${file}:${i + 1}`);
    if (!target || target.root !== 'feature') {
      unresolved++;
      continue;
    }
    emitter.claim({file, feature: target.id, rung: 1, source: 'marker', mode: owns ? undefined : 'participates', props: {}, line: i + 1});
  }
  return unresolved;
}

/** `HelpUrl.Name` -> the url it stands for, read once from the table in grok_shared. */
function helpConstants(repoRoot: string): Map<string, string> {
  const out = new Map<string, string>();
  const full = path.join(repoRoot, HELP_TABLE);
  if (!fs.existsSync(full)) return out;
  for (const m of fs.readFileSync(full, 'utf8').matchAll(HELP_CONST)) out.set(m[1], m[2]);
  return out;
}

/** The help pages a file names, through the table, a `/help/...` literal or a `public/help/...` path in a doc comment.
 * A page that resolves is recorded for membership to draw `documents` from; one that does not is a stale reference. */
function helpRefs(emitter: Emitter, repoRoot: string, file: string, text: string, lines: string[], constants: Map<string, string>): void {
  const urls = new Set<string>();
  for (const m of text.matchAll(HELP_REF)) {
    const url = constants.get(m[1]);
    if (url) urls.add(url);
  }
  for (const m of text.matchAll(HELP_LITERAL)) urls.add(m[1]);
  for (const line of lines) {
    const m = HELP_DOC_PATH.exec(line);
    if (m) urls.add(m[1].slice('public'.length));
  }
  const pages = new Set<string>();
  for (const raw of urls) {
    const url = raw.startsWith('/') ? raw : `/${raw}`;
    const page = helpPage(repoRoot, url);
    if (!page) emitter.problem('unresolved_ids', `${file}: help-url ${url} names no page under ${HELP_DIR}`);
    else if (!pages.has(page)) {
      pages.add(page);
      emitter.stub(docId(page), 'doc-page', path.posix.basename(page), 'ast', {path: page, kind: docKind(page)});
      emitter.helpRef(file, page);
    }
  }
}

/** What the pass covered, for the manifest: every Dart package of the checkout with the number of files it holds. */
function coverage(repoRoot: string, perPackage: Map<string, number>): Record<string, number> {
  const known = PUBSPECS.flatMap((p) => globSync(p, {cwd: repoRoot, posix: true, windowsPathsNoEscape: true})).map(packageOf);
  const out: Record<string, number> = {};
  for (const name of [...new Set([...perPackage.keys(), ...known])].filter((n): n is string => n !== undefined).sort())
    out[name] = perPackage.get(name) ?? 0;
  return out;
}

/** `core/shared/ddt/lib/x.dart` and `core/shared/ddt/pubspec.yaml` are both `ddt`; `core/server/libs/shelf/...` is `shelf`. */
function packageOf(file: string): string | undefined {
  const segments = file.split('/');
  if (segments[0] !== 'core' || segments.length < 3) return undefined;
  return segments[2] === 'libs' ? segments[3] : segments[2];
}
