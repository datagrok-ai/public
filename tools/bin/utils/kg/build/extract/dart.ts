/// The Dart sources of `core/` (build-plan.md WO-6): one lexical pass per file — the file itself, its
/// top-level declarations, the tests of a test file and the `~id` markers of conventions.md §6. No
/// analyzer and no AST: no members, no imports, no calls, no heritage, and nothing a regex cannot see.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../emitter';
import {BuildContext, Extractor} from '../context';
import {HomeSet} from '../../homes';
import {fileId, declId, testId, suiteId, docId, docKind, HELP_DIR, helpPage, countLines, sourceFileRow} from '../../ids';
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
const TEST_FILE = /(?:^|\/)test\/|_test\.dart$/;
const TEST_CALL = /(?<![\w.$])(group|test)\s*\(\s*(['"])((?:\\.|(?!\2).)*)\2/g;

/** The table of help urls every other file reaches a page through; one of its constants lacks the leading slash. */
const HELP_TABLE = 'core/shared/grok_shared/lib/src/help_url.dart';
const HELP_CONST = /static\s+const\s+String\s+([A-Za-z_$][\w$]*)\s*=\s*'(\/?help\/[^']+)'/g;
const HELP_REF = /\bHelpUrl\.([A-Za-z_$][\w$]*)/g;
const HELP_LITERAL = /'(\/help\/[^'\s]+)'/g;
/** A help page a doc comment names by its repo path, as the viewer cores do. */
const HELP_DOC_PATH = /^\s*\/\/\/.*?(public\/help\/[\w./-]+\.mdx?)/;

export const dartExtractor: Extractor = {
  name: 'dart',
  describes: {dart: 'Dart source files, their top-level types, tests and markers (a lexical pass)'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const homes = homesOf(ctx);
    const helpUrls = helpConstants(ctx.repoRoot);
    const files = globSync(SOURCES, {cwd: ctx.repoRoot, ignore: SOURCE_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
    const perPackage = new Map<string, number>();
    let unresolved = 0;
    for (const file of files) {
      const full = path.join(ctx.repoRoot, file);
      const text = fs.readFileSync(full, 'utf8');
      const lines = text.split(/\r?\n/);
      const generated = file.endsWith('.g.dart');
      if (!emitter.node(sourceFileRow(file, {loc: countLines(text), generated: generated ? true : undefined})).accepted) continue;
      const pkg = packageOf(file);
      if (pkg) perPackage.set(pkg, (perPackage.get(pkg) ?? 0) + 1);
      declarations(emitter, file, lines, generated);
      if (TEST_FILE.test(file)) tests(emitter, file, text);
      unresolved += markers(emitter, homes, file, lines);
      helpRefs(emitter, ctx.repoRoot, file, text, lines, helpUrls);
    }
    emitter.manifest('dart_packages', coverage(ctx.repoRoot, perPackage));
    emitter.manifest('dart_depth', 'lexical');
    emitter.source('dart', unresolved ? 'partial' : 'ok');
  },
};

/** The top-level types of a file, with the doc comment and the annotations above each of them. */
function declarations(emitter: Emitter, file: string, lines: string[], generated: boolean): void {
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
  }
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

/** One suite per test file, and its tests under the `group` titles that enclose them. */
function tests(emitter: Emitter, file: string, text: string): void {
  const found = parseDartTests(text);
  if (!found.length) return;
  const suite = suiteId('dart', file);
  emitter.node({type: 'test-suite', id: suite, name: path.posix.basename(file), framework: 'dart', path: file,
    provenance: 'filesystem', source_layer: 'core'});
  for (const t of found)
    emitter.node({type: 'test', id: testId('dart', file, t.category, t.name), name: t.name, path: file, framework: 'dart',
      level: 'unit', category: t.category || undefined, suite, provenance: 'ast', source_layer: 'core'});
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
