/// Automated tests under `public/` (build-plan.md WO-3c): DG `category()`/`test()` calls in package sources
/// and Playwright `describe`/`test` titles in playwright-public, package playwright folders and Test Track
/// specs, each in its suite, with a `tests` edge when a category, tag or title names a feature.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../../registry';
import {pkgId, testId, suiteId} from '../../ids';
import {HomeIndex, homesOf, leadingId, resolveMention} from '../markers';
import {listPackages} from './packages';

export interface DgTest {
  category: string;
  name: string;
  skipReason?: string;
  benchmark: boolean;
  tags: string[];
}

export interface PlaywrightTest {
  /** Enclosing `describe` titles, outermost first. */
  describes: string[];
  title: string;
  skipped: boolean;
}

const SOURCE_IGNORE = ['**/node_modules/**', '**/dist/**'];
const PLAYWRIGHT_GLOBS = ['public/playwright-public/**/*.test.ts', 'public/packages/*/playwright/**/*.test.ts',
  'public/packages/UsageAnalysis/files/TestTrack/**/*{-spec,.test}.ts'];
const API_TESTS_PACKAGE = 'ApiTests';
const DG_CALL = /(?<![\w.$])(category|test)\s*\(\s*(['"`])((?:\\.|(?!\2).)*)\2/g;
const PLAYWRIGHT_CALL = /(?<![\w.$])(test\.describe(?:\.(?:serial|parallel|only|skip|fixme))*|describe(?:\.(?:only|skip))?|test(?:\.(?:only|skip|fixme|fail))?)\s*\(\s*(['"`])((?:\\.|(?!\2).)*)\2/g;
/** A trailing `{...}` argument of a call, the DG test options when its keys are TestOptions keys (utils/src/test.ts). */
const TRAILING_OBJECT = /[\w)\]}'"`]\s*,\s*\{([^{}]*)\}\s*\)/g;
const OPTION_KEYS = /\b(?:timeout|skipReason|benchmark|tags|stressTest|owner|isAggregated|benchmarkTimeout|benchmarkWarnTimeout|unhandledExceptionTimeout)\s*:/;
const STRING = /^(['"`])((?:\\.|(?!\1).)*)\1$/;

export const testsExtractor: Extractor = {
  name: 'ts-tests',
  layer: 'public',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const layer = new TestLayer(ctx.repoRoot, emitter, new HomeIndex(homesOf(ctx)));
    for (const pkg of listPackages(ctx.repoRoot))
      for (const file of layer.glob(`${pkg.dir}/src/**/*.ts`))
        if (!file.endsWith('.d.ts')) layer.emitDgFile(pkg.folder, file);
    for (const pattern of PLAYWRIGHT_GLOBS)
      for (const file of layer.glob(pattern)) layer.emitPlaywrightFile(file);
  },
};

class TestLayer {
  private suites = new Set<string>();

  constructor(private repoRoot: string, private emitter: Emitter, private index: HomeIndex) {}

  emitDgFile(pkg: string, file: string): void {
    const level = pkg === API_TESTS_PACKAGE ? 'api' : 'unit';
    for (const t of parseDgTests(this.read(file))) {
      const id = testId('dg', file, t.category, t.name);
      const row: Row = {type: 'test', id, name: t.name, path: file, framework: 'dg', level, category: t.category, skipped: t.skipReason !== undefined ? true : undefined,
        skip_reason: t.skipReason, benchmark: t.benchmark ? true : undefined, tags: t.tags.length ? t.tags : undefined, provenance: 'ast', source_layer: 'public'};
      this.emitter.node(row);
      const suite = suiteId('dg', pkg, t.category);
      if (!this.suites.has(suite)) {
        this.suites.add(suite);
        this.emitter.node({type: 'test-suite', id: suite, name: t.category, framework: 'dg', package: pkgId(pkg), provenance: 'ast', source_layer: 'public'});
      }
      this.emitter.edge({type: 'in-suite', from: id, to: suite, derived_by: 'ast', confidence: 1, evidence: [file]});
      const feature = leadingId(t.category) ?? (t.tags[0] === undefined ? undefined : leadingId(t.tags[0]));
      if (feature) this.tests(id, feature, level, file);
    }
  }

  emitPlaywrightFile(file: string): void {
    const tests = parsePlaywrightTests(this.read(file));
    if (!tests.length) return;
    const suite = suiteId('playwright', file);
    const pkg = /^public\/packages\/([^/]+)\//.exec(file)?.[1];
    this.emitter.node({type: 'test-suite', id: suite, name: path.posix.basename(file), framework: 'playwright', path: file, package: pkg ? pkgId(pkg) : undefined,
      provenance: 'filesystem', source_layer: 'public'});
    for (const t of tests) {
      const chain = t.describes.join(' > ');
      const id = playwrightTestId(file, t);
      this.emitter.node({type: 'test', id, name: chain ? `${chain} > ${t.title}` : t.title, path: file, framework: 'playwright', level: 'e2e', category: chain || undefined,
        skipped: t.skipped ? true : undefined, provenance: 'ast', source_layer: 'public'});
      this.emitter.edge({type: 'in-suite', from: id, to: suite, derived_by: 'ast', confidence: 1, evidence: [file]});
      const feature = leadingId(t.title);
      if (feature) this.tests(id, feature, 'e2e', file);
    }
  }

  /** The `~id` of a category, a tag or a title is a marker like any other: one no home declares is counted, not drawn. */
  private tests(test: string, token: string, level: string, file: string): void {
    const feature = resolveMention(this.emitter, this.index, token, file);
    if (feature) this.emitter.edge({type: 'tests', from: test, to: feature.id, kind: level, derived_by: 'annotation', confidence: 1, evidence: [file]});
  }

  glob(pattern: string): string[] {
    return globSync(pattern, {cwd: this.repoRoot, ignore: SOURCE_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
  }

  private read(file: string): string {
    return fs.readFileSync(path.join(this.repoRoot, file), 'utf8');
  }
}

/** `test:playwright:<file>#<describe chain>/<title>`; a test outside any describe sits under the file stem. */
export function playwrightTestId(file: string, t: PlaywrightTest): string {
  return testId('playwright', file, t.describes.join(' > ') || path.posix.basename(file).replace(/\.[^.]+$/, ''), t.title);
}

/** The DG tests of a source, in order: `category('X'` sets the category of every `test('Y'` after it, as the framework does; a
 * test before any category is not registered and is left out. A title built from a concatenation (`test('Correctness: ' + name`)
 * keeps only its literal head, trimmed, so the id of a row and the name it carries stay the same text. */
export function parseDgTests(source: string): DgTest[] {
  const text = blankComments(source);
  const calls = [...text.matchAll(DG_CALL)];
  const out: DgTest[] = [];
  let category: string | undefined;
  calls.forEach((m, i) => {
    if (m[1] === 'category') {
      category = m[3].trim();
      return;
    }
    if (category === undefined) return;
    const span = text.slice(m.index!, calls[i + 1]?.index ?? text.length);
    const options = [...span.matchAll(TRAILING_OBJECT)].map((o) => o[1]).filter((o) => OPTION_KEYS.test(o)).pop() ?? '';
    const skip = /\bskipReason\s*:\s*((['"`])(?:\\.|(?!\2).)*\2|[^,}]+)/.exec(options);
    const tags = /\btags\s*:\s*\[([^\]]*)\]/.exec(options);
    out.push({
      category, name: m[3].trim(),
      skipReason: skip ? (STRING.exec(skip[1].trim())?.[2] ?? skip[1].trim()) : undefined,
      benchmark: /\bbenchmark\s*:\s*true\b/.test(options),
      tags: tags ? tags[1].split(',').map((t) => STRING.exec(t.trim())?.[2]).filter((t): t is string => !!t) : [],
    });
  });
  return out;
}

/** The Playwright tests of a spec with their enclosing describes; `test.skip('x'`, `test.fixme('x'` and a `describe.skip` make a test skipped. */
export function parsePlaywrightTests(source: string): PlaywrightTest[] {
  const text = blankComments(source);
  const describes: {title: string, start: number, end: number, skipped: boolean}[] = [];
  const out: PlaywrightTest[] = [];
  for (const m of text.matchAll(PLAYWRIGHT_CALL)) {
    const call = m[1];
    const at = m.index!;
    if (call.includes('describe')) {
      const open = callbackBrace(text, at + m[0].length);
      if (open < 0) continue;
      describes.push({title: m[3].trim(), start: open, end: matchBrace(text, open), skipped: /\.skip\b/.test(call)});
      continue;
    }
    const enclosing = describes.filter((d) => d.start < at && at < d.end);
    out.push({describes: enclosing.map((d) => d.title), title: m[3].trim(), skipped: /\.(skip|fixme)\b/.test(call) || enclosing.some((d) => d.skipped)});
  }
  return out;
}

/** The `{` opening the callback of a call whose title ends at [from]: the first brace after `=>` or `function(...)`. */
function callbackBrace(text: string, from: number): number {
  const arrow = /=>|function\b[^(]*\([^)]*\)/g;
  arrow.lastIndex = from;
  const m = arrow.exec(text);
  return m ? text.indexOf('{', m.index + m[0].length) : -1;
}

/** Index of the `}` closing the brace at [open], strings skipped; the end of the text when unbalanced. */
function matchBrace(text: string, open: number): number {
  let depth = 0;
  let quote = '';
  for (let i = open; i < text.length; i++) {
    const c = text[i];
    if (quote) {
      if (c === '\\') i++;
      else if (c === quote) quote = '';
      continue;
    }
    if (c === '\'' || c === '"' || c === '`') quote = c;
    else if (c === '{') depth++;
    else if (c === '}' && --depth === 0) return i;
  }
  return text.length;
}

/** [source] with every line and block comment blanked out (line breaks kept), so a commented-out test is not a test. */
export function blankComments(source: string): string {
  let out = '';
  let quote = '';
  for (let i = 0; i < source.length; i++) {
    const c = source[i];
    if (quote) {
      out += c;
      if (c === '\\' && i + 1 < source.length) out += source[++i];
      else if (c === quote || (c === '\n' && quote !== '`')) quote = '';
      continue;
    }
    if (c === '\'' || c === '"' || c === '`') {
      quote = c;
      out += c;
      continue;
    }
    if (c === '/' && source[i + 1] === '/') {
      const end = source.indexOf('\n', i);
      const stop = end < 0 ? source.length : end;
      out += ' '.repeat(stop - i);
      i = stop - 1;
      continue;
    }
    if (c === '/' && source[i + 1] === '*') {
      const end = source.indexOf('*/', i + 2);
      const stop = end < 0 ? source.length : end + 2;
      out += source.slice(i, stop).replace(/[^\n]/g, ' ');
      i = stop - 1;
      continue;
    }
    out += c;
  }
  return out;
}
