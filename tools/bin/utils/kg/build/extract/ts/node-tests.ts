/// vitest, jest and node:test cases (change-tests/plan.md): `describe`/`test`/`it` titles in the `*.test.ts` and
/// `*.spec.ts` sources the ts pass parses (the CLI, the libraries, the packages), one suite per file, run by
/// `npx vitest run <file>`. A spec the Playwright or the DG runner owns is left to its own extractor.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../../emitter';
import {BuildContext, Extractor} from '../../context';
import {fileId, testId, suiteId} from '../../../ids';
import {tsSources} from './declarations';
import {blankComments, callbackBrace, matchBrace} from './tests';

export interface NodeTest {
  /** Enclosing `describe` titles, outermost first. */
  describes: string[];
  title: string;
  skipped: boolean;
  /** The runner fills the title in (a `.each` table, a template, a concatenation); [title] is its literal head with an ellipsis. */
  dynamic: boolean;
}

const TEST_FILE = /\.(test|spec)\.(?:tsx?|[cm]?js)$/;
const OTHER_RUNNERS = /^(@playwright\/test|playwright(-core)?|@datagrok-libraries\/utils\/src\/test)$/;
/** `.each([...])` or `.each\`table\`` between the name and the title, one level of parentheses inside. */
const EACH = String.raw`\.each\s*(?:\((?:[^()]|\([^()]*\))*\)|\`[^\`]*\`)`;
const nodeCall = (names: string[]) => new RegExp(String.raw`(?<![\w.$])((?:${names.join('|')})(?:\.(?:only|skip|todo|concurrent|sequential|shuffle|fails))*(?:${EACH})?)\s*\(\s*(['"\`])((?:\\.|(?!\2).)*)\2`, 'g');
const CALLS = ['describe', 'test', 'it'];
/** A local `function smoke(name, body) { test(name, …` (the u2 suites): its calls register tests under another spelling. */
const WRAPPER = /(?<![\w.$])function\s+([A-Za-z_$][\w$]*)\s*\(\s*([A-Za-z_$][\w$]*)\b[^)]*\)\s*\{\s*(?:test|it)\s*\(\s*\2\s*,/g;
const CONCATENATED = /^\s*\+/;
const INTERPOLATION = '${';
/** A `%s` printf or a `$name` placeholder an `.each` title carries. */
const PLACEHOLDER = /[%$]/;

export const nodeTestsExtractor: Extractor = {
  name: 'ts-node-tests',
  describes: {'ts-node-tests': 'vitest, jest and node:test cases and their suites'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    for (const file of tsSources(ctx, emitter).files) {
      if (!TEST_FILE.test(file.path) || file.imports.some((i) => OTHER_RUNNERS.test(i.specifier))) continue;
      const tests = parseNodeTests(fs.readFileSync(path.join(ctx.repoRoot, file.path), 'utf8'));
      if (!tests.length) continue;
      const suite = suiteId('node', file.path);
      emitter.node({type: 'test-suite', id: suite, name: path.posix.basename(file.path), framework: 'node', path: file.path,
        package: file.unit.id.startsWith('pkg:') ? file.unit.id : undefined, provenance: 'filesystem', source_layer: 'public'});
      for (const t of tests) {
        const category = t.describes.join(' > ') || undefined;
        const id = testId('node', file.path, category ?? '', t.title);
        if (!emitter.node({type: 'test', id, name: t.title, path: file.path, framework: 'node', level: 'unit', category, suite,
          dynamic: t.dynamic ? true : undefined, skipped: t.skipped ? true : undefined, provenance: 'ast', source_layer: 'public'}).accepted) continue;
        if (emitter.has(fileId(file.path))) emitter.edge({type: 'declares', from: fileId(file.path), to: id, derived_by: 'ast', confidence: 1, evidence: [file.path]});
        if (t.dynamic) emitter.problem('dynamic_tests', `${file.path}: ${category ? `${category}/` : ''}${t.title} names a registration site, not a runnable test: the title is built at run time`);
      }
    }
  },
};

/** The `test`/`it` calls of a suite with the `describe` chain around each; `.skip` and `.todo` on the call or a describe make a
 * test skipped, `.each` on either makes it dynamic; a commented-out test is not a test. */
export function parseNodeTests(source: string): NodeTest[] {
  const text = blankComments(source);
  const describes: {title: string, start: number, end: number, skipped: boolean, dynamic: boolean}[] = [];
  const out: NodeTest[] = [];
  const wrappers = [...text.matchAll(WRAPPER)].map((m) => m[1]).filter((n) => !CALLS.includes(n));
  for (const m of text.matchAll(nodeCall([...CALLS, ...wrappers]))) {
    const call = m[1];
    const at = m.index!;
    const skipped = /\.(skip|todo)\b/.test(call);
    const {title, dynamic} = titleOf(call, m[2], m[3], text.slice(at + m[0].length));
    if (call.startsWith('describe')) {
      const open = callbackBrace(text, at + m[0].length);
      if (open >= 0) describes.push({title, start: open, end: matchBrace(text, open), skipped, dynamic});
      continue;
    }
    const enclosing = describes.filter((d) => d.start < at && at < d.end);
    out.push({describes: enclosing.map((d) => d.title), title, skipped: skipped || enclosing.some((d) => d.skipped), dynamic: dynamic || enclosing.some((d) => d.dynamic)});
  }
  return out;
}

/** The title as written, or its literal head with an ellipsis when the runner builds it (the DG rule, plus `.each` placeholders). */
function titleOf(call: string, quote: string, raw: string, after: string): {title: string, dynamic: boolean} {
  const each = call.includes('.each');
  const interpolated = quote === '`' && raw.includes(INTERPOLATION);
  if (!each && !interpolated && !CONCATENATED.test(after)) return {title: raw.trim(), dynamic: false};
  const cut = each ? raw.search(PLACEHOLDER) : interpolated ? raw.indexOf(INTERPOLATION) : -1;
  return {title: `${(cut < 0 ? raw : raw.slice(0, cut)).trim()}…`, dynamic: true};
}
