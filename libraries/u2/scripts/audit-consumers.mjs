/* The u2 consumer ratchet: per tree, the counts of the misuse patterns u2's consumers move off.
   node scripts/audit-consumers.mjs [--json] [--check] [--update-baseline [--allow-growth "<reason>"]] [--board <dir>]
   --check fails when a ratcheted count grew over consumers-baseline.json; --update-baseline refuses growth
   unless --allow-growth says why (the reason is kept in the baseline). */
import {existsSync, readdirSync, readFileSync, statSync, writeFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {join, parse, relative} from 'node:path';
import {fileURLToPath} from 'node:url';

const u2Root = fileURLToPath(new URL('..', import.meta.url));
const publicRoot = join(u2Root, '..', '..');
const baselineFile = join(u2Root, 'consumers-baseline.json');
const ts = createRequire(import.meta.url)('typescript');

const args = process.argv.slice(2);
const flag = (name) => args.includes(name);
const option = (name) => {
  const at = args.indexOf(name);
  return at >= 0 ? args[at + 1] : undefined;
};

const NEUTRAL = ['files', 'u2Files', 'cssFiles', 'controlBuild', 'commentLines'];
const SKIP_DIRS = new Set(['node_modules', 'dist', 'generated', 'test-results', '.test-build', '.turbo', 'playwright']);
const U2_IMPORT = /@datagrok-libraries\/u2|datagrok-api\/u2core|@datagrok-libraries\/ems|DG\.U2\b/;

const PATTERNS = {
  valueValue: /\.value\.value\b/g,
  valuePeek: /\.value\.peek\(\)/g,
  runInScope: /\brunInScope\(/g,
  runWith: /Scope\.runWith\(/g,
  runWithUndefined: /Scope\.runWith\(undefined/g,
  newScope: /new Scope\(/g,
  viewZoo: /\b(ownView|viewSignals|showInView|tableViewChrome|appView)\(/g,
  freeSubscribe: /(^|[^.\w])subscribe\(/gm,
  enabledOpt: /\benabled:/g,
  okEnabled: /\.okEnabled\(/g,
  deepDgIndex: /u2\/src\/dg\/index\.js/g,
  deepImports: /['"]@datagrok-libraries\/u2\/src\//g,
  controlBuild: /Control\.build\(/g,
  newAppHeader: /new AppHeader\(/g,
};
const CONSUMER_ONLY = {
  tsInternalSelectors: /(querySelector(All)?|closest|matches)\(\s*['"`][^'"`]*(\.u2-|\[data-u2)|classList\.(add|remove|toggle|contains)\(\s*['"`]u2-/g,
};

/** The trees audited: u2 itself (code patterns only), its engine, the companions, and every package that depends
 * on u2 or ems, plus ApiTests (the engine through DG.U2). */
function trees() {
  const out = [
    {key: 'libraries/u2/src', dirs: ['libraries/u2/src'], own: true},
    {key: 'libraries/u2/tests', dirs: ['libraries/u2/tests'], own: true},
    {key: 'libraries/u2/docs', dirs: ['libraries/u2/docs', 'libraries/u2/README.md'], own: true},
    {key: 'js-api/src/u2core', dirs: ['js-api/src/u2core'], own: true},
    {key: 'libraries/ems/src', dirs: ['libraries/ems/src'], css: ['libraries/ems/css']},
    {key: 'libraries/ems/tests', dirs: ['libraries/ems/tests']},
    {key: 'libraries/bdd', dirs: ['libraries/bdd/src', 'libraries/bdd/bindings']},
  ];
  const packages = join(publicRoot, 'packages');
  for (const name of readdirSync(packages).sort()) {
    const pkg = join(packages, name, 'package.json');
    if (!existsSync(pkg))
      continue;
    const json = JSON.parse(readFileSync(pkg, 'utf8'));
    const deps = {...json.dependencies, ...json.devDependencies};
    const uses = '@datagrok-libraries/u2' in deps || '@datagrok-libraries/ems' in deps || name === 'ApiTests';
    if (uses)
      out.push({key: `packages/${name}`, dirs: [`packages/${name}/src`, `packages/${name}/bdd/bindings`],
        css: [`packages/${name}/css`, `packages/${name}/src`]});
  }
  const board = option('--board');
  if (board)
    out.push({key: 'board-samples', dirs: [board], css: [board], external: true});
  return out;
}

function walk(path, accept, out = []) {
  if (!existsSync(path))
    return out;
  if (statSync(path).isFile()) {
    if (accept(path))
      out.push(path);
    return out;
  }
  for (const e of readdirSync(path, {withFileTypes: true}).sort((a, b) => a.name.localeCompare(b.name))) {
    if (e.isDirectory() && SKIP_DIRS.has(e.name))
      continue;
    walk(join(path, e.name), accept, out);
  }
  return out;
}

const isCode = (f) => /\.(ts|md)$/.test(f) && !/\.d\.ts$|\.g\.ts$/.test(f) || /\.test\.js$|\.mjs$/.test(f);
const isCss = (f) => f.endsWith('.css');
const resolve = (p) => p.startsWith('/') ? p : join(publicRoot, p);

/** File headers (a first block comment that is not a doc block) and doc blocks on private or protected members. */
function commentCounts(text) {
  const lines = text.split(/\r?\n/);
  let comment = 0;
  let header = 0;
  let privateDocs = 0;
  let inBlock = false;
  let first = true;
  for (let i = 0; i < lines.length; i++) {
    const t = lines[i].trim();
    if (t === '')
      continue;
    if (inBlock) {
      comment++;
      if (t.includes('*/'))
        inBlock = false;
      continue;
    }
    if (t.startsWith('/*')) {
      comment++;
      if (first && !t.startsWith('/**'))
        header++;
      if (t.startsWith('/**')) {
        let j = i;
        while (j < lines.length && !lines[j].includes('*/'))
          j++;
        if (/^(static\s+)?(private|protected)\b/.test((lines[j + 1] ?? '').trim()))
          privateDocs++;
      }
      if (!t.includes('*/'))
        inBlock = true;
      first = false;
      continue;
    }
    if (t.startsWith('//')) {
      comment++;
      continue;
    }
    first = false;
  }
  return {commentLines: comment, fileHeaders: header, privateDocs};
}

/** The u2 and ems classes that are controls (the closure under Component and Control) and the exported functions
 * whose declared return type is UI: what "builds u2 UI" means below. */
function uiNames() {
  const parents = new Map();
  const factories = new Set();
  for (const dir of ['libraries/u2/src', 'libraries/ems/src', 'js-api/src/u2core']) {
    for (const f of walk(resolve(dir), (p) => p.endsWith('.ts') && !p.endsWith('.d.ts'))) {
      const text = readFileSync(f, 'utf8');
      for (const m of text.matchAll(/class\s+(\w+)(?:<[^>{]*>)?\s+extends\s+(\w+)/g))
        parents.set(m[1], m[2]);
      for (const m of text.matchAll(/^export\s+(?:async\s+)?function\s+([a-z]\w*)[^\n{]*\)\s*:\s*([^{\n]+)\{/gm)) {
        if (/HTML|Element\b|Control|Input|Form|Dialog|Menu|View\b|Chrome|TypeAhead|Grid|Table|Tree|List/.test(m[2]))
          factories.add(m[1]);
      }
    }
  }
  const controls = new Set(['Component', 'Control']);
  for (let grew = true; grew;) {
    grew = false;
    for (const [child, parent] of parents) {
      if (controls.has(parent) && !controls.has(child)) {
        controls.add(child);
        grew = true;
      }
    }
  }
  return {controls, factories};
}

const UI = uiNames();

/** `rawSubscribe` is the count of `u2/no-raw-subscribe` findings, made by that rule, so the two cannot differ. */
function rawSubscribes() {
  const plugin = createRequire(join(publicRoot, 'libraries', 'eslint-plugin-u2', 'package.json'));
  const {Linter} = plugin('eslint');
  const linter = new Linter({configType: 'flat', cwd: parse(publicRoot).root});
  const config = [{files: ['**/*.ts', '**/*.js', '**/*.mjs'], plugins: {u2: plugin('./index.js')},
    languageOptions: {parser: plugin('@typescript-eslint/parser'), ecmaVersion: 'latest', sourceType: 'module'},
    linterOptions: {reportUnusedDisableDirectives: false}, rules: {'u2/no-raw-subscribe': 'warn'}}];
  return (file, text) => text.includes('subscribe') ?
    linter.verify(text, config, {filename: file}).filter((m) => m.ruleId === 'u2/no-raw-subscribe').length : 0;
}

const rawSubscribe = rawSubscribes();

/** Classes that extend a u2 or ems class (an identifier imported from either, or `DG.U2.<X>`), and other classes whose
 * constructor or field initializers build u2 UI (`new <a u2/ems control>(…)` or a u2/ems UI factory call). */
function classFacts(file, text) {
  const source = ts.createSourceFile(file, text, ts.ScriptTarget.Latest, true,
    file.endsWith('.ts') ? ts.ScriptKind.TS : ts.ScriptKind.JS);
  const imported = new Map();
  for (const st of source.statements) {
    if (!ts.isImportDeclaration(st) || !ts.isStringLiteral(st.moduleSpecifier))
      continue;
    const from = st.moduleSpecifier.text;
    const lib = /^@datagrok-libraries\/ems/.test(from) ? 'ems' :
      /^@datagrok-libraries\/u2|^datagrok-api\/u2core$/.test(from) ? 'u2' : null;
    const names = st.importClause?.namedBindings;
    if (lib === null || names === undefined || !ts.isNamedImports(names))
      continue;
    for (const el of names.elements)
      imported.set(el.name.text, {lib, name: (el.propertyName ?? el.name).text});
  }
  const found = {subclassU2: [], subclassEms: [], classBuildsUi: []};
  const libOf = (expr) => {
    if (ts.isIdentifier(expr))
      return imported.get(expr.text)?.lib ?? null;
    if (ts.isPropertyAccessExpression(expr) && ts.isPropertyAccessExpression(expr.expression) &&
      expr.expression.name.text === 'U2')
      return 'u2';
    return null;
  };
  const buildsUi = (node) => {
    let hit = false;
    const visit = (n) => {
      if (hit)
        return;
      // a callback runs later, except the builder handed to runInScope / build / runWith, which runs now
      const fn = ts.isArrowFunction(n) || ts.isFunctionExpression(n);
      const now = fn && ts.isCallExpression(n.parent) && /(^|\.)(runInScope|build|runWith)$/.test(n.parent.expression.getText());
      if (fn && n !== node && !now)
        return;
      const callee = (ts.isNewExpression(n) || ts.isCallExpression(n)) && ts.isIdentifier(n.expression) ?
        imported.get(n.expression.text) : undefined;
      if (callee !== undefined && (ts.isNewExpression(n) ? UI.controls.has(callee.name) : UI.factories.has(callee.name)))
        hit = true;
      else
        ts.forEachChild(n, visit);
    };
    visit(node);
    return hit;
  };
  const visit = (node) => {
    if (ts.isClassDeclaration(node) || ts.isClassExpression(node)) {
      const name = node.name?.text ?? '(anonymous)';
      const line = source.getLineAndCharacterOfPosition(node.getStart()).line + 1;
      const heritage = node.heritageClauses?.find((h) => h.token === ts.SyntaxKind.ExtendsKeyword)?.types[0];
      const lib = heritage === undefined ? null : libOf(heritage.expression);
      const where = `${name} (line ${line})`;
      if (lib === 'u2')
        found.subclassU2.push(where);
      else if (lib === 'ems')
        found.subclassEms.push(where);
      else {
        const parts = node.members.filter((m) => ts.isConstructorDeclaration(m) && m.body !== undefined ||
          ts.isPropertyDeclaration(m) && m.initializer !== undefined && !m.modifiers?.some((x) =>
            x.kind === ts.SyntaxKind.StaticKeyword));
        if (parts.some((m) => buildsUi(ts.isConstructorDeclaration(m) ? m.body : m.initializer)))
          found.classBuildsUi.push(where);
      }
    }
    ts.forEachChild(node, visit);
  };
  visit(source);
  return found;
}

function audit(tree) {
  const metrics = {files: 0, u2Files: 0};
  for (const k of [...Object.keys(PATTERNS), 'rawSubscribe', ...Object.keys(CONSUMER_ONLY), 'subclassU2', 'subclassEms', 'classBuildsUi',
    'cssFiles', 'cssLines', 'cssInternalLines', 'cssImportant', 'commentLines', 'fileHeaders', 'privateDocs'])
    metrics[k] = 0;
  const classes = {subclassU2: [], subclassEms: [], classBuildsUi: []};
  const codeFiles = tree.dirs.flatMap((d) => walk(resolve(d), isCode));
  for (const f of codeFiles) {
    const text = readFileSync(f, 'utf8');
    metrics.files++;
    const u2 = U2_IMPORT.test(text);
    if (u2)
      metrics.u2Files++;
    for (const [k, re] of Object.entries(PATTERNS))
      metrics[k] += (text.match(re) ?? []).length;
    if (f.endsWith('.md'))
      continue;
    const comments = commentCounts(text);
    metrics.commentLines += comments.commentLines;
    metrics.fileHeaders += comments.fileHeaders;
    metrics.privateDocs += comments.privateDocs;
    if (tree.own || !u2)
      continue;
    for (const [k, re] of Object.entries(CONSUMER_ONLY))
      metrics[k] += (text.match(re) ?? []).length;
    metrics.rawSubscribe += rawSubscribe(f, text);
    const facts = classFacts(f, text);
    const rel = tree.external ? f : relative(publicRoot, f);
    for (const k of Object.keys(classes)) {
      metrics[k] += facts[k].length;
      classes[k].push(...facts[k].map((c) => `${rel}: ${c}`));
    }
  }
  const cssFiles = [...new Set((tree.css ?? []).flatMap((d) => walk(resolve(d), isCss)))];
  for (const f of cssFiles) {
    const text = readFileSync(f, 'utf8');
    metrics.cssFiles++;
    metrics.cssLines += text.split('\n').length - (text.endsWith('\n') ? 1 : 0);
    metrics.cssInternalLines += text.split('\n').filter((l) => /\.u2-|\[data-u2/.test(l)).length;
    metrics.cssImportant += (text.match(/!important/g) ?? []).length;
  }
  return {metrics, classes};
}

const results = {};
for (const tree of trees())
  results[tree.key] = audit(tree);

function table() {
  const keys = Object.keys(Object.values(results)[0].metrics);
  const rows = Object.entries(results).map(([k, r]) => [k, ...keys.map((m) => String(r.metrics[m]))]);
  const head = ['tree', ...keys];
  const widths = head.map((h, i) => Math.max(h.length, ...rows.map((r) => r[i].length)));
  const line = (cells) => cells.map((c, i) => i === 0 ? c.padEnd(widths[i]) : c.padStart(widths[i])).join('  ');
  console.log(line(head));
  for (const r of rows)
    console.log(line(r));
  for (const [k, r] of Object.entries(results)) {
    for (const [what, list] of Object.entries(r.classes))
      for (const c of list)
        console.log(`  ${k} ${what}: ${c}`);
  }
}

function snapshot() {
  const out = {};
  for (const [k, r] of Object.entries(results)) {
    if (k === 'board-samples')
      continue;
    // check-consumer's lint counts ride along unchanged; only check-consumer --update-baseline rewrites them
    out[k] = {...r.metrics, classes: r.classes, ...baseline?.trees[k]?.lint ? {lint: baseline.trees[k].lint} : {}};
  }
  return out;
}

function grown(baseline) {
  const out = [];
  for (const [k, r] of Object.entries(results)) {
    const before = baseline.trees[k];
    if (k === 'board-samples')
      continue;
    for (const [m, v] of Object.entries(r.metrics)) {
      if (NEUTRAL.includes(m))
        continue;
      const was = before?.[m] ?? 0;
      if (v > was)
        out.push(`${k} ${m}: ${was} → ${v}`);
    }
  }
  return out;
}

const baseline = existsSync(baselineFile) ? JSON.parse(readFileSync(baselineFile, 'utf8')) : null;
if (flag('--json'))
  console.log(JSON.stringify(results, null, 2));
else
  table();

if (flag('--check')) {
  if (baseline === null) {
    console.error('audit-consumers: no consumers-baseline.json; run with --update-baseline first');
    process.exit(1);
  }
  const growth = grown(baseline);
  if (growth.length > 0) {
    console.error(`audit-consumers: ${growth.length} count(s) grew over the baseline:\n  ${growth.join('\n  ')}\n` +
      'Remove the new uses, or run --update-baseline --allow-growth "<why>" in the same change.');
    process.exit(1);
  }
}

if (flag('--update-baseline')) {
  const growth = baseline === null ? [] : grown(baseline);
  const reason = option('--allow-growth');
  if (growth.length > 0 && !reason) {
    console.error(`audit-consumers: refusing to raise the baseline:\n  ${growth.join('\n  ')}\n` +
      'Pass --allow-growth "<why>" to record the reason.');
    process.exit(1);
  }
  const log = [...baseline?.growth ?? [], ...growth.length > 0 ? [{date: new Date().toISOString().slice(0, 10), reason,
    grown: growth}] : []];
  const out = {
    description: 'Shrink-only counts per u2 consumer (scripts/audit-consumers.mjs). Neutral, not ratcheted: ' +
      NEUTRAL.join(', ') + '.',
    trees: snapshot(),
    growth: log,
  };
  writeFileSync(baselineFile, JSON.stringify(out, null, 2) + '\n');
  console.log(`audit-consumers: wrote ${relative(publicRoot, baselineFile)}`);
}
