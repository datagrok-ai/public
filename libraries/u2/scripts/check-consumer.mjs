import {spawnSync} from 'node:child_process';
import {existsSync, readdirSync, readFileSync, realpathSync, statSync, writeFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {dirname, isAbsolute, join, relative, resolve, sep} from 'node:path';
import {fileURLToPath} from 'node:url';

const u2Root = realpathSync(fileURLToPath(new URL('..', import.meta.url)));
const publicRoot = join(u2Root, '..', '..');
const pluginDir = join(publicRoot, 'libraries', 'eslint-plugin-u2');
const auditScript = join(u2Root, 'scripts', 'audit-consumers.mjs');
const baselineFile = join(u2Root, 'consumers-baseline.json');
const {INTERNAL_SELECTOR} = createRequire(join(pluginDir, 'package.json'))('./lib/u2-imports.js');

const USAGE = 'node scripts/check-consumer.mjs <package-dir>... | --all  [--json] [--verbose] ' +
  '[--update-baseline [--allow-growth "<why>"]]';

const args = process.argv.slice(2);
const flag = (name) => args.includes(name);
const option = (name) => {
  const at = args.indexOf(name);
  return at >= 0 ? args[at + 1] : undefined;
};
const targets = args.filter((a, i) => !a.startsWith('--') && args[i - 1] !== '--allow-growth');

/** Audit metrics each lint rule points at: a grown metric lists that rule's findings. */
const RULE_OF = {
  deepImports: 'u2/import-entry', tsInternalSelectors: 'u2/no-internal-selectors', rawSubscribe: 'u2/no-raw-subscribe',
  runWithUndefined: 'u2/no-unowned-root', newScope: 'u2/no-unowned-root', subclassU2: 'u2/no-subclass',
  subclassEms: 'u2/no-subclass', classBuildsUi: 'u2/no-ui-in-class', cssInternalLines: 'css-internal-selectors',
  fileHeaders: 'u2/comments', privateDocs: 'u2/comments',
};
const SHOWN = ['subclassU2', 'subclassEms', 'classBuildsUi', 'deepImports', 'tsInternalSelectors', 'rawSubscribe',
  'runWithUndefined', 'cssInternalLines'];
/** Report-mode rules whose finding count is itself ratcheted (`lint` of a tree in consumers-baseline.json). */
const LINT_RATCHET = ['u2/comments', 'u2/no-snapshot-in-live-slot', 'u2/reason-text', 'u2/no-raw-column-list',
  'u2/build-in-owner', 'u2/no-unowned-root'];
const CLASS_METRICS = ['subclassU2', 'subclassEms', 'classBuildsUi'];
const ALLOW = {'u2/no-subclass': ['packages/ApiTests/src/u2/u2.ts']};
const SKIP_DIRS = new Set(['node_modules', 'dist', 'generated', 'test-results', '.test-build', '.turbo', 'playwright']);
const UI_TYPE = /HTML|Element\b|Control|Input|Form|Dialog|Menu|View\b|Chrome|TypeAhead|Grid|Table|Tree|List/;

if (targets.length === 0 && !flag('--all') && !flag('--update-baseline')) {
  console.error(`check-consumer: name a package directory or pass --all\n  ${USAGE}`);
  process.exit(2);
}

const posix = (p) => p.split(sep).join('/');
const rel = (f) => {
  const r = posix(relative(publicRoot, f));
  return r.startsWith('../') ? posix(f) : r;
};

function walk(path, accept, out = []) {
  if (!existsSync(path))
    return out;
  if (statSync(path).isFile()) {
    if (accept(path))
      out.push(path);
    return out;
  }
  for (const e of readdirSync(path, {withFileTypes: true}).sort((a, b) => a.name.localeCompare(b.name))) {
    if (!(e.isDirectory() && SKIP_DIRS.has(e.name)))
      walk(join(path, e.name), accept, out);
  }
  return out;
}

// Two runs: the audit exits at once on growth, which can cut a long --json write to a pipe short.
function audit() {
  const run = (flagName) => spawnSync(process.execPath, [auditScript, flagName],
    {encoding: 'utf8', maxBuffer: 256 << 20});
  const json = run('--json');
  const check = run('--check');
  if (json.status !== 0 || check.status !== 0 && !/count\(s\) grew/.test(check.stderr))
    throw new Error(`check-consumer: the audit failed: ${(json.stderr + check.stderr).trim()}`);
  const grown = [...check.stderr.matchAll(/^\s+(\S+) (\w+): (\d+) → (\d+)$/gm)]
    .map((m) => ({tree: m[1], metric: m[2], baseline: Number(m[3]), now: Number(m[4])}));
  return {results: JSON.parse(json.stdout), grown};
}

/** The u2 and ems control classes and UI factories (the audit's method; type parameters may nest one level). */
function libraryNames() {
  const parents = new Map();
  const factories = new Set();
  const controlFactories = new Set();
  const returnTypes = new Map();
  const roots = new Set();
  for (const dir of ['libraries/u2/src', 'libraries/ems/src', 'js-api/src/u2core']) {
    for (const f of walk(join(publicRoot, dir), (p) => p.endsWith('.ts') && !p.endsWith('.d.ts'))) {
      const text = readFileSync(f, 'utf8');
      for (const m of text.matchAll(/class\s+(\w+)\s*(?:<(?:[^<>{]|<[^<>{]*>)*>)?\s+extends\s+(\w+)/g))
        parents.set(m[1], m[2]);
      let cls = null;
      for (const m of text.matchAll(/\bclass\s+(\w+)|\bstatic\s+isRoot\s*=\s*true/g)) {
        if (m[1] !== undefined)
          cls = m[1];
        else if (cls !== null)
          roots.add(cls);
      }
      for (const m of text.matchAll(/^export\s+(?:async\s+)?function\s+([a-z]\w*)[^\n{]*\)\s*:\s*([^{\n]+)\{/gm)) {
        if (!UI_TYPE.test(m[2]))
          continue;
        factories.add(m[1]);
        if (!/HTML|Element\b/.test(m[2]))
          controlFactories.add(m[1]);
        returnTypes.set(m[1], m[2]);
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
      if (roots.has(parent) && !roots.has(child)) {
        roots.add(child);
        grew = true;
      }
    }
  }
  // the factories that build a control (`funcForm(…): FuncCallForm`), not a status (`state(…): IControlState`) or a
  // lookup (`controlAt(…): Control | undefined`)
  const builders = [...controlFactories].filter((name) => {
    const types = returnTypes.get(name).match(/\w+/g) ?? [];
    return !types.includes('undefined') && !types.includes('null') && types.some((t) => controls.has(t) && !roots.has(t));
  });
  return {controls: [...controls], roots: [...roots], factories: [...factories],
    controlFactories: [...controlFactories], builders};
}

/** The slots u2 and ems keep live, read off their declarations with the type checker: per exported function, class
 * constructor and method, the parameters typed as a signal (or a union holding one, as `Text` and `Child` do) and the
 * options of that kind an options parameter takes (`value`, an initial value, aside). */
function liveSlots(ts) {
  const dirs = ['libraries/u2/src', 'libraries/ems/src'].map((d) => join(publicRoot, d) + sep);
  const files = dirs.flatMap((d) => walk(d, (p) => p.endsWith('.ts') && !p.endsWith('.d.ts')));
  const program = ts.createProgram(files, {noEmit: true, skipLibCheck: true, types: [], target: ts.ScriptTarget.ES2022,
    module: ts.ModuleKind.ESNext, moduleResolution: ts.ModuleResolutionKind.Bundler});
  const checker = program.getTypeChecker();
  const LIVE = new Set(['Signal', 'ReadonlySignal', 'Computed']);
  const own = (symbol) => symbol.declarations?.some((d) =>
    dirs.some((dir) => d.getSourceFile().fileName.startsWith(dir)));
  const typeOf = (symbol, at) => checker.getNonNullableType(checker.getTypeOfSymbolAtLocation(symbol, at));
  const live = (t, depth = 0) => depth < 4 && (LIVE.has(t.aliasSymbol?.name) || LIVE.has(t.getSymbol()?.name) ||
    t.isUnion() && t.types.some((x) => live(x, depth + 1)) ||
    checker.isArrayType(t) && checker.getTypeArguments(t).some((x) => live(x, depth + 1)));
  const slots = (signatures, at) => {
    const out = {at: [], options: {}};
    for (const sig of signatures) {
      sig.parameters.forEach((p, i) => {
        const t = typeOf(p, at);
        if (live(t))
          out.at.push(i);
        else if (t.getCallSignatures().length === 0) {
          const names = t.getProperties().filter((q) => q.name !== 'value' && own(q) && live(typeOf(q, at)))
            .map((q) => q.name);
          if (names.length > 0)
            out.options[i] = [...new Set([...out.options[i] ?? [], ...names])];
        }
      });
    }
    out.at = [...new Set(out.at)];
    return out.at.length > 0 || Object.keys(out.options).length > 0 ? out : undefined;
  };
  const functions = {};
  const classes = {};
  for (const sf of program.getSourceFiles().filter((f) => dirs.some((d) => f.fileName.startsWith(d)))) {
    const module = checker.getSymbolAtLocation(sf);
    for (const exported of module ? checker.getExportsOfModule(module) : []) {
      const symbol = exported.flags & ts.SymbolFlags.Alias ? checker.getAliasedSymbol(exported) : exported;
      if (!own(symbol))
        continue;
      if (symbol.flags & ts.SymbolFlags.Function) {
        const sigs = typeOf(symbol, sf).getCallSignatures();
        const entry = slots(sigs, sf) ?? {};
        const returned = sigs[0] && checker.getReturnTypeOfSignature(sigs[0]).getSymbol();
        if (returned && returned.flags & ts.SymbolFlags.Class && own(returned))
          entry.returns = returned.name;
        if (Object.keys(entry).length > 0)
          functions[exported.name] = entry;
      }
      else if (symbol.flags & ts.SymbolFlags.Class && !classes[exported.name]) {
        const entry = {methods: {}};
        const ctor = slots(typeOf(symbol, sf).getConstructSignatures(), sf);
        if (ctor)
          entry.ctor = ctor;
        const instance = checker.getDeclaredTypeOfSymbol(symbol);
        for (const member of instance.getProperties().filter(own)) {
          if (member.name === 'value')
            entry.valueAccessor = member.declarations.some((d) => d.kind === ts.SyntaxKind.GetAccessor);
          else if (member.flags & ts.SymbolFlags.Method) {
            const method = slots(typeOf(member, sf).getCallSignatures(), sf);
            if (method)
              entry.methods[member.name] = method;
          }
        }
        classes[exported.name] = entry;
      }
    }
  }
  return {functions, classes};
}

function lintScope(key) {
  if (key.startsWith('packages/') || isAbsolute(key))
    return {kind: 'consumer', dirs: [`${key}/src`, `${key}/bdd/bindings`], css: [`${key}/css`, `${key}/src`]};
  if (key === 'libraries/u2/src' || key === 'libraries/ems/src')
    return {kind: 'library', dirs: [key], css: []};
  return null;
}

function linter() {
  const pluginRequire = createRequire(join(pluginDir, 'package.json'));
  const plugin = pluginRequire('./index.js');
  const {Linter} = pluginRequire('eslint');
  const parser = pluginRequire('@typescript-eslint/parser');
  const {U2Imports} = pluginRequire('./lib/u2-imports.js');
  const names = libraryNames();
  const live = liveSlots(pluginRequire('typescript'));
  const configs = {};
  for (const kind of ['consumer', 'library']) {
    const rules = {};
    for (const [rule, severity] of Object.entries(plugin.configs[kind].rules)) {
      rules[rule] = rule === 'u2/no-ui-in-class' ? [severity, {controls: names.controls, factories: names.factories}] :
        rule === 'u2/build-in-owner' ? [severity, {controls: names.controls, factories: names.controlFactories,
          methods: kind === 'consumer'}] :
          // a library's controls build under their own scopes, which build-in-owner holds
          rule === 'u2/no-unowned-root' ? [severity, kind === 'consumer' ? {controls: names.controls,
            roots: names.roots, factories: names.builders} : {controls: [], factories: []}] :
          rule === 'u2/no-snapshot-in-live-slot' ? [severity, live] : severity;
    }
    configs[kind] = [{
      files: ['**/*.ts', '**/*.tsx', '**/*.js', '**/*.mjs'],
      languageOptions: {parser, ecmaVersion: 'latest', sourceType: 'module'},
      linterOptions: {reportUnusedDisableDirectives: false},
      plugins: plugin.configs[kind].plugins,
      rules,
    }];
  }
  // a flat config matches no file outside its linter's folder: a tree outside public/ is linted from its own folder
  const engines = new Map();
  const engineFor = (file) => {
    const cwd = file.startsWith(publicRoot + sep) ? publicRoot : dirname(file);
    if (!engines.has(cwd))
      engines.set(cwd, new Linter({configType: 'flat', cwd}));
    return engines.get(cwd);
  };
  return (file, kind) => {
    const text = readFileSync(file, 'utf8');
    if (kind === 'consumer' && !U2Imports.usesU2(text))
      return [];
    return engineFor(file).verify(text, configs[kind], {filename: file})
      .filter((m) => m.fatal || m.ruleId?.startsWith('u2/'))
      .filter((m) => !ALLOW[m.ruleId]?.includes(rel(file)))
      .map((m) => ({rule: m.fatal ? 'parse' : m.ruleId, mode: m.severity === 2 ? 'error' : 'report', file: rel(file),
        line: m.line, column: m.column, message: m.message}));
  };
}

const lintable = (f) => /\.(ts|tsx|mjs)$/.test(f) && !/\.d\.ts$|\.g\.ts$/.test(f) ||
  f.endsWith('.js') && !existsSync(f.slice(0, -3) + '.ts');

/** The lines of a tree's CSS that name a u2 class or `data-u2` attribute: notes while their count stays at what
 * consumers-baseline.json allows the tree (`cssInternalLines`; none for a package it does not list), errors once it
 * is more, so a new internal selector fails the check. */
function cssFindings(key, dirs, baseline) {
  const out = [];
  for (const f of [...new Set(dirs.flatMap((d) => walk(resolve(publicRoot, d), (p) => p.endsWith('.css'))))]) {
    const lines = readFileSync(f, 'utf8').split('\n');
    for (let i = 0; i < lines.length; i++) {
      if (INTERNAL_SELECTOR.test(lines[i]))
        out.push({rule: 'css-internal-selectors', mode: 'report', file: rel(f), line: i + 1, column: 1,
          message: `\`${lines[i].trim()}\` styles u2 internals: style your own classes and set \`--u2-*\` ` +
            'properties.'});
    }
  }
  const allowed = baseline.trees[key]?.cssInternalLines ?? 0;
  if (out.length > allowed) {
    for (const f of out) {
      f.mode = 'error';
      f.message += ` (${out.length} in this tree's CSS; consumers-baseline.json allows ${allowed})`;
    }
  }
  return out;
}

const usesU2 = (dir) => {
  const pkg = join(dir, 'package.json');
  if (!existsSync(pkg))
    return false;
  const json = JSON.parse(readFileSync(pkg, 'utf8'));
  const deps = {...json.dependencies, ...json.devDependencies};
  return '@datagrok-libraries/u2' in deps || '@datagrok-libraries/ems' in deps;
};

/** The audited trees under `dir`; a package outside this repository's packages/ that depends on u2 or ems is a tree
 * of its own, keyed by its absolute path, linted against no baseline. */
function treeKeys(dir, results) {
  if (!existsSync(dir))
    throw new Error(`check-consumer: ${dir} does not exist`);
  const real = realpathSync(resolve(dir));
  const at = posix(relative(publicRoot, real));
  const keys = Object.keys(results).filter((k) => k === at || k.startsWith(at + '/'));
  if (keys.length > 0)
    return keys;
  if (posix(relative(join(publicRoot, 'packages'), real)).startsWith('../') && usesU2(real))
    return [posix(real)];
  throw new Error(`check-consumer: ${dir} is not a tree the u2 audit knows (a package that depends on u2 or ems)`);
}

const readBaseline = () => existsSync(baselineFile) ? JSON.parse(readFileSync(baselineFile, 'utf8')) : {trees: {}};

/** Per tree with a lint scope, the count of each ratcheted rule's findings. */
function lintCounts(findings) {
  return Object.fromEntries(LINT_RATCHET.map((rule) => [rule, findings.filter((f) => f.rule === rule).length]));
}

/** A package's count never recorded is 0, so a new package starts clean; elsewhere a count never recorded is not
 * growth, and its first `--update-baseline` records it. */
function lintGrowth(key, counts, baseline) {
  const unrecorded = key.startsWith('packages/') ? 0 : undefined;
  return Object.entries(counts).map(([rule, now]) => ({tree: key, metric: rule,
    baseline: baseline.trees[key]?.lint?.[rule] ?? unrecorded, now}))
    .filter((g) => g.baseline !== undefined && g.now > g.baseline);
}

function check(baseline) {
  const {results, grown} = audit();
  const keys = flag('--all') || flag('--update-baseline') ? Object.keys(results) :
    [...new Set(targets.flatMap((d) => treeKeys(d, results)))];
  const lint = linter();
  const trees = [];
  for (const key of keys) {
    const scope = lintScope(key);
    const findings = scope === null ? [] : [
      ...scope.dirs.flatMap((d) => walk(resolve(publicRoot, d), lintable)).flatMap((f) => lint(f, scope.kind)),
      ...cssFindings(key, scope.css, baseline),
    ];
    const lintNow = scope === null ? null : lintCounts(findings);
    const growth = [...grown.filter((g) => g.tree === key),
      ...lintNow === null ? [] : lintGrowth(key, lintNow, baseline)];
    const counts = results[key] === undefined ? {} : Object.fromEntries(SHOWN.map((m) => [m,
      {baseline: baseline.trees[key]?.[m] ?? 0, now: results[key].metrics[m]}]));
    for (const [rule, now] of Object.entries(lintNow ?? {}))
      counts[rule] = {baseline: baseline.trees[key]?.lint?.[rule] ?? '–', now};
    trees.push({tree: key, kind: scope?.kind ?? null, grown: growth, counts, findings, lint: lintNow,
      ...describe(key, scope?.kind ?? null, growth, counts, findings)});
  }
  return trees;
}

/** The audit writes its counts (refusing their growth without a reason), then the lint counts join each tree; lint
 * growth needs the same reason and is logged with it. */
function updateBaseline() {
  const reason = option('--allow-growth');
  const before = readBaseline();
  const trees = check(before);
  const growth = trees.flatMap((t) => t.grown.filter((g) => LINT_RATCHET.includes(g.metric)));
  if (growth.length > 0 && !reason) {
    const lines = growth.map((g) => `${g.tree} ${g.metric}: ${g.baseline} → ${g.now}`).join('\n  ');
    console.error(`check-consumer: refusing to raise the baseline:\n  ${lines}\n` +
      'Pass --allow-growth "<why>" to record the reason.');
    return 1;
  }
  const run = spawnSync(process.execPath,
    [auditScript, '--update-baseline', ...reason ? ['--allow-growth', reason] : []], {stdio: 'inherit'});
  if (run.status !== 0)
    return run.status ?? 1;
  const after = readBaseline();
  for (const t of trees) {
    if (t.lint !== null && after.trees[t.tree])
      after.trees[t.tree].lint = t.lint;
  }
  if (growth.length > 0) {
    after.growth.push({date: new Date().toISOString().slice(0, 10), reason,
      grown: growth.map((g) => `${g.tree} ${g.metric}: ${g.baseline} → ${g.now}`)});
  }
  writeFileSync(baselineFile, JSON.stringify(after, null, 2) + '\n');
  console.log(`check-consumer: lint counts recorded in ${rel(baselineFile)}`);
  return 0;
}

function byRule(findings) {
  const out = new Map();
  for (const f of findings)
    out.set(f.rule, (out.get(f.rule) ?? 0) + 1);
  return [...out].map(([rule, n]) => `${rule} ${n}`).join(', ');
}

/** The lines a reader gets: `summary` always, `problems` fail the check, `notes` are the report-mode findings. */
function describe(key, kind, grown, counts, findings) {
  const name = key.startsWith('packages/') ? key.slice('packages/'.length) : isAbsolute(key) ? key.split('/').pop() :
    key;
  const at = (f) => `${f.file}:${f.line}:${f.column}  ${f.rule}  ${f.message}`;
  const shown = Object.entries(counts).filter(([, c]) => (Number(c.baseline) || 0) + c.now > 0)
    .map(([m, c]) => `${m} ${c.now}/${c.baseline}`).join(', ');
  const report = findings.filter((f) => f.mode === 'report');
  const errors = findings.filter((f) => f.mode === 'error');
  const summary = [isAbsolute(key) ? `u2 ratchet ${name}: not in consumers-baseline.json (outside packages/)` :
    grown.length === 0 ? `u2 ratchet ${name}: ok, at or under consumers-baseline.json${shown ? ` (${shown})` : ''}` :
      `u2 ratchet ${name}: ${grown.length} count(s) grew over consumers-baseline.json`];
  if (kind !== null) {
    summary.push(`u2 lint ${name}: ${errors.length ? `${errors.length} error(s) (${byRule(errors)})` : 'no errors'}; ` +
      `report mode: ${report.length ? byRule(report) : 'none'}`);
  }
  const problems = [...errors.map(at), ...grown.map((g) => {
    const rule = RULE_OF[g.metric] ?? g.metric;
    const where = report.filter((f) => f.rule === rule).map((f) => `\n    ${at(f)}`).join('');
    const remedy = CLASS_METRICS.includes(g.metric) ?
      'Function components are the documented form: rewrite the class as one, or, where a class is truly needed, ' +
      'record why' : 'Remove the new uses, or record why';
    return `${name}: ${g.metric} grew ${g.baseline} → ${g.now} over consumers-baseline.json. ${remedy} in the same ` +
      `change: node libraries/u2/scripts/check-consumer.mjs --update-baseline --allow-growth "<why>"${where}`;
  })];
  return {ok: problems.length === 0, summary, problems, notes: report.map(at)};
}

try {
  if (flag('--update-baseline'))
    process.exit(updateBaseline());
  const trees = check(readBaseline());
  if (flag('--json'))
    console.log(JSON.stringify({ok: trees.every((t) => t.ok), trees}, null, 2));
  else {
    for (const t of trees) {
      for (const line of [...t.summary, ...t.problems.map((p) => `  ${p}`), ...flag('--verbose') ?
        t.notes.map((n) => `  ${n}`) : []])
        console.log(line);
    }
  }
  process.exitCode = trees.every((t) => t.ok) ? 0 : 1;
}
catch (e) {
  console.error(e.message);
  process.exitCode = 2;
}
