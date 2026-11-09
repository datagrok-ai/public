import {spawn} from 'node:child_process';
import {copyFileSync, existsSync, lstatSync, mkdirSync, readFileSync, realpathSync, renameSync, rmdirSync, rmSync,
  statSync, symlinkSync, unlinkSync, writeFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {basename, dirname, join, relative, resolve, sep} from 'node:path';
import {compileFeature, CompiledFeature, Diagnostic, outFileFor} from './compile.js';
import {listFiles, loadBindings} from './discover.js';
import {GherkinParseError, parseFeature} from './gherkin.js';
import {scaffold} from './init.js';
import {StepMatcher} from './match.js';
import {availableTiers, LIB_DIR, LIB_ROOT, loadProject, PACKAGE_NAME, Project} from './project.js';
import {allContexts, allDatasets, allElements, allKinds, steps} from './registry.js';
import {slugOf} from './runtime/guide.js';

function rel(project: Project, file: string): string {
  return relative(project.root, file).split(sep).join('/');
}

async function compileAll(project: Project): Promise<{compiled: CompiledFeature[]; diagnostics: Diagnostic[]}> {
  const bindings = await loadBindings(project.sources);
  const matcher = new StepMatcher();
  const compiled: CompiledFeature[] = [];
  const diagnostics: Diagnostic[] = [];
  for (const file of listFiles(project.featuresDir, '.feature')) {
    let feature;
    try {
      feature = parseFeature(file, readFileSync(file, 'utf8'));
    } catch (e) {
      if (e instanceof GherkinParseError) {
        diagnostics.push({file: rel(project, file), line: 0, level: 'error', message: e.message});
        continue;
      }
      throw e;
    }
    const result = compileFeature(feature, {root: project.root, matcher, bindings});
    compiled.push(result);
    diagnostics.push(...result.diagnostics);
  }
  return {compiled, diagnostics};
}

/** Errors and warnings always; the notes (how each phrase resolves) only when asked — they are an
 * authoring aid, and hundreds of them bury the one line that matters in a run. */
function report(diagnostics: Diagnostic[], notes: boolean): number {
  for (const d of diagnostics) {
    if (notes || d.level !== 'note')
      console.log(`${d.file}:${d.line} ${d.level}: ${d.message}`);
  }
  const errors = diagnostics.filter((d) => d.level === 'error').length;
  const noted = diagnostics.length - errors;
  console.log(`${errors} error(s), ${noted} note(s)${noted > 0 && !notes ? ' (--verbose shows them)' : ''}`);
  return errors;
}

const eol = (s: string) => s.replace(/\r\n/g, '\n');

/** `only`: the features the command is about — only their errors fail it, and a check reads only their specs (a
 * guide films a few, a run names its specs; another feature being written beside them does not stop either). */
async function compile(project: Project, check: boolean, notes = false, only?: Set<string>): Promise<number> {
  const {compiled, diagnostics} = await compileAll(project);
  const counted = report(diagnostics, notes);
  const errors = only ? diagnostics.filter((d) => d.level === 'error' && only.has(resolve(project.root, d.file))).length : counted;
  const expected = new Set(compiled.map((c) => c.outFile));
  const orphans = only ? [] : listFiles(project.generatedDir, '.test.ts').filter((f) => !expected.has(f));
  let drift = 0;
  for (const c of compiled) {
    if (check && only && !only.has(c.feature.path))
      continue;
    const current = existsSync(c.outFile) ? eol(readFileSync(c.outFile, 'utf8')) : undefined;
    if (current === eol(c.code))
      continue;
    drift++;
    if (check) {
      console.log(`${current === undefined ? 'missing' : 'stale'}: ${rel(project, c.outFile)}`);
      continue;
    }
    mkdirSync(dirname(c.outFile), {recursive: true});
    writeFileSync(c.outFile, c.code, 'utf8');
    console.log(`wrote ${rel(project, c.outFile)}`);
  }
  for (const f of orphans) {
    drift++;
    if (check)
      console.log(`orphan: ${rel(project, f)} (no feature produces it)`);
    else {
      rmSync(f);
      console.log(`removed ${rel(project, f)}`);
    }
  }
  if (check) {
    console.log(drift === 0 ? `${compiled.length} spec(s) up to date` : `${drift} spec(s) out of date — run \`grok-bdd compile\``);
    return errors > 0 || drift > 0 ? 1 : 0;
  }
  console.log(`${compiled.length} spec(s), ${drift} written`);
  return errors > 0 ? 1 : 0;
}

async function listSteps(project: Project): Promise<void> {
  const bindings = await loadBindings(project.sources);
  for (const module of bindings.modules.filter((m) => m.stepDefs.length > 0)) {
    console.log(`\n${module.specifier.startsWith('@') ? module.specifier : rel(project, module.file)}`);
    for (const def of module.stepDefs)
      console.log(`  ${def.keyword.padEnd(5)} ${def.expression}${def.meta.tier ? `  [${def.meta.tier}]` : ''}${def.meta.enters ? `  → ${def.meta.enters}` : ''}`);
  }
  const local = allContexts().reduce((n, c) => n + c.elements().length, 0);
  console.log(`\n  ${steps.length} step(s), ${allElements().length} element(s) + ${local} on ${allContexts().length} context(s), ` +
    `${allKinds().length} kind(s), ${allDatasets().length} dataset(s)`);
}

/** The package's own Playwright when it has one (the spec files resolve `@playwright/test` from
 * there, and Playwright insists on a single copy), else the library's. */
function playwrightCli(project: Project): string {
  for (const from of [project.packageDir, LIB_ROOT]) {
    try {
      const pkg = createRequire(join(from, 'package.json')).resolve('@playwright/test/package.json');
      return join(dirname(pkg), 'cli.js');
    } catch {
      continue;
    }
  }
  throw new Error('@playwright/test is not installed — add it to the package devDependencies');
}

async function run(project: Project, args: string[]): Promise<number> {
  // a run of named specs or folders checks the features behind them: one being written beside them does not stop it
  const paths = args.filter((a) => !a.startsWith('-')).map((a) => resolve(project.root, a));
  const only = paths.length === 0 ? undefined : new Set(listFiles(project.featuresDir, '.feature')
    .filter((f) => paths.some((p) => outFileFor(project.root, f).startsWith(p))));
  const code = await compile(project, true, false, only);
  return code !== 0 ? code : launch(project, args);
}

/** Playwright over the generated specs. Bindings load once per process (the registry is global),
 * so a caller compiles or checks before this, never both. */
function launch(project: Project, args: string[], env: Record<string, string> = {}): Promise<number> {
  const config = ['playwright.config.js', 'playwright.config.ts'].map((f) => join(LIB_DIR, f)).find(existsSync);
  if (!config)
    throw new Error(`no Playwright config in ${LIB_DIR} — run \`npm run build\` in the library`);
  const report = withJsonReport(project, args);
  const child = spawn(process.execPath, [playwrightCli(project), 'test', '--config', config, ...report.args],
    {cwd: project.root, stdio: 'inherit',
      env: {...process.env, ...report.env, ...env, BDD_ROOT: project.root, BDD_PROJECT: project.name}});
  return new Promise((resolve) => child.on('exit', (exit) => resolve(exit ?? 1)));
}

/** A `--reporter` on the command line replaces the config's reporters, the JSON report every run
 * leaves in test-results/report.json among them: json joins the reporters given, to the same file. */
function withJsonReport(project: Project, args: string[]): {args: string[]; env: Record<string, string>} {
  const i = args.findIndex((a) => a === '--reporter' || a.startsWith('--reporter='));
  const spaced = args[i] === '--reporter';
  const value = i < 0 ? '' : spaced ? args[i + 1] ?? '' : args[i].slice('--reporter='.length);
  if (i < 0 || value.split(',').includes('json'))
    return {args, env: {}};
  const named = process.env.PLAYWRIGHT_JSON_OUTPUT_FILE ?? process.env.PLAYWRIGHT_JSON_OUTPUT_NAME;
  return {args: [...args.slice(0, i), `--reporter=${value},json`, ...args.slice(i + (spaced ? 2 : 1))],
    env: named ? {} : {PLAYWRIGHT_JSON_OUTPUT_FILE: join(project.root, 'test-results', 'report.json')}};
}

function exec(command: string, args: string[], cwd: string): Promise<number> {
  const child = spawn(command, args, {cwd, stdio: 'inherit', shell: process.platform === 'win32'});
  return new Promise((resolve) => child.on('exit', (exit) => resolve(exit ?? 1)).on('error', () => resolve(127)));
}

/** The `help/` tree of the public repository, from the package upwards; `BDD_HELP_ROOT` names another. */
function helpRoot(from: string): string | undefined {
  if (process.env.BDD_HELP_ROOT)
    return process.env.BDD_HELP_ROOT;
  let dir = from;
  for (let i = 0; i < 6; i++) {
    const candidate = join(dir, 'help');
    if (existsSync(join(candidate, 'access')) && existsSync(join(candidate, 'visualize')))
      return candidate;
    const parent = dirname(dir);
    if (parent === dir)
      break;
    dir = parent;
  }
  return undefined;
}

/** `guide <feature|spec paths…> [--gif] [--demo] [--help-mode] [--help-pages] [--out <dir>] [--settle <ms>] [GIF options]
 * [playwright args]`: the scenarios run on one worker with `BDD_GUIDE` set, then every manifest the run wrote
 * goes through `tool/guide-render.py`. A `.feature` path names its generated spec. `--help-mode` films a GIF
 * for a help page (an 800x500 window unless BDD_GUIDE_VIEWPORT says otherwise, filmed at 2x; the bare shell,
 * no captions, no set-up); `--help-pages` is that for every feature tagged `@help:<page dir>`, and puts each
 * scenario's GIF and thumb into `<help root>/<page dir>/img/<scenario slug>.gif` — beside the one there, as
 * `<slug>.new.gif`, for the reviewer to swap in. `--demo`, or a feature's `@demo` tag, films a demo clip: a
 * 1080p video without captions or set-up, opened by a title card, of a 1600x900 window at twice its density (the
 * page a little larger than on a 1080p screen). A feature's `@viewport:<w>x<h>` films it in that window instead,
 * at twice its density below 1920 px. `--fps`, `--colors`,
 * `--size` and `--hold` set the GIF's rate, palette, width and the hold of its last frame. */
async function guide(project: Project, flags: string[]): Promise<number> {
  const args = [...flags];
  const take = (name: string): boolean => {
    const i = args.indexOf(name);
    if (i >= 0)
      args.splice(i, 1);
    return i >= 0;
  };
  const option = (name: string): string | undefined => {
    const i = args.indexOf(name);
    if (i < 0)
      return undefined;
    return args.splice(i, 2)[1];
  };
  const helpPages = take('--help-pages');
  const help = take('--help-mode') || helpPages;
  const gif = take('--gif') || help;
  const demoAll = take('--demo');
  const out = resolve(project.root, option('--out') ?? 'guides');
  const settle = option('--settle');
  const render = [['--fps', '--gif-fps'], ['--colors', '--colors'], ['--size', '--gif-width'], ['--hold', '--end-hold']]
    .flatMap(([name, as]) => {
      const value = option(name);
      return value === undefined ? [] : [as, value];
    });
  const tagged = helpPages ?
    listFiles(project.featuresDir, '.feature').map((file) => parseFeature(file, readFileSync(file, 'utf8')))
      .filter((f) => f.tags.some((t) => t.startsWith('@help:'))) :
    [];
  const named = args.filter((a) => a.endsWith('.feature')).map((a) => resolve(project.root, a));
  const models = [...named.map((f) => parseFeature(f, readFileSync(f, 'utf8'))), ...tagged];
  const rest = args.filter((a) => !a.endsWith('.feature'));
  if (models.length === 0 && !rest.some((a) => !a.startsWith('-'))) {
    console.error('guide: name the features to film, or --help-pages for every @help: tagged one');
    return 2;
  }
  // a guide is often a feature written a minute ago: its spec is generated here, not demanded, and only
  // its own errors stop the film
  const compiled = await compile(project, false, false, new Set(models.map((m) => m.path)));
  if (compiled !== 0)
    return compiled;
  // one run per way of filming: a demo and the window a feature asks for are the browser's settings
  const runs = new Map<string, {demo: boolean; view?: string; specs: string[]}>();
  for (const model of models) {
    const demo = demoAll || model.tags.includes('@demo');
    const view = model.tags.find((t) => t.startsWith('@viewport:'))?.slice('@viewport:'.length);
    const key = `${demo} ${view ?? ''}`;
    if (!runs.has(key))
      runs.set(key, {demo, view, specs: []});
    // Playwright reads its file arguments as regular expressions: a Windows path's backslashes
    // would escape the letters after them
    runs.get(key)!.specs.push(rel(project, outFileFor(project.root, model.path)));
  }
  if (runs.size === 0)
    runs.set('', {demo: demoAll, specs: []});
  const started = Date.now();
  let code = 0;
  for (const [i, run] of [...runs.values()].entries()) {
    const view = run.view ?? process.env.BDD_GUIDE_VIEWPORT ?? (help ? '800x500' : run.demo ? '1600x900' : undefined);
    const width = Number(/^(\d+)x\d+$/.exec(view ?? '')?.[1] ?? 1920);
    // each run keeps its results apart, so two films of one project can go on at once
    const results = join(project.root, 'test-results', `guide-${process.pid}-${i}`);
    const exit = await launch(project, ['--workers=1', ...run.specs, ...rest], {BDD_GUIDE: out,
      BDD_OUTPUT_DIR: results, PLAYWRIGHT_JSON_OUTPUT_FILE: join(results, 'report.json'),
      ...(settle ? {BDD_GUIDE_SETTLE: settle} : {}),
      ...(view ? {BDD_GUIDE_VIEWPORT: view} : {}),
      ...(help ? {BDD_GUIDE_HELP: '1'} : {}),
      ...(run.demo ? {BDD_GUIDE_DEMO: '1', BDD_GUIDE_SCALE: width < 1920 ? '2' : '1'} : {})});
    code ||= exit;
    if (exit === 0)
      rmSync(results, {recursive: true, force: true});
  }
  // the scenarios of the features filmed: another film writing into the same folder meanwhile is its own
  const filmed = new Set(models.flatMap((m) => m.scenarios.map((sc) => join(out, slugOf(m.name), slugOf(sc.name)))));
  const manifests = listFiles(out, '.json')
    .filter((f) => basename(f) === 'steps.json' && statSync(f).mtimeMs >= started && (filmed.size === 0 || filmed.has(dirname(f))));
  if (manifests.length === 0) {
    console.error(`guide: the run wrote no manifest under ${out}`);
    return code || 1;
  }
  const python = process.env.BDD_PYTHON ?? (process.platform === 'win32' ? 'py' : 'python3');
  const script = join(LIB_ROOT, 'tool', 'guide-render.py');
  let failed = 0;
  for (const manifest of manifests) {
    const dir = dirname(manifest);
    const filmed = JSON.parse(readFileSync(manifest, 'utf8'));
    // what the title card says under the title: the scenario's own description, else the feature's first paragraph
    const model = models.find((m) => slugOf(m.name) === slugOf(filmed.feature));
    const scenario = model?.scenarios.find((s) => slugOf(s.name) === slugOf(filmed.scenario));
    const description = (scenario?.description || model?.description.split(/\n\s*\n/)[0] || '').replace(/\s*\n\s*/g, ' ');
    if (description !== filmed.description)
      writeFileSync(manifest, JSON.stringify({...filmed, description}, null, 2));
    const demo = demoAll || (filmed.tags ?? []).includes('@demo');
    console.log(`guide: ${relative(project.root, dir)}`);
    const exit = await exec(python, [script, dir, ...(gif ? ['--gif'] : []), ...(help ? ['--help'] : []),
      ...(demo ? ['--demo'] : []), ...render], project.root);
    if (exit !== 0) {
      console.error(`guide: rendering failed (${exit}) — ${python} ${script}`);
      failed++;
    }
  }
  if (helpPages && failed === 0) {
    const helpDir = helpRoot(project.packageDir);
    if (!helpDir) {
      console.error('guide: no help/ tree found above the package (BDD_HELP_ROOT names one)');
      return 1;
    }
    for (const feature of tagged) {
      const page = feature.tags.find((t) => t.startsWith('@help:'))!.slice('@help:'.length);
      const img = join(helpDir, ...page.split('/'), 'img');
      mkdirSync(img, {recursive: true});
      for (const scenario of feature.scenarios) {
        const from = join(out, slugOf(feature.name), slugOf(scenario.name));
        if (!existsSync(join(from, 'guide.gif')))
          continue;
        // a GIF the page already has is replaced only once the reviewer approves the new one
        const name = existsSync(join(img, `${slugOf(scenario.name)}.gif`)) ? `${slugOf(scenario.name)}.new` : slugOf(scenario.name);
        copyFileSync(join(from, 'guide.gif'), join(img, `${name}.gif`));
        copyFileSync(join(from, 'guide-thumb.png'), join(img, `${name}-thumb.png`));
        console.log(`guide: ${relative(helpDir, img)}/${name}.gif${name.endsWith('.new') ? ' — beside the one the page shows, for review' : ''}`);
      }
    }
  }
  return code || (failed > 0 ? 1 : 0);
}

function libraryResolvesFrom(dir: string): boolean {
  try {
    createRequire(join(dir, 'package.json')).resolve(`${PACKAGE_NAME}/package.json`);
    return true;
  } catch {
    return false;
  }
}

/** Scaffolds the project, then compiles the smoke feature when the library is already installed in
 * the package (a fresh package compiles after `npm i`). */
async function init(cwd: string): Promise<number> {
  const result = scaffold(cwd);
  for (const file of result.created)
    console.log(`created ${file}`);
  for (const file of result.skipped)
    console.log(`exists  ${file}`);
  for (const note of result.notes)
    console.log(`note: ${note}`);
  const tiers = availableTiers();
  if (tiers.length > 0)
    console.log(`tiers a package can opt into in bdd/bdd.config.json: ${tiers.join(', ')}`);
  if (!libraryResolvesFrom(cwd)) {
    console.log(`next: npm i (installs ${PACKAGE_NAME} and @playwright/test), then \`grok-bdd compile\` and \`grok-bdd run\``);
    return 0;
  }
  return compile(loadProject(cwd), false);
}

/** Wires the package in `cwd` to this checkout of the library — `node_modules/@datagrok-libraries/bdd`
 * to the library, and `node_modules/@playwright/test` to the library's copy, because Playwright
 * refuses to be loaded twice in one run and the library's runtime resolves it from its own path.
 * For a library that is linked or not on npm; the package's own copies go to
 * `node_modules/.bdd-link-backup`, `--undo` (or an `npm ci`) puts them back. */
function link(cwd: string, undo: boolean): number {
  const nm = join(cwd, 'node_modules');
  // a project's bdd/ has a package.json of its own ({"type": "module"}) but no node_modules
  if (!existsSync(join(cwd, 'package.json')) || !existsSync(nm)) {
    console.error('grok-bdd link: run it in the package directory (the one with node_modules)');
    return 2;
  }
  const backup = join(nm, '.bdd-link-backup');
  const same = (a: string, b: string) => {
    try { return realpathSync(a) === realpathSync(b); }
    catch { return false; }
  };
  const kindOf = (path: string): 'link' | 'dir' | 'none' => {
    try { return lstatSync(path).isSymbolicLink() ? 'link' : 'dir'; }
    catch { return 'none'; }
  };
  const unlink = (path: string) => {
    try { rmdirSync(path); }
    catch { unlinkSync(path); }
  };
  const show = (path: string) => relative(cwd, path).split(sep).join('/');
  const links: [string, string, string][] = [
    [join(nm, '@datagrok-libraries', 'bdd'), LIB_ROOT, 'bdd'],
    [join(nm, '@playwright', 'test'), join(LIB_ROOT, 'node_modules', '@playwright', 'test'), 'test'],
  ];
  for (const [path, target, name] of links) {
    const kind = kindOf(path);
    if (undo) {
      if (kind === 'link')
        unlink(path);
      if (existsSync(join(backup, name))) {
        renameSync(join(backup, name), path);
        console.log(`restored ${show(path)}`);
      }
      else if (kind === 'link')
        console.log(`removed ${show(path)} (nothing to restore: npm i puts the package's own copy back)`);
      continue;
    }
    if (!existsSync(target)) {
      console.log(`skipped ${show(path)}: ${target} is missing — npm ci && npm run build in the library first`);
      continue;
    }
    if (kind === 'link' && same(path, target)) {
      console.log(`already linked ${show(path)} → ${target}`);
      continue;
    }
    if (kind === 'dir') {
      mkdirSync(backup, {recursive: true});
      rmSync(join(backup, name), {recursive: true, force: true});
      renameSync(path, join(backup, name));
    }
    else if (kind === 'link')
      unlink(path);
    mkdirSync(dirname(path), {recursive: true});
    symlinkSync(target, path, 'junction');
    console.log(`linked ${show(path)} → ${target}`);
  }
  if (!undo)
    console.log(`next: npx playwright install chromium in ${LIB_ROOT} (once per Playwright version), then grok-bdd run`);
  return 0;
}

/** The library's sources newer than their build. A package's compile and run load the library from
 * `dist/`, so after a pull without a rebuild a sound feature fails on a kind or a step only the
 * sources have. */
function staleBuild(): string[] {
  // a published copy ships no src/, and its install gives the files no meaningful order in time
  if (LIB_DIR === LIB_ROOT || !existsSync(join(LIB_ROOT, 'src')))
    return [];
  const sources = [...listFiles(join(LIB_ROOT, 'src'), '.ts'), ...listFiles(join(LIB_ROOT, 'bindings'), '.ts'),
    join(LIB_ROOT, 'playwright.config.ts')];
  return sources.filter((file) => {
    const built = join(LIB_DIR, relative(LIB_ROOT, file)).replace(/\.ts$/, '.js');
    return !existsSync(built) || statSync(built).mtimeMs < statSync(file).mtimeMs;
  });
}

async function main(): Promise<number> {
  const [command = 'compile', ...flags] = process.argv.slice(2);
  if (command === '--help' || command === '-h' || command === 'help') {
    console.log('grok-bdd init | link [--undo] | compile [--check] [--verbose] | lint | list-steps | ' +
      'run [playwright args] | guide <features> [--gif] [--demo] [--help-mode] [--help-pages] [--out <dir>] [--settle <ms>]');
    return 0;
  }
  if (command === 'init')
    return init(process.cwd());
  if (command === 'link')
    return link(process.cwd(), flags.includes('--undo'));
  const stale = staleBuild();
  if (stale.length > 0) {
    const first = relative(LIB_ROOT, stale[0]).split(sep).join('/');
    console.error(`grok-bdd: the library's build is older than its sources (${first}` +
      `${stale.length > 1 ? ` and ${stale.length - 1} more` : ''}) — run \`npm run build\` in ${LIB_ROOT}`);
    return 2;
  }
  const project = loadProject(process.cwd());
  console.log(`bdd: ${project.name} at ${project.root}${project.tiers.length > 0 ? ` (tiers: ${project.tiers.join(', ')})` : ''}`);
  switch (command) {
    case 'compile':
      return compile(project, flags.includes('--check'), flags.includes('--verbose'));
    case 'check':
      return compile(project, true, flags.includes('--verbose'));
    case 'lint':
      return report((await compileAll(project)).diagnostics, true) > 0 ? 1 : 0;
    case 'list-steps':
      await listSteps(project);
      return 0;
    case 'run':
      return run(project, flags);
    case 'guide':
      return guide(project, flags);
    default:
      console.error(`unknown command "${command}" — init | link [--undo] | compile [--check] [--verbose] | lint | ` +
        'list-steps | run [playwright args] | guide <features> [--gif] [--demo] [--help-mode] [--help-pages]');
      return 2;
  }
}

process.exitCode = await main();
