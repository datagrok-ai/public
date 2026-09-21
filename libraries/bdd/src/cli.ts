/* grok-bdd <command> — run inside a package (features under bdd/features) or inside the library.
     init                  bootstrap bdd/ in the package: config, starter bindings, a smoke feature, manifest entries
     compile [--check] [--verbose]  features/** → generated/**.test.ts; --check reports drift instead of writing, --verbose prints how every phrase resolves
     lint                  diagnostics only, notes included
     list-steps            the vocabulary this project sees (base tiers + its tiers + its own bindings)
     run [playwright args] compile --check, then Playwright over generated/ with the library's config
     guide <features> [--gif] [--help-pages] [--out <dir>] [--settle <ms>] [playwright args]
                           run the scenarios in guide mode, render each into a how-to video (tool/guide-render.py)
     link [--undo]         wire the package to this checkout of the library (one Playwright per run) */
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

async function compile(project: Project, check: boolean, notes = false): Promise<number> {
  const {compiled, diagnostics} = await compileAll(project);
  const errors = report(diagnostics, notes);
  const expected = new Set(compiled.map((c) => c.outFile));
  const orphans = listFiles(project.generatedDir, '.test.ts').filter((f) => !expected.has(f));
  let drift = 0;
  for (const c of compiled) {
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
  const code = await compile(project, true);
  return code !== 0 ? code : launch(project, args);
}

/** Playwright over the generated specs. Bindings load once per process (the registry is global),
 * so a caller compiles or checks before this, never both. */
function launch(project: Project, args: string[], env: Record<string, string> = {}): Promise<number> {
  const config = ['playwright.config.js', 'playwright.config.ts'].map((f) => join(LIB_DIR, f)).find(existsSync);
  if (!config)
    throw new Error(`no Playwright config in ${LIB_DIR} — run \`npm run build\` in the library`);
  const child = spawn(process.execPath, [playwrightCli(project), 'test', '--config', config, ...args],
    {cwd: project.root, stdio: 'inherit', env: {...process.env, ...env, BDD_ROOT: project.root}});
  return new Promise((resolve) => child.on('exit', (exit) => resolve(exit ?? 1)));
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

/** `guide <feature|spec paths…> [--gif] [--help-pages] [--out <dir>] [--settle <ms>] [playwright args]`:
 * the scenarios run on one worker with `BDD_GUIDE` set, then every manifest the run wrote goes
 * through `tool/guide-render.py`. A `.feature` path names its generated spec. `--help-pages` adds
 * every feature tagged `@help:<page dir>` and copies each scenario's GIF and thumb into
 * `<help root>/<page dir>/img/<scenario slug>.gif` — the docs' own layout. */
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
  const gif = take('--gif') || helpPages;
  const out = resolve(project.root, option('--out') ?? 'guides');
  const settle = option('--settle');
  const tagged = helpPages ?
    listFiles(project.featuresDir, '.feature').map((file) => parseFeature(file, readFileSync(file, 'utf8')))
      .filter((f) => f.tags.some((t) => t.startsWith('@help:'))) :
    [];
  const named = args.filter((a) => a.endsWith('.feature')).map((a) => resolve(project.root, a));
  const features = [...named, ...tagged.map((f) => f.path)];
  // Playwright reads its file arguments as regular expressions: a Windows path's backslashes
  // would escape the letters after them
  const specs = features.map((f) => rel(project, outFileFor(project.root, f)));
  const rest = args.filter((a) => !a.endsWith('.feature'));
  if (specs.length === 0 && !rest.some((a) => !a.startsWith('-'))) {
    console.error('guide: name the features to film, or --help-pages for every @help: tagged one');
    return 2;
  }
  // a guide is often a feature written a minute ago: its spec is generated here, not demanded
  const compiled = await compile(project, false);
  if (compiled !== 0)
    return compiled;
  const started = Date.now();
  const code = await launch(project, ['--workers=1', ...specs, ...rest],
    {BDD_GUIDE: out, ...(settle ? {BDD_GUIDE_SETTLE: settle} : {})});
  const manifests = listFiles(out, '.json')
    .filter((f) => basename(f) === 'steps.json' && statSync(f).mtimeMs >= started);
  if (manifests.length === 0) {
    console.error(`guide: the run wrote no manifest under ${out}`);
    return code || 1;
  }
  const python = process.env.BDD_PYTHON ?? (process.platform === 'win32' ? 'py' : 'python3');
  const script = join(LIB_ROOT, 'tool', 'guide-render.py');
  let failed = 0;
  for (const manifest of manifests) {
    const dir = dirname(manifest);
    console.log(`guide: ${relative(project.root, dir)}`);
    const exit = await exec(python, [script, dir, ...(gif ? ['--gif'] : [])], project.root);
    if (exit !== 0) {
      console.error(`guide: rendering failed (${exit}) — ${python} ${script}`);
      failed++;
    }
  }
  if (helpPages && failed === 0) {
    const help = helpRoot(project.packageDir);
    if (!help) {
      console.error('guide: no help/ tree found above the package (BDD_HELP_ROOT names one)');
      return 1;
    }
    for (const feature of tagged) {
      const page = feature.tags.find((t) => t.startsWith('@help:'))!.slice('@help:'.length);
      const img = join(help, ...page.split('/'), 'img');
      mkdirSync(img, {recursive: true});
      for (const scenario of feature.scenarios) {
        const from = join(out, slugOf(feature.name), slugOf(scenario.name));
        if (!existsSync(join(from, 'guide.gif')))
          continue;
        copyFileSync(join(from, 'guide.gif'), join(img, `${slugOf(scenario.name)}.gif`));
        copyFileSync(join(from, 'guide-thumb.png'), join(img, `${slugOf(scenario.name)}-thumb.png`));
        console.log(`guide: ${relative(help, img)}/${slugOf(scenario.name)}.gif`);
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

async function main(): Promise<number> {
  const [command = 'compile', ...flags] = process.argv.slice(2);
  if (command === '--help' || command === '-h' || command === 'help') {
    console.log('grok-bdd init | link [--undo] | compile [--check] [--verbose] | lint | list-steps | ' +
      'run [playwright args] | guide <features> [--gif] [--help-pages] [--out <dir>] [--settle <ms>]');
    return 0;
  }
  if (command === 'init')
    return init(process.cwd());
  if (command === 'link')
    return link(process.cwd(), flags.includes('--undo'));
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
        'list-steps | run [playwright args] | guide <features> [--gif] [--help-pages]');
      return 2;
  }
}

process.exitCode = await main();
