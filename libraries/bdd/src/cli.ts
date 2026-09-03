/* grok-bdd <command> — run inside a package (features under bdd/features) or inside the library.
     init                  bootstrap bdd/ in the package: config, starter bindings, a smoke feature, manifest entries
     compile [--check]     features/** → generated/**.test.ts; --check reports drift instead of writing
     lint                  diagnostics only
     list-steps            the vocabulary this project sees (base tiers + its tiers + its own bindings)
     run [playwright args] compile --check, then Playwright over generated/ with the library's config */
import {spawn} from 'node:child_process';
import {existsSync, mkdirSync, readFileSync, rmSync, writeFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {dirname, join, relative, sep} from 'node:path';
import {compileFeature, CompiledFeature, Diagnostic} from './compile.js';
import {listFiles, loadBindings} from './discover.js';
import {GherkinParseError, parseFeature} from './gherkin.js';
import {scaffold} from './init.js';
import {StepMatcher} from './match.js';
import {availableTiers, LIB_DIR, LIB_ROOT, loadProject, PACKAGE_NAME, Project} from './project.js';
import {allContexts, allDatasets, allElements, allKinds, steps} from './registry.js';

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

function report(diagnostics: Diagnostic[]): number {
  for (const d of diagnostics)
    console.log(`${d.file}:${d.line} ${d.level}: ${d.message}`);
  const errors = diagnostics.filter((d) => d.level === 'error').length;
  console.log(`${errors} error(s), ${diagnostics.length - errors} note(s)`);
  return errors;
}

const eol = (s: string) => s.replace(/\r\n/g, '\n');

async function compile(project: Project, check: boolean): Promise<number> {
  const {compiled, diagnostics} = await compileAll(project);
  const errors = report(diagnostics);
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
  if (code !== 0)
    return code;
  const config = ['playwright.config.js', 'playwright.config.ts'].map((f) => join(LIB_DIR, f)).find(existsSync);
  if (!config)
    throw new Error(`no Playwright config in ${LIB_DIR} — run \`npm run build\` in the library`);
  const child = spawn(process.execPath, [playwrightCli(project), 'test', '--config', config, ...args],
    {cwd: project.root, stdio: 'inherit', env: {...process.env, BDD_ROOT: project.root}});
  return new Promise((resolve) => child.on('exit', (exit) => resolve(exit ?? 1)));
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

async function main(): Promise<number> {
  const [command = 'compile', ...flags] = process.argv.slice(2);
  if (command === '--help' || command === '-h' || command === 'help') {
    console.log('grok-bdd init | compile [--check] | lint | list-steps | run [playwright args]');
    return 0;
  }
  if (command === 'init')
    return init(process.cwd());
  const project = loadProject(process.cwd());
  console.log(`bdd: ${project.name} at ${project.root}${project.tiers.length > 0 ? ` (tiers: ${project.tiers.join(', ')})` : ''}`);
  switch (command) {
    case 'compile':
      return compile(project, flags.includes('--check'));
    case 'check':
      return compile(project, true);
    case 'lint':
      return report((await compileAll(project)).diagnostics) > 0 ? 1 : 0;
    case 'list-steps':
      await listSteps(project);
      return 0;
    case 'run':
      return run(project, flags);
    default:
      console.error(`unknown command "${command}" — init | compile [--check] | lint | list-steps | run [playwright args]`);
      return 2;
  }
}

process.exitCode = await main();
