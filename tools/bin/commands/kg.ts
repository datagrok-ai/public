/// `grok kg check|gen|build|query|impact|tests-for|explain|find` — the knowledge-graph type files,
/// home documents, the built graph and the index over it (conventions.md §11.1, build-plan.md).
import * as fs from 'fs';
import * as path from 'path';
import {loadTypeSystem, Issue, TypeSystem} from '../utils/kg/types';
import {loadHomes, makeReport, CheckReport} from '../utils/kg/homes';
import {generate, writeOutputs} from '../utils/kg/gen';
import {Emitter} from '../utils/kg/build/emitter';
import {selectExtractors, runExtractors, EXTRACTORS, Mode} from '../utils/kg/build/registry';
import {writeBuild, projectPublic, gitRevisions, batchId, toolsVersion, Manifest} from '../utils/kg/build/write';
import {loadKuzu, load as loadIndex, open, run, MISSING_KUZU, LoadResult, TableRows} from '../utils/kg/kuzu';
import {impact, testsFor, explain, find, printOps, resolveTarget, coverageNote, OpsResult, DEFAULT_LIMIT} from '../utils/kg/ops';
import {OutputFormat, printOutput} from '../utils/server-output';
import {HELP_KG} from './help';

const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const PLANNED_VERBS = ['report'];
const OPS = ['impact', 'tests-for', 'explain', 'find'];

export async function kg(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1).map(String);
  const verb: string | undefined = args[0];
  if (!verb || verb === 'help' || argv.help) {
    console.log(HELP_KG);
    return true;
  }
  if (PLANNED_VERBS.includes(verb)) return fail(`grok kg ${verb} is not implemented yet (conventions.md §11.1); check, gen, build, query and the bounded operations are`);
  const output: string = argv.output ?? argv.o ?? 'table';
  if (verb === 'query' || OPS.includes(verb)) {
    if (output !== 'table' && output !== 'json' && output !== 'csv') return fail(`--output must be table, json or csv, got '${output}'`);
    return graph(verb, args.slice(1), argv, output as OutputFormat);
  }
  if (verb !== 'check' && verb !== 'gen' && verb !== 'build') {
    console.error(`unknown verb '${verb}'`);
    return false;
  }
  if (args.length > 1) return fail(`unexpected argument '${args[1]}': grok kg ${verb} takes options only`);
  if (output !== 'table' && output !== 'json') return fail(`--output must be table or json, got '${output}'`);
  const quiet = argv.quiet === true;
  const typesOnly = argv['types-only'] === true;
  if (verb === 'gen' && typesOnly) return fail('--types-only cannot be combined with gen: feature-tree.md is generated from the home documents');

  const kgRoot = argv.kg ? path.resolve(String(argv.kg)) : findKgRoot(process.cwd());
  if (!kgRoot || !fs.existsSync(path.join(kgRoot, 'schema.yaml')))
    return fail(kgRoot ? `${slashes(kgRoot)}: no schema.yaml` : 'grok kg currently needs the monorepo: run it inside a checkout with core/docs/knowledge-graph/schema.yaml and public/, or pass --kg <dir>');
  const repoRoot = path.resolve(kgRoot, '..', '..', '..');
  if (verb === 'build') return build(argv, kgRoot, repoRoot, output);

  const system = loadTypeSystem(kgRoot);
  const homes = typesOnly ? null : loadHomes(system, repoRoot);
  const report = makeReport(system, homes);
  if (verb === 'gen') {
    const generated = generate(system, kgRoot, repoRoot, homes);
    report.errors.push(...generated.errors);
    if (report.errors.length)
      report.warnings.push({file: kgRoot, code: 'gen-skipped', message: 'nothing generated while check reports errors'});
    else {
      const result = writeOutputs(generated.outputs, argv.check === true);
      report.written = result.written;
      report.stale = result.stale;
      for (const file of result.stale)
        report.errors.push({file, code: 'stale', message: 'stale: run grok kg gen'});
    }
  }
  for (const issue of [...report.errors, ...report.warnings]) issue.file = rel(issue.file, repoRoot);
  report.written = report.written?.map((f) => rel(f, repoRoot));
  report.stale = report.stale?.map((f) => rel(f, repoRoot));
  print(report, output, quiet);
  if (report.errors.length) process.exitCode = 1;
  return true;
}

/** `grok kg build`: the extractors into JSONL and a manifest under `.kg/` (public/.kg/ with --public). */
async function build(argv: any, kgRoot: string, repoRoot: string, output: string): Promise<boolean> {
  const mode: Mode = argv.public === true ? 'public' : 'full';
  const only = argv.only === undefined ? undefined : String(argv.only).split(',').map((s) => s.trim()).filter(Boolean);
  const {selected, unknown} = selectExtractors(mode, only);
  if (unknown.length) return fail(`--only names unknown extractors: ${unknown.join(', ')} (known: ${EXTRACTORS.map((e) => e.name).join(', ')})`);
  const system = loadTypeSystem(kgRoot);
  if (system.errors.length) {
    for (const e of system.errors) console.error(`${rel(e.file, repoRoot)}: ${e.message}`);
    return fail(`${system.errors.length} type-system error${system.errors.length === 1 ? '' : 's'}; nothing built (run grok kg check)`);
  }
  const builder = toolsVersion();
  const revisions = gitRevisions(repoRoot);
  const batch = batchId(revisions, system.schemaVersion, builder);
  const emitter = new Emitter(system, batch);
  const backlogDir = argv.backlog === undefined ? undefined : path.resolve(String(argv.backlog));
  await runExtractors(selected, {system, kgRoot, repoRoot, mode, backlogDir}, emitter);
  let graph = emitter.finalize();
  if (mode === 'public') graph = projectPublic(graph, system);
  const outRoot = argv.out === undefined ? path.join(repoRoot, ...(mode === 'public' ? ['public', '.kg'] : ['.kg'])) : path.resolve(String(argv.out));
  const manifest = writeBuild(graph, outRoot, {mode, batch, builder, schemaVersion: system.schemaVersion, revisions});
  if (output === 'json') console.log(JSON.stringify(manifest, null, 2));
  else console.log(summary(manifest, rel(outRoot, repoRoot)));
  return index(argv, system, outRoot, repoRoot);
}

/** The JSONL is canonical: without the binding the build says so in one line and still succeeds. */
async function index(argv: any, system: TypeSystem, outRoot: string, repoRoot: string): Promise<boolean> {
  if (argv.db === false) return true;
  if (!loadKuzu()) {
    console.log(`index not built: ${MISSING_KUZU}`);
    return true;
  }
  try {
    const loaded = await loadIndex(outRoot, system);
    console.log(indexSummary(loaded, rel(loaded.db, repoRoot)));
  }
  catch (e: any) {
    return fail(`index not built: ${e.message}`);
  }
  return true;
}

/** `grok kg query` and the bounded operations: read-only, over `<repoRoot>/.kg/kg.kuzu`. */
async function graph(verb: string, args: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const dir = graphDir(argv);
  if (!dir) return fail('no graph found: run grok kg build inside the monorepo, or pass --kg <folder>');
  const cypher = verb !== 'query' ? '' : argv.file ? fs.readFileSync(path.resolve(String(argv.file)), 'utf8') : args.join(' ');
  const text = args.join(' ').trim();
  if (verb === 'query' && !cypher.trim()) return fail('grok kg query needs a Cypher statement, or --file <path>');
  if (verb !== 'query' && !text) return fail(`grok kg ${verb} needs ${verb === 'find' ? 'a text to search for' : 'a path or a ~id'}`);
  if (!fs.existsSync(path.join(dir, 'kg.kuzu'))) return fail(`${slashes(path.join(dir, 'kg.kuzu'))}: no index yet; run grok kg build`);
  const opened = await open(dir, true);
  if (!opened) {
    console.error(MISSING_KUZU);
    process.exitCode = 2;
    return true;
  }
  try {
    if (verb === 'query') {
      const {rows} = await run(opened.conn, cypher);
      printOutput(rows, output);
      return true;
    }
    const limit = Number(argv.limit) > 0 ? Number(argv.limit) : DEFAULT_LIMIT;
    const sources = readManifest(dir)?.sources;
    let result: OpsResult;
    if (verb === 'find') result = await find(opened.conn, text, {limit});
    else {
      const target = await resolveTarget(opened.conn, text);
      if (!target) return fail(`${text}: no such node; try grok kg find ${text.replace(/^~/, '')}`);
      result = verb === 'impact' ? await impact(opened.conn, target, {limit})
        : verb === 'tests-for' ? await testsFor(opened.conn, target, {limit})
          : await explain(opened.conn, target, {limit});
    }
    result.note = coverageNote(sources);
    printOps(result, output);
    return true;
  }
  catch (e: any) {
    return fail(e.message ?? String(e));
  }
  finally {
    await opened.conn.close();
    await opened.db.close();
  }
}

/** Full mode: `.kg` beside the type files. Public mode: the nearest committed `.kg/manifest.json`. */
function graphDir(argv: any): string | null {
  if (argv.out !== undefined) return path.resolve(String(argv.out));
  const kgRoot = argv.kg ? path.resolve(String(argv.kg)) : findKgRoot(process.cwd());
  if (kgRoot) return path.join(path.resolve(kgRoot, '..', '..', '..'), '.kg');
  let dir = path.resolve(process.cwd());
  for (;;) {
    if (fs.existsSync(path.join(dir, '.kg', 'manifest.json'))) return path.join(dir, '.kg');
    const parent = path.dirname(dir);
    if (parent === dir) return null;
    dir = parent;
  }
}

function readManifest(dir: string): Manifest | undefined {
  const file = path.join(dir, 'manifest.json');
  return fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) as Manifest : undefined;
}

function indexSummary(r: LoadResult, db: string): string {
  const total = (tables: TableRows[]) => tables.reduce((sum, t) => sum + t.rows, 0);
  const list = (tables: TableRows[]) => tables.filter((t) => t.rows).map((t) => `${t.table} ${t.rows}`).join(', ');
  return `loaded ${db} in ${(r.ms / 1000).toFixed(1)}s: ${total(r.nodes)} nodes (${list(r.nodes)}), ${total(r.rels)} edges (${list(r.rels)}); ` +
    `${(r.bytes / 1048576).toFixed(1)} MB on disk${r.parameterized ? `; ${r.parameterized} row${r.parameterized === 1 ? '' : 's'} inserted one by one` : ''}`;
}

function summary(m: Manifest, out: string): string {
  const total = (counts: Record<string, number>) => Object.values(counts).reduce((a, b) => a + b, 0);
  const list = (counts: Record<string, number | string>) => Object.entries(counts).map(([k, n]) => `${k} ${n}`).join(', ');
  const problems = Object.entries(m.problems).filter(([, n]) => n > 0).map(([k, n]) => `${k} ${n}`).join(', ') || 'none';
  return `wrote ${out}: ${total(m.counts.nodes)} nodes (${list(m.counts.nodes)}), ${total(m.counts.edges)} edges (${list(m.counts.edges)}); ` +
    `sources: ${list(m.sources) || 'none'}; problems: ${problems}; batch ${m.batch} (${m.mode})`;
}

function fail(message: string): boolean {
  console.error(message);
  process.exitCode = 1;
  return true;
}

function findKgRoot(from: string): string | null {
  let dir = path.resolve(from);
  for (;;) {
    if (fs.existsSync(path.join(dir, KG_DIR, 'schema.yaml')) && fs.existsSync(path.join(dir, 'public')))
      return path.join(dir, KG_DIR);
    const parent = path.dirname(dir);
    if (parent === dir) return null;
    dir = parent;
  }
}

function print(report: CheckReport, output: string, quiet: boolean): void {
  if (output === 'json') {
    console.log(JSON.stringify(report, null, 2));
    return;
  }
  const where = (issue: Issue) => `${issue.file}${issue.line ? `:${issue.line}` : ''}`;
  for (const e of report.errors) console.log(`${where(e)}: ${e.message}`);
  if (quiet) return;
  for (const w of report.warnings) console.log(`warning: ${where(w)}: ${w.message}`);
  for (const f of report.written ?? []) console.log(`wrote ${f}`);
  const homes = Object.entries(report.homes).map(([t, n]) => `${t} ${n}`).join(', ');
  const total = Object.values(report.homes).reduce((a, b) => a + b, 0);
  console.log(`${report.types.nodes} node types, ${report.types.edges} edge types, ${report.types.prefixes} prefixes; ` +
    `${report.scanned} files scanned, ${total} home${total === 1 ? '' : 's'}${homes ? ` (${homes})` : ''}, ` +
    `${report.annotatedPages} annotated page${report.annotatedPages === 1 ? '' : 's'}; ` +
    `${report.citations.doc} doc links and ${report.citations.code} code citations checked; ` +
    `${report.unresolvedExternal.length} unresolved external reference${report.unresolvedExternal.length === 1 ? '' : 's'}; ` +
    `${report.stubs.length} stub${report.stubs.length === 1 ? '' : 's'} needed; ` +
    `${report.errors.length} error${report.errors.length === 1 ? '' : 's'}, ${report.warnings.length} warning${report.warnings.length === 1 ? '' : 's'}`);
}

function rel(file: string, repoRoot: string): string {
  const r = path.isAbsolute(file) ? path.relative(repoRoot, file) : file;
  return slashes(r.startsWith('..') || path.isAbsolute(r) ? file : r);
}

function slashes(p: string): string {
  return p.replace(/\\/g, '/');
}
