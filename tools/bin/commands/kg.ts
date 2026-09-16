/// `grok kg check|gen|build|report|query|impact|tests-for|explain|find` — the knowledge-graph type
/// files, home documents, the built graph, the maintainer reports and the index over it
/// (conventions.md §11.1, build-plan.md).
import * as fs from 'fs';
import * as path from 'path';
import {loadTypeSystem, edgeGroups, Issue, TypeSystem} from '../utils/kg/types';
import {loadHomes, makeReport, CheckReport} from '../utils/kg/homes';
import {generate, writeOutputs} from '../utils/kg/gen';
import {Emitter} from '../utils/kg/build/emitter';
import {selectExtractors, runExtractors, provides, EXTRACTORS} from '../utils/kg/build/registry';
import {Mode} from '../utils/kg/build/context';
import {writeBuild, writeManifest, projectPublic, gitRevisions, buildInputs, batchId, toolsVersion} from '../utils/kg/build/write';
import {readManifest, generationDir, newGeneration, currentDir, readCurrent, publish, generations, gc, Manifest} from '../utils/kg/generation';
import {loadKuzu, load as loadIndex, open, run, memoryMb, MISSING_KUZU, BUILD_MEMORY_MB, LoadResult, TableRows} from '../utils/kg/kuzu';
import {impact, testsFor, explain, find, resolveTarget, DEFAULT_LIMIT} from '../utils/kg/ops';
import {Answer, caveatNotes} from '../utils/kg/answer';
import {printAnswer, printReport} from '../utils/kg/print';
import {readGraph, fromGraph, makeReport as buildReport, writeReports, REPORT_NAMES, ReportName, ReportFormat} from '../utils/kg/report';
import {OutputFormat, printOutput} from '../utils/server-output';
import {exportVis, hasVis} from '../utils/kg/vis';
import {serve as listen} from '../utils/kg/serve';
import {loadQuestions, ask, Question} from '../utils/kg/questions';
import {openBrowser} from '../utils/utils';
import {HELP_KG} from './help';

const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const OPS = ['impact', 'tests-for', 'explain', 'find'];
const VERBS = ['check', 'gen', 'build', 'report', 'gc', 'serve', 'ask'];
const SERVE_PORT = 7475;
/** How many generations `grok kg gc` keeps beside the current one. */
const KEEP_GENERATIONS = 2;
const BACKLOG_FALLBACK = 'C:/dg/backlog';

export async function kg(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1).map(String);
  const verb: string | undefined = args[0];
  if (!verb || verb === 'help' || argv.help) {
    console.log(HELP_KG);
    return true;
  }
  const output: string = argv.output ?? argv.o ?? 'table';
  if (verb === 'query' || OPS.includes(verb)) {
    if (output !== 'table' && output !== 'json' && output !== 'csv') return fail(`--output must be table, json or csv, got '${output}'`);
    return graph(verb, args.slice(1), argv, output as OutputFormat);
  }
  if (!VERBS.includes(verb)) {
    console.error(`unknown verb '${verb}'`);
    return false;
  }
  const takes = verb === 'report' || verb === 'ask' ? 2 : 1;
  if (args.length > takes) return fail(`unexpected argument '${args[takes]}': grok kg ${verb} takes ${takes === 2 ? `one ${verb === 'report' ? 'report' : 'question'} name and ` : ''}options only`);
  const formats = verb === 'report' ? ['table', 'json', 'md'] : verb === 'ask' ? ['table', 'json', 'csv'] : ['table', 'json'];
  if (!formats.includes(output)) return fail(`--output must be ${formats.slice(0, -1).join(', ')} or ${formats[formats.length - 1]}, got '${output}'`);
  if (verb === 'gc') return collect(argv, output);
  const quiet = argv.quiet === true;
  const typesOnly = argv['types-only'] === true;
  if (verb === 'gen' && typesOnly) return fail('--types-only cannot be combined with gen: feature-tree.md is generated from the home documents');

  const kgRoot = argv.kg ? path.resolve(String(argv.kg)) : findKgRoot(process.cwd());
  if (!kgRoot || !fs.existsSync(path.join(kgRoot, 'schema.yaml')))
    return fail(kgRoot ? `${slashes(kgRoot)}: no schema.yaml` : 'grok kg currently needs the monorepo: run it inside a checkout with core/docs/knowledge-graph/schema.yaml and public/, or pass --kg <dir>');
  const repoRoot = path.resolve(kgRoot, '..', '..', '..');
  if (verb === 'build') return build(argv, kgRoot, repoRoot, output);
  if (verb === 'report') return reportVerb(argv, kgRoot, repoRoot, args[1], output as ReportFormat);
  if (verb === 'serve') return serve(argv, kgRoot, repoRoot);
  if (verb === 'ask') return askVerb(argv, kgRoot, repoRoot, args[1], output as OutputFormat);

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

/** `grok kg build`: the extractors into JSONL and a manifest in one immutable generation under `.kg/gen/<batch>-<suffix>/`
 * (public/.kg/ with --public), which `<out>/current` names once everything, the index included, is complete. */
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
  const backlogDir = backlogRoot(repoRoot, argv.backlog);
  const root = argv.out === undefined ? path.join(repoRoot, ...(mode === 'public' ? ['public', '.kg'] : ['.kg'])) : path.resolve(String(argv.out));
  const batch = batchId(buildInputs({repoRoot, mode, schemaVersion: system.schemaVersion, builder, revisions,
    extractors: selected.map((e) => e.name), backlogDir, outRoot: root}));
  const emitter = new Emitter(system, batch);
  await runExtractors(selected, {system, kgRoot, repoRoot, mode, backlogDir}, emitter);
  let graph = emitter.finalize();
  if (mode === 'public') graph = projectPublic(graph, system);
  Object.assign(graph.manifest, {edge_groups: edgeGroups(system), provides: provides(selected)});
  const genDir = newGeneration(root, batch);
  const manifest = writeBuild(graph, genDir, {mode, batch, builder, schemaVersion: system.schemaVersion, revisions});
  if (mode !== 'public') writeReports(genDir, fromGraph(graph, repoRoot, manifest.sources, manifest.revisions), {system, repoRoot});
  const failure = await index(argv, system, genDir, repoRoot, manifest);
  if (failure) return fail(`index not built: ${failure}; ${slashes(genDir)} stays unpublished`);
  writeManifest(genDir, manifest);
  if (output === 'json') console.log(JSON.stringify(manifest, null, 2));
  else console.log(summary(manifest, rel(genDir, repoRoot)));
  return promote(root, path.basename(genDir), manifest);
}

/** The pointer moves only to a generation at least as usable as the one it leaves: a `--no-db` rebuild does not
 * take the index away from `query` and the operations. */
function promote(root: string, name: string, manifest: Manifest): boolean {
  const previous = readCurrent(root);
  const prev = previous ? readManifest(generationDir(root, previous)) : undefined;
  const indexed = prev !== undefined && prev.indexed_batch === prev.batch;
  if (!manifest.indexed_batch && indexed) {
    console.log(`current stays at ${previous}: generation ${name} has no index (--no-db); run grok kg build to index it`);
    return true;
  }
  publish(root, name);
  return true;
}

/** `grok kg gc`: the older generations go, the current one and the newest `--keep` stay. */
function collect(argv: any, output: string): boolean {
  const root = findOutRoot(argv);
  if (!root) return fail('no graph found: run grok kg build inside the monorepo, or pass --out <folder>');
  const keep = Number(argv.keep) > 0 ? Math.round(Number(argv.keep)) : KEEP_GENERATIONS;
  const before = generations(root).length;
  const result = gc(root, keep);
  if (output === 'json') console.log(JSON.stringify(result, null, 2));
  else {
    console.log(`${slashes(root)}: ${before} generation${before === 1 ? '' : 's'}, kept ${result.kept.join(', ') || 'none'}` +
      `${result.removed.length ? `, removed ${result.removed.join(', ')}` : ', removed none'}`);
    for (const batch of result.locked) console.log(`${batch}: in use, left alone`);
  }
  return true;
}

/** `grok kg report <name>`: one maintainer report over the JSONL a build already wrote (build-plan.md WO-8). */
function reportVerb(argv: any, kgRoot: string, repoRoot: string, name: string | undefined, output: ReportFormat): boolean {
  if (!name) return fail(`grok kg report needs a report name: ${REPORT_NAMES.join(', ')}`);
  if (!REPORT_NAMES.includes(name as ReportName)) return fail(`unknown report '${name}': ${REPORT_NAMES.join(', ')}`);
  const base = argv.diff === undefined ? undefined : String(argv.diff);
  if ((name === 'diff') !== (base !== undefined))
    return fail(name === 'diff' ? 'grok kg report diff needs the revision to compare against: --diff <ref>' : `--diff is for grok kg report diff, not ${name}`);
  const root = argv.out === undefined ? path.join(repoRoot, '.kg') : path.resolve(String(argv.out));
  const outRoot = currentDir(root);
  if (!outRoot) return fail(`${slashes(root)}: nothing built yet; run grok kg build`);
  const system = loadTypeSystem(kgRoot);
  const data = readGraph(outRoot, repoRoot, system, name as ReportName);
  printReport(buildReport(name as ReportName, data, {system, repoRoot, base}), output);
  return true;
}

/** The JSONL is canonical: without the binding the build says so in one line and still succeeds. A load that
 * fails leaves the generation without a manifest, so `current` keeps naming the last complete one. */
async function index(argv: any, system: TypeSystem, genDir: string, repoRoot: string, manifest: Manifest): Promise<string | undefined> {
  if (argv.db === false) return undefined;
  if (!loadKuzu()) {
    console.log(`index not built: ${MISSING_KUZU}`);
    return undefined;
  }
  try {
    const loaded = await loadIndex(genDir, system, memoryMb(argv.memory, BUILD_MEMORY_MB));
    manifest.indexed_batch = manifest.batch;
    manifest.index_memory_mb = loaded.memoryMb;
    manifest.index_platform = loaded.platform;
    console.log(indexSummary(loaded, rel(loaded.db, repoRoot)));
    return undefined;
  }
  catch (e: any) {
    return e.message;
  }
}

/** `grok kg query` and the bounded operations: read-only, over the current generation's `kg.kuzu`. */
async function graph(verb: string, args: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const root = findOutRoot(argv);
  const dir = root ? currentDir(root) : null;
  if (!dir) return fail('no graph found: run grok kg build inside the monorepo, or pass --kg <folder>');
  const cypher = verb !== 'query' ? '' : argv.file ? fs.readFileSync(path.resolve(String(argv.file)), 'utf8') : args.join(' ');
  const text = args.join(' ').trim();
  if (verb === 'query' && !cypher.trim()) return fail('grok kg query needs a Cypher statement, or --file <path>');
  if (verb !== 'query' && !text) return fail(`grok kg ${verb} needs ${verb === 'find' ? 'a text to search for' : 'a path or a ~id'}`);
  if (!fs.existsSync(path.join(dir, 'kg.kuzu'))) return fail(`${slashes(path.join(dir, 'kg.kuzu'))}: no index yet; run grok kg build`);
  const mismatch = indexMismatch(dir);
  if (mismatch) return fail(mismatch);
  const opened = await open(dir, true, argv.memory);
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
    const manifest = readManifest(dir);
    const options = {limit, groups: manifest?.edge_groups};
    let result: Answer;
    if (verb === 'find') result = await find(opened.conn, text, options);
    else {
      const target = await resolveTarget(opened.conn, text);
      if (!target) return fail(`${text}: no such node; try grok kg find ${text.replace(/^~/, '')}`);
      result = verb === 'impact' ? await impact(opened.conn, target, options)
        : verb === 'tests-for' ? await testsFor(opened.conn, target, options)
          : await explain(opened.conn, target, options);
    }
    result.notes = caveatNotes(manifest);
    printAnswer(result, output);
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

/** `grok kg serve`: the browser over the current generation; the render tier is exported into `<gen>/vis/` on
 * first start and kept there (a derived cache, never part of the batch). Runs until Ctrl-C. */
async function serve(argv: any, kgRoot: string, repoRoot: string): Promise<boolean> {
  const root = argv.out === undefined ? path.join(repoRoot, '.kg') : path.resolve(String(argv.out));
  const dir = currentDir(root);
  if (!dir) return fail(`${slashes(root)}: nothing built yet; run grok kg build`);
  if (!fs.existsSync(path.join(dir, 'kg.kuzu'))) return fail(`${slashes(path.join(dir, 'kg.kuzu'))}: no index yet; run grok kg build`);
  const mismatch = indexMismatch(dir);
  if (mismatch) return fail(mismatch);
  const manifest = readManifest(dir)!;
  const system = loadTypeSystem(kgRoot);
  if (system.errors.length) return fail(`${system.errors.length} type-system error${system.errors.length === 1 ? '' : 's'}; run grok kg check`);
  const loaded = loadQuestions(kgRoot);
  for (const e of loaded.errors) console.error(`questions: ${e}`);
  if (!hasVis(dir)) {
    const started = Date.now();
    const summary = exportVis(dir, system, manifest.batch);
    console.log(`exported ${rel(summary.dir, repoRoot)}: ${summary.nodes} nodes, ${summary.edges} edges` +
      `${summary.dropped ? `, ${summary.dropped} dangling dropped` : ''}, ${(summary.bytes / 1048576).toFixed(1)} MB in ${((Date.now() - started) / 1000).toFixed(1)}s`);
  }
  const opened = await open(dir, true, argv.memory);
  if (!opened) {
    console.error(MISSING_KUZU);
    process.exitCode = 2;
    return true;
  }
  const port = argv.port === undefined ? SERVE_PORT : Number(argv.port);
  if (!Number.isInteger(port) || port < 0 || port > 65535) return fail(`--port must be a port number, got '${argv.port}'`);
  const served = await listen({genDir: dir, repoRoot, manifest, system, questions: loaded.questions, db: opened.db, conn: opened.conn, port});
  console.log(`grok kg serve: ${served.url} over ${rel(dir, repoRoot)} (batch ${manifest.batch}); Ctrl-C stops it`);
  if (argv.open === true) openBrowser(served.url);
  await new Promise<void>((resolve) => {
    for (const signal of ['SIGINT', 'SIGTERM'] as const) process.once(signal, () => resolve());
  });
  await served.close();
  await opened.conn.close();
  await opened.db.close();
  return true;
}

/** `grok kg ask [<question>] [--set k=v]`: one of the questions under core/docs/knowledge-graph/questions/ against the
 * current generation, or the list of them. */
async function askVerb(argv: any, kgRoot: string, repoRoot: string, id: string | undefined, output: OutputFormat): Promise<boolean> {
  const loaded = loadQuestions(kgRoot);
  for (const e of loaded.errors) console.error(`questions: ${e}`);
  if (!id) {
    const rows = loaded.questions.map((q) => ({id: q.id, question: q.question, params: Object.keys(q.params).join(', '), status: q.status}));
    printOutput(rows, output);
    return true;
  }
  const question = loaded.questions.find((q) => q.id === id);
  if (!question) return fail(`no question '${id}'; grok kg ask lists them`);
  const given: Record<string, string> = {};
  for (const pair of [].concat(argv.set ?? []).map(String)) {
    const at = pair.indexOf('=');
    if (at < 1) return fail(`--set takes name=value, got '${pair}'`);
    given[pair.slice(0, at)] = pair.slice(at + 1);
  }
  const root = argv.out === undefined ? path.join(repoRoot, '.kg') : path.resolve(String(argv.out));
  const dir = currentDir(root);
  if (!dir || !fs.existsSync(path.join(dir, 'kg.kuzu'))) return fail('no index yet; run grok kg build');
  const mismatch = indexMismatch(dir);
  if (mismatch) return fail(mismatch);
  const opened = await open(dir, true, argv.memory);
  if (!opened) {
    console.error(MISSING_KUZU);
    process.exitCode = 2;
    return true;
  }
  try {
    const answer = await ask(opened.conn, question, given, readManifest(dir));
    if (output === 'json') console.log(JSON.stringify({id: question.id, question: question.question, params: answer.params, notes: answer.notes, columns: answer.columns, rows: answer.rows, ms: answer.ms}, null, 2));
    else {
      for (const note of answer.notes) console.log(`note: ${note}`);
      printOutput(answer.rows, output);
      if (output === 'table') console.log(`${answer.rows.length} row${answer.rows.length === 1 ? '' : 's'} in ${answer.ms} ms` +
        `${Object.keys(answer.params).length ? ` (${Object.entries(answer.params).map(([k, v]) => `${k}=${v}`).join(', ')})` : ''}`);
    }
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

/** Full mode: `.kg` beside the type files. Public mode: the nearest committed `.kg` snapshot. */
/** The backlog snapshot every consumer reads: `--backlog`, else `<repo>/../backlog`, else the dev-box clone; none of
 * the three means no ticket layer and no taxonomy. Resolved once, so the batch identity and the extractors agree. */
function backlogRoot(repoRoot: string, flag: unknown): string | undefined {
  if (flag !== undefined) return path.resolve(String(flag));
  return [path.resolve(repoRoot, '..', 'backlog'), BACKLOG_FALLBACK].find((d) => fs.existsSync(path.join(d, 'index.jsonl')));
}

function findOutRoot(argv: any): string | null {
  if (argv.out !== undefined) return path.resolve(String(argv.out));
  const kgRoot = argv.kg ? path.resolve(String(argv.kg)) : findKgRoot(process.cwd());
  if (kgRoot) return path.join(path.resolve(kgRoot, '..', '..', '..'), '.kg');
  let dir = path.resolve(process.cwd());
  for (;;) {
    if (currentDir(path.join(dir, '.kg'))) return path.join(dir, '.kg');
    const parent = path.dirname(dir);
    if (parent === dir) return null;
    dir = parent;
  }
}

/** The index in a generation must have been loaded from that generation; anything else is a mixed graph. */
function indexMismatch(dir: string): string | undefined {
  const manifest = readManifest(dir);
  if (manifest?.indexed_batch === manifest?.batch) return undefined;
  return `${slashes(path.join(dir, 'kg.kuzu'))}: this index was loaded from batch ${manifest?.indexed_batch ?? 'unknown'}, ` +
    `and the data beside it is ${manifest?.batch ?? 'unknown'}; run grok kg build`;
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
  const observed = m.inventory?.observed_files;
  const problems = Object.entries(m.problems).filter(([, n]) => n > 0)
    .map(([k, n]) => k === 'orphans' && observed ? `orphans ${n} of ${observed} observed files` : `${k} ${n}`).join(', ') || 'none';
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
