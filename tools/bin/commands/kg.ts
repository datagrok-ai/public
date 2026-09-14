/// `grok kg check|gen|build` — the knowledge-graph type files, home documents and the built graph
/// (core/docs/knowledge-graph/conventions.md §11.1, build-plan.md).
import * as fs from 'fs';
import * as path from 'path';
import {loadTypeSystem, Issue} from '../utils/kg/types';
import {loadHomes, makeReport, CheckReport} from '../utils/kg/homes';
import {generate, writeOutputs} from '../utils/kg/gen';
import {Emitter} from '../utils/kg/build/emitter';
import {selectExtractors, runExtractors, EXTRACTORS, Mode} from '../utils/kg/build/registry';
import {writeBuild, projectPublic, gitRevisions, batchId, toolsVersion, Manifest} from '../utils/kg/build/write';
import {HELP_KG} from './help';

const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const PLANNED_VERBS = ['query', 'report'];

export async function kg(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1).map(String);
  const verb: string | undefined = args[0];
  if (!verb || verb === 'help' || argv.help) {
    console.log(HELP_KG);
    return true;
  }
  if (PLANNED_VERBS.includes(verb)) return fail(`grok kg ${verb} is not implemented yet (conventions.md §11.1); check, gen and build are`);
  if (verb !== 'check' && verb !== 'gen' && verb !== 'build') {
    console.error(`unknown verb '${verb}'`);
    return false;
  }
  if (args.length > 1) return fail(`unexpected argument '${args[1]}': grok kg ${verb} takes options only`);
  const output: string = argv.output ?? argv.o ?? 'table';
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
  return true;
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
