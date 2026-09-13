/// `grok kg check|gen` — the knowledge-graph type files and home documents
/// (core/docs/knowledge-graph/CONVENTIONS.md §11.1).
import * as fs from 'fs';
import * as path from 'path';
import {loadTypeSystem, Issue} from '../utils/kg/types';
import {loadHomes, makeReport, Report} from '../utils/kg/homes';
import {generate, writeOutputs} from '../utils/kg/gen';
import {HELP_KG} from './help';

const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const PLANNED_VERBS = ['build', 'query', 'report'];

export async function kg(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1).map(String);
  const verb: string | undefined = args[0];
  if (!verb || verb === 'help' || argv.help) {
    console.log(HELP_KG);
    return true;
  }
  if (PLANNED_VERBS.includes(verb)) return fail(`grok kg ${verb} is not implemented yet (CONVENTIONS.md §11.1); check and gen are`);
  if (verb !== 'check' && verb !== 'gen') {
    console.error(`unknown verb '${verb}'`);
    return false;
  }
  if (args.length > 1) return fail(`unexpected argument '${args[1]}': grok kg ${verb} takes options only`);
  const output: string = argv.output ?? argv.o ?? 'table';
  if (output !== 'table' && output !== 'json') return fail(`--output must be table or json, got '${output}'`);
  const quiet = argv.quiet === true;
  const typesOnly = argv['types-only'] === true;
  if (verb === 'gen' && typesOnly) return fail('--types-only cannot be combined with gen: FEATURES.md is generated from the home documents');

  const kgRoot = argv.kg ? path.resolve(String(argv.kg)) : findKgRoot(process.cwd());
  if (!kgRoot || !fs.existsSync(path.join(kgRoot, 'schema.yaml')))
    return fail(kgRoot ? `${slashes(kgRoot)}: no schema.yaml` : 'grok kg currently needs the monorepo: run it inside a checkout with core/docs/knowledge-graph/schema.yaml and public/, or pass --kg <dir>');
  const repoRoot = path.resolve(kgRoot, '..', '..', '..');

  const system = loadTypeSystem(kgRoot);
  const homes = typesOnly ? null : loadHomes(system, repoRoot);
  const report = makeReport(system, homes);
  if (verb === 'gen') {
    const generated = generate(system, kgRoot, repoRoot, homes);
    report.errors.push(...generated.errors);
    if (report.errors.length)
      report.warnings.push({file: kgRoot, message: 'nothing generated while check reports errors'});
    else {
      const result = writeOutputs(generated.outputs, argv.check === true);
      report.written = result.written;
      report.stale = result.stale;
      for (const file of result.stale)
        report.errors.push({file, message: 'stale: run grok kg gen'});
    }
  }
  for (const issue of [...report.errors, ...report.warnings]) issue.file = rel(issue.file, repoRoot);
  report.written = report.written?.map((f) => rel(f, repoRoot));
  report.stale = report.stale?.map((f) => rel(f, repoRoot));
  print(report, output, quiet);
  if (report.errors.length) process.exitCode = 1;
  return true;
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

function print(report: Report, output: string, quiet: boolean): void {
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
    `${report.scanned} markdown files scanned, ${total} home document${total === 1 ? '' : 's'}${homes ? ` (${homes})` : ''}; ` +
    `${report.unresolvedExternal} unresolved external reference${report.unresolvedExternal === 1 ? '' : 's'}; ` +
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
