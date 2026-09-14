/// The Dart batch (build-plan.md WO-6): the node, edge and claim lines `prop_gen` writes to
/// `.kg/batches/kg-dart.jsonl` before a build, and the freshness the manifest reports for them.
/// Dart is never parsed here; an absent batch is `missing`, not an error, and every op says so.
/// "Fresh and readable" is not "complete": one versioned header, an explicit revision and package
/// coverage, and claim identities that resolve are required before a batch may call itself ok
/// (kg-codex-review-3.md #6), and the rows the emitter refuses make it partial after finalize.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {Emitter} from '../emitter';
import {Row} from '../normalize';
import {BuildContext, Extractor} from '../registry';
import {gitRevisions} from '../write';
import {HomeIndex, homesOf, resolveMention} from './markers';

/** Where the generator writes, relative to the monorepo root. */
export const BATCH_FILE = '.kg/batches/kg-dart.jsonl';
const STALE_DAYS = 7;
const RECORDS = ['batch', 'node', 'edge', 'claim'];
/** What the header must say about itself before its payload may be trusted. */
const HEADER_KEYS = ['revision', 'generator', 'packages'];
/** A Dart package of the checkout: a folder under `core/` with a pubspec.yaml, `core/<area>/libs/<pkg>` included. */
const PUBSPECS = ['core/*/*/pubspec.yaml', 'core/*/*/*/pubspec.yaml'];

interface Header {
  built_at?: string;
  revision?: string;
  generator?: string;
  packages?: unknown;
}

export const dartExtractor: Extractor = {
  name: 'dart',
  layer: 'dart',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const file = path.join(ctx.repoRoot, ...BATCH_FILE.split('/'));
    if (!fs.existsSync(file)) {
      emitter.source('dart', 'missing');
      return;
    }
    const index = new HomeIndex(homesOf(ctx));
    const rows = new Map<string, number>();
    let header: Header | undefined;
    let payload = false;
    let invalid = 0;
    const refuse = (at: number, message: string) => {
      invalid++;
      emitter.problem('invalid_rows', `${BATCH_FILE}:${at}: ${message}`);
    };
    const lines = fs.readFileSync(file, 'utf8').split('\n');
    for (let i = 0; i < lines.length; i++) {
      const text = lines[i].trim();
      if (!text) continue;
      const row = parse(text);
      if (!row || !RECORDS.includes(String(row.record))) {
        refuse(i + 1, row ? `unknown record '${row.record}'` : 'not JSON');
        continue;
      }
      const {record, ...rest} = row;
      if (record === 'batch') {
        if (header || payload) refuse(i + 1, header ? 'a second batch record' : 'a batch record after the payload');
        else header = rest as Header;
        continue;
      }
      payload = true;
      if (record === 'node') {
        if (emitter.node(rest).accepted) countRow(rows, rest.path);
      }
      else if (record === 'edge') emitter.edge(rest);
      else {
        const outcome = claim(emitter, index, rest, i + 1);
        if (outcome === 'malformed') refuse(i + 1, 'a claim needs a file, a feature, rung 1 and no mode but participates');
        else if (outcome === 'unresolved') invalid++;
      }
    }
    if (!header) {
      invalid++;
      emitter.problem('invalid_rows', `${BATCH_FILE}: no batch record`);
    }
    for (const key of header ? HEADER_KEYS : [])
      if (header![key as keyof Header] === undefined) {
        invalid++;
        emitter.problem('invalid_rows', `${BATCH_FILE}: the batch record has no ${key}`);
      }
    if (header) emitter.manifest('dart_packages', coverage(ctx.repoRoot, header, rows));
    const stale = header && staleReason(header, ctx.repoRoot);
    if (stale) emitter.problem('stale_sources', `${BATCH_FILE}: ${stale}`);
    emitter.source('dart', invalid ? 'partial' : stale ? 'stale' : 'ok');
  },
};

function parse(text: string): Row | null {
  try {
    const row = JSON.parse(text);
    return row && typeof row === 'object' && !Array.isArray(row) ? row as Row : null;
  }
  catch {
    return null;
  }
}

/**
 * A Dart file marker is rung 1 and nothing else: the generator sees one file at a time and cannot speak for a home
 * document's roots or its prose. `/// ~id` owns, `// ~id` on its own line only participates, and the feature it names
 * has to be one the home documents declare, or the claim is a typo and is counted, never drawn.
 */
function claim(emitter: Emitter, index: HomeIndex, row: Row, at: number): 'ok' | 'malformed' | 'unresolved' {
  if (typeof row.file !== 'string' || typeof row.feature !== 'string') return 'malformed';
  if ((row.rung !== undefined && row.rung !== 1) || (row.mode !== undefined && row.mode !== 'participates')) return 'malformed';
  const target = resolveMention(emitter, index, row.feature, `${BATCH_FILE}:${at}`);
  if (!target || target.root !== 'feature') return 'unresolved';
  emitter.claim({
    file: row.file, feature: target.id, rung: 1, source: 'marker',
    mode: row.mode === 'participates' ? 'participates' : undefined,
    props: row.props && typeof row.props === 'object' ? row.props as Record<string, unknown> : {},
    line: typeof row.line === 'number' ? row.line : at,
  });
  return 'ok';
}

/** A batch generated from another revision of the core repo, or more than a week ago, is still loaded and reported. */
function staleReason(header: Header, repoRoot: string): string | undefined {
  const head = gitRevisions(repoRoot).reddata;
  if (typeof header.revision === 'string' && head !== 'unknown' && header.revision !== head)
    return `generated from ${header.revision.slice(0, 12)}, HEAD is ${head.slice(0, 12)}`;
  const builtAt = Date.parse(String(header.built_at ?? ''));
  const days = (Date.now() - builtAt) / 86400000;
  if (Number.isNaN(builtAt)) return 'no built_at';
  return days > STALE_DAYS ? `generated ${Math.floor(days)} days ago` : undefined;
}

/**
 * What the batch covers, for the manifest: a package the header lists with the rows it carried for it, so `0` reads as
 * "the generator walked it and it holds nothing"; a Dart package of this checkout the header does not list is `omitted`,
 * which reads as "coverage unknown". Without the distinction an absent package and an empty one look alike.
 */
function coverage(repoRoot: string, header: Header, rows: Map<string, number>): Record<string, number | string> {
  const listed = (Array.isArray(header.packages) ? header.packages : []).map(String);
  const known = PUBSPECS.flatMap((p) => globSync(p, {cwd: repoRoot, posix: true, windowsPathsNoEscape: true})).map(packageOf);
  const out: Record<string, number | string> = {};
  for (const name of [...new Set([...listed, ...rows.keys(), ...known])].sort())
    out[name] = listed.includes(name) ? rows.get(name) ?? 0 : 'omitted';
  return out;
}

function countRow(rows: Map<string, number>, file: unknown): void {
  const name = typeof file === 'string' ? packageOf(file) : undefined;
  if (name) rows.set(name, (rows.get(name) ?? 0) + 1);
}

/** `core/shared/ddt/lib/x.dart` and `core/shared/ddt/pubspec.yaml` are both `ddt`; `core/server/libs/shelf/...` is `shelf`. */
function packageOf(file: string): string | undefined {
  const segments = file.split('/');
  if (segments[0] !== 'core' || segments.length < 3) return undefined;
  return segments[2] === 'libs' ? segments[3] : segments[2];
}
