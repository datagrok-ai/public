/// The Dart batch (build-plan.md WO-6): the node, edge and claim lines `prop_gen` writes to
/// `.kg/batches/kg-dart.jsonl` before a build, and the freshness the manifest reports for them.
/// Dart is never parsed here; an absent batch is `missing`, not an error, and every op says so.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../emitter';
import {Row} from '../normalize';
import {BuildContext, Extractor} from '../registry';
import {gitRevisions} from '../write';

/** Where the generator writes, relative to the monorepo root. */
export const BATCH_FILE = '.kg/batches/kg-dart.jsonl';
const STALE_DAYS = 7;
const RECORDS = ['batch', 'node', 'edge', 'claim'];

interface Header {
  built_at?: string;
  revision?: string;
  generator?: string;
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
    let header: Header | undefined;
    let invalid = 0;
    const lines = fs.readFileSync(file, 'utf8').split('\n');
    for (let i = 0; i < lines.length; i++) {
      const text = lines[i].trim();
      if (!text) continue;
      const row = parse(text);
      if (!row || !RECORDS.includes(String(row.record))) {
        invalid++;
        emitter.problem('invalid_rows', `${BATCH_FILE}:${i + 1}: ${row ? `unknown record '${row.record}'` : 'not JSON'}`);
        continue;
      }
      const {record, ...rest} = row;
      if (record === 'batch') header = rest as Header;
      else if (record === 'node') emitter.node(rest);
      else if (record === 'edge') emitter.edge(rest);
      else if (!claim(emitter, rest, i + 1)) {
        invalid++;
        emitter.problem('invalid_rows', `${BATCH_FILE}:${i + 1}: a claim needs file and feature`);
      }
    }
    if (!header) {
      invalid++;
      emitter.problem('invalid_rows', `${BATCH_FILE}: no batch record`);
    }
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

/** `/// ~id` on a file or class owns (rung 1); `// ~id` on its own line only participates. */
function claim(emitter: Emitter, row: Row, at: number): boolean {
  if (typeof row.file !== 'string' || typeof row.feature !== 'string') return false;
  const rung = row.rung === 2 || row.rung === 3 ? row.rung : 1;
  emitter.claim({
    file: row.file, feature: row.feature, rung, source: 'marker',
    mode: row.mode === 'participates' ? 'participates' : undefined,
    props: row.props && typeof row.props === 'object' ? row.props as Record<string, unknown> : {},
    line: typeof row.line === 'number' ? row.line : at,
  });
  return true;
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
