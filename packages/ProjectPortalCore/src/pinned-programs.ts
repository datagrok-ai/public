/** The "My programs" tiles: one per pinned program the user can see, across portfolios. The pins
 * come from the favorites (./pinning); stage, readiness for the next gate and attention from the
 * portfolio table's facts, read for the pinned programs only, in the same batched queries. */
import * as DG from 'datagrok-api/dg';
import {portalDb, ProgramRow} from './generated/db';
import {inBatches, portfolioFacts} from './portfolio-table';
import {myPins} from './pinning';
import {assembleTiles, describeTile, ProgramTile, TileInput} from './model/pinning';

const LIMIT = 10000;

/** The tiles of every pinned program the user can see, across portfolios. The folders come from
 * the favorites; the programs from their space ids (rows the user cannot see give no tile); stage,
 * readiness and attention from the portfolio table's facts, read for those programs only. */
export async function pinnedTiles(): Promise<ProgramTile[]> {
  const pins = await myPins();
  if (!pins.size)
    return [];
  // `space_id in (...)` fails on the server (2026-10-02: a string column of uuid-shaped values is
  // compared as uuid, "operator does not exist: text = uuid"); `=` joined by `or` does not
  const keys = [...pins.keys()];
  const programs: ProgramRow[] = [];
  for (let i = 0; i < keys.length; i += 50) {
    const filter = keys.slice(i, i + 50).map((k) => `space_id = "${k}"`).join(' or ');
    programs.push(...await portalDb.programs.query({filter, limit: LIMIT}));
  }
  if (!programs.length)
    return [];
  const ids = programs.map((p) => p.id);
  const me = DG.User.current().id;
  const [{stages, facts}, follows] = await Promise.all([
    portfolioFacts(undefined, Date.now(), ids),
    inBatches(ids, (b) => portalDb.follows.query({filter: `user_id = "${me}" and program_id in (${b})`, limit: LIMIT})),
  ]);
  const str = (v: any): string | null => v ? String(v.toISOString?.() ?? v) : null;
  const inputs: TileInput[] = facts.map((f) => ({
    programId: f.program.id, code: f.program.code, name: f.program.name, spaceId: f.program.space_id ?? null,
    portfolioCode: f.portfolio.code, portfolioName: f.portfolio.name, status: f.program.status ?? null,
    stages: (stages.get(f.program.stage_model_id) ?? stages.get(f.portfolio.stage_model_id) ?? [])
      .map((s) => ({code: s.code, name: s.name})),
    stage: f.stage?.code ?? null,
    next: f.next ? {code: f.next.code, name: f.next.name} : null,
    readiness: f.readiness ? {score: f.readiness.score, requiredMet: f.readiness.requiredMet,
      requiredTotal: f.readiness.requiredTotal, gaps: f.readiness.gaps, dataError: f.readiness.dataError} : null,
    attention: f.attention?.count ?? null,
    follows: follows.filter((x) => x.program_id === f.program.id && x.active !== false)
      .map((x) => ({newsOn: str(x.news_on), seenOn: str(x.seen_on), lastNews: x.last_news ?? null})),
  }));
  return assembleTiles(inputs, pins);
}

/** One row per pinned program: what its tile shows. */
export function pinnedFrame(tiles: ProgramTile[]): DG.DataFrame {
  const S = DG.COLUMN_TYPE.STRING;
  const col = (type: DG.ColumnType, n: string, f: (t: ProgramTile) => any) => DG.Column.fromList(type, n, tiles.map(f));
  const df = DG.DataFrame.fromColumns([
    col(S, 'id', (t) => t.programId),
    col(S, 'program', (t) => t.code),
    col(S, 'name', (t) => t.name),
    col(S, 'portfolio', (t) => t.portfolioCode),
    col(S, 'stage', (t) => t.stage),
    col(DG.COLUMN_TYPE.INT, 'step', (t) => t.step),
    col(DG.COLUMN_TYPE.INT, 'steps', (t) => t.steps),
    col(S, 'status', (t) => t.status),
    col(S, 'next stage', (t) => t.next),
    col(DG.COLUMN_TYPE.FLOAT, 'readiness', (t) => t.readiness),
    col(DG.COLUMN_TYPE.BOOL, 'ready', (t) => t.ready),
    col(DG.COLUMN_TYPE.INT, 'attention', (t) => t.attention),
    col(S, 'new data', (t) => t.newData),
    col(DG.COLUMN_TYPE.BOOL, 'pinned by me', (t) => t.pinnedByMe),
    col(S, 'pinned for', (t) => t.pinnedFor.join(', ') || null),
    col(S, 'address', (t) => t.address),
    col(S, 'summary', (t) => describeTile(t)),
  ]);
  df.name = 'Pinned programs';
  df.col('readiness')!.meta.format = '0%';
  return df;
}
