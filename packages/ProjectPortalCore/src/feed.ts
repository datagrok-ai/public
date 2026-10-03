/** The activity feed against the platform (F-21): reads one program's records (gate log, evidence,
 * runs, selections, documents, snapshots, compound role changes, the portfolio's sync history and
 * the follows' last news) and turns them into entries; "my feed" does the same for every program the
 * current user follows. Nothing new is stored. The rules live in ./model/feed. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb, PortfolioRow, ProgramRow, SyncRunRow} from './generated/db';
import {portfolioByCode, programData, stagesOf} from './portal';
import {FeedEntry, FeedFilter, FeedRecords, filterFeed, mergeFeeds, programFeedEntries} from './model/feed';
import {ProgramChange} from './model/sync-runs';
import {columnValue} from './model/conditions';
import {rowsOf} from './model/new-data';
import {myFollowsAll} from './new-data';
import {filterValue, ROW_LIMIT} from './rows';

const iso = (d: any): string | null => d ? dayjs(d).toISOString() : null;
const inList = (ids: string[]) => ids.map(filterValue).join(', ');

/** Sync runs read per portfolio: the feed looks this far back in the sync history. */
const SYNC_RUNS = 50;

type SyncChanges = {programs?: ProgramChange[], newData?: {program: string, line: string}[]} | undefined;

/** The portfolio's recent sync runs, read once per feed. */
async function syncRunsOf(portfolioId: string): Promise<SyncRunRow[]> {
  return portalDb.syncRuns.query({filter: `portfolio_id = ${filterValue(portfolioId)}`, sort: '!started_on', limit: SYNC_RUNS})
    .catch(() => []);
}

/** Everything one program's feed is made of. [syncRuns] spares a second read of the portfolio's history. */
export async function feedRecordsOf(program: ProgramRow, portfolio: PortfolioRow, syncRuns?: SyncRunRow[]): Promise<FeedRecords> {
  const by = `program_id = "${program.id}"`;
  const [stages, decisions, evidence, runs, artifacts, documents, snapshots, roleChanges, follows, predictions, syncs] = await Promise.all([
    stagesOf(program.stage_model_id),
    portalDb.stageTransitions.query({filter: by, limit: ROW_LIMIT}),
    portalDb.evidences.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.modelRuns.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.criterionArtifacts.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.documents.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.programSnapshots.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.compoundRoleChanges.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.follows.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.predictions.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    syncRuns ? Promise.resolve(syncRuns) : syncRunsOf(portfolio.id),
  ]);
  const statusIds = [...new Set(artifacts.map((a) => a.criterion_status_id))];
  const statuses = statusIds.length ? await portalDb.criterionStatuses.query({filter: `id in (${inList(statusIds)})`, limit: ROW_LIMIT}).catch(() => []) : [];
  const criterionIds = [...new Set(statuses.map((s) => s.gate_criterion_id))];
  const criteria = criterionIds.length ? await portalDb.gateCriterions.query({filter: `id in (${inList(criterionIds)})`, limit: ROW_LIMIT}).catch(() => []) : [];
  const compoundIds = [...new Set(roleChanges.map((c) => c.compound_id))];
  const compounds = compoundIds.length ? await portalDb.compounds.query({filter: `id in (${inList(compoundIds)})`, limit: ROW_LIMIT}).catch(() => []) : [];
  const runTitle = (id: string) => {
    const r = runs.find((x) => x.id === id);
    return r ? r.title ?? r.run_id ?? r.model_name : 'a run';
  };
  return {
    program: {id: program.id, code: program.code, portfolio: portfolio.code},
    stages: stages.map((s) => ({id: s.id, code: s.code, name: s.name})),
    decisions: decisions.map((t) => ({id: t.id, decision: t.decision, from: t.from_stage_id ?? null, to: t.to_stage_id ?? null,
      by: t.decided_by ?? (t.source === 'approver' ? t.author_id : null), on: iso(t.decided_on ?? t.created_on),
      source: t.source ?? null, rationale: t.rationale ?? null})),
    evidence: evidence.map((e) => ({id: e.id, decisionId: e.stage_transition_id ?? null, kind: e.kind, name: e.name,
      by: e.added_by ?? e.author_id, on: iso(e.added_on ?? e.created_on)})),
    runs: runs.map((x) => ({id: x.id, title: x.title ?? x.run_id ?? x.model_name, model: x.model_name.split(':').pop()!,
      by: x.run_author ?? x.author_id, on: iso(x.created_on)})),
    selections: artifacts.filter((a) => a.selected_on).map((a) => {
      const s = statuses.find((x) => x.id === a.criterion_status_id);
      const c = criteria.find((x) => x.id === s?.gate_criterion_id);
      return {id: a.id, run: runTitle(a.model_run_id), criterion: c?.name ?? 'a criterion', by: a.selected_by ?? null,
        on: iso(a.selected_on), rationale: a.rationale ?? null};
    }),
    predictions: [...new Set(predictions.map((p) => p.model_run_id).filter((x): x is string => !!x))].map((runId) => {
      const mine = predictions.filter((p) => p.model_run_id === runId);
      return {runId, model: mine[0].model_code.split(':').pop()!, version: mine[0].model_version,
        count: mine.length, by: mine[0].predicted_by ?? mine[0].author_id, on: iso(mine[0].predicted_on ?? mine[0].created_on)};
    }),
    documents: documents.map((d) => ({id: d.id, title: d.title, version: d.doc_version ?? null, status: d.status,
      linkedBy: d.linked_by ?? d.author_id, linkedOn: iso(d.linked_on ?? d.created_on), approvedBy: d.approved_by ?? null,
      approvedOn: d.status !== 'draft' ? iso(d.approved_on) : null,
      ...((d.evidence_check as any)?.state === 'stale' && d.status !== 'superseded' ?
        {staleOn: iso((d.evidence_check as any).checkedOn), stale: ((d.evidence_check as any).changes ?? []).join('; ')} : {})})),
    snapshots: snapshots.map((s) => ({id: s.id, label: s.label, by: s.taken_by ?? s.author_id, on: iso(s.taken_on)})),
    roleChanges: roleChanges.map((c) => ({id: c.id, compound: compounds.find((x) => x.id === c.compound_id)?.code ?? c.compound_id,
      from: c.from_role ?? null, to: c.to_role, by: c.changed_by ?? null, on: iso(c.changed_on ?? c.created_on)})),
    syncs: syncs.map((s) => {
      const ch = s.changes as SyncChanges;
      return {id: s.id, on: iso(s.finished_on ?? s.started_on), by: s.run_by ?? null,
        change: ch?.programs?.find((p) => p.program === program.code) ?? null,
        newData: (ch?.newData ?? []).filter((n) => n.program === program.code).map((n) => n.line)};
    }).filter((s) => s.change || s.newData.length),
    news: follows.filter((f) => f.last_news && f.news_on).map((f) => ({id: f.id, line: f.last_news!, on: iso(f.news_on)})),
  };
}

/** One program's feed, newest first. */
export async function programFeedOf(portfolio: PortfolioRow, program: ProgramRow, filter: FeedFilter = {}): Promise<FeedEntry[]> {
  const all = mergeFeeds(programFeedEntries(await feedRecordsOf(program, portfolio)));
  if (filter.series && !filter.seriesCompounds)
    filter = {...filter, seriesCompounds: await seriesCompounds(portfolio, program, filter.series)};
  return filterFeed(all, filter);
}

/** The compounds of one series, from the program dataset's `series` and `compound` columns. */
export async function seriesCompounds(portfolio: PortfolioRow, program: ProgramRow, series: string): Promise<string[]> {
  const rows = await programData(portfolio, program.code).catch(() => null);
  return [...new Set(rowsOf(rows ?? [], series).map((r) => String(columnValue(r, 'compound') ?? '')).filter((c) => c.length > 0))];
}

/** Every program the current user follows, or one of whose work units they follow: their feeds as
 * one, newest first. A series follow narrows its program's feed to the series unless the program is
 * followed as a whole too. */
export async function myFeedOf(filter: FeedFilter = {}): Promise<{entries: FeedEntry[], programs: number}> {
  const follows = await myFollowsAll();
  const ids = [...new Set(follows.map((f) => f.program_id))];
  if (!ids.length)
    return {entries: [], programs: 0};
  const programs = await portalDb.programs.query({filter: `id in (${inList(ids)})`, limit: ROW_LIMIT});
  const portfolios = await portalDb.portfolios.query({filter: `id in (${inList([...new Set(programs.map((p) => p.portfolio_id))])})`, limit: ROW_LIMIT});
  const syncs = new Map<string, Promise<SyncRunRow[]>>();
  const lists = await Promise.all(programs.map(async (p) => {
    const pf = portfolios.find((x) => x.id === p.portfolio_id);
    // a portfolio not shared with the user: its programs' feeds are not theirs to read
    if (!pf)
      return [];
    if (!syncs.has(pf.id))
      syncs.set(pf.id, syncRunsOf(pf.id));
    const entries = programFeedEntries(await feedRecordsOf(p, pf, await syncs.get(pf.id)));
    const mine = follows.filter((f) => f.program_id === p.id);
    if (mine.some((f) => !f.work_unit_code))
      return entries;
    // only series follows: what is about those series
    const parts = await Promise.all(mine.filter((f) => f.work_unit_kind === 'series').map(async (f) =>
      filterFeed(entries, {series: f.work_unit_code, seriesCompounds: await seriesCompounds(pf, p, f.work_unit_code!)})));
    return parts.flat();
  }));
  return {entries: filterFeed(mergeFeeds(...lists), filter), programs: programs.length};
}

/** Display names of the users in [entries], by id; a user that cannot be read keeps the id. */
export async function feedUserNames(entries: FeedEntry[]): Promise<Map<string, string>> {
  const ids = [...new Set(entries.map((e) => e.who).filter((x): x is string => !!x))];
  const names = await Promise.all(ids.map((id) => grok.dapi.users.find(id).then((u) => u.friendlyName || u.login).catch(() => id)));
  return new Map(ids.map((id, i) => [id, names[i]]));
}

/** The feed as a table: the stable contract of programFeed and myFeed. */
export function feedFrame(entries: FeedEntry[], names: Map<string, string>, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'at', entries.map((e) => new Date(e.at))),
    DG.Column.fromStrings('portfolio', entries.map((e) => e.portfolio ?? '')),
    DG.Column.fromStrings('program', entries.map((e) => e.program)),
    DG.Column.fromStrings('kind', entries.map((e) => e.kind)),
    DG.Column.fromStrings('what', entries.map((e) => e.text)),
    DG.Column.fromStrings('who', entries.map((e) => e.who ? names.get(e.who) ?? e.who : '')),
    DG.Column.fromStrings('work unit', entries.map((e) => e.unit ?? '')),
    DG.Column.fromStrings('record', entries.map((e) => e.ref ? `${e.ref.kind}:${e.ref.id}` : '')),
  ]);
  df.name = name;
  return df;
}

/** A program by portfolio and program code, under row security. */
export async function programOf(portfolioCode: string, programCode: string): Promise<{portfolio: PortfolioRow, program: ProgramRow}> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(programCode)}`, limit: ROW_LIMIT}))[0];
  if (!program)
    throw new Error(`${portfolioCode} has no program ${programCode}, or it is not shared with you`);
  return {portfolio, program};
}
