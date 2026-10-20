/** The activity feed against the platform (F-21): reads one program's records (gate log, evidence,
 * runs, selections, documents, snapshots, compound role changes, the portfolio's sync history and
 * the follows' last news, the discussion's posts) and turns them into entries; "my feed" does the same for every program the
 * current user follows. Nothing new is stored. The rules live in ./model/feed. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {textColumn} from './columns';
import dayjs from 'dayjs';
import {FollowingRow, portalDb, PortfolioRow, ProgramRow, SyncLogRow} from './generated/db';
import {portfolioByCode, programData, stagesOf} from './portal';
import {FeedEntry, FeedFilter, FeedRecords, filterFeed, mergeFeeds, programFeedEntries} from './model/feed';
import {ProgramChange} from './model/sync-runs';
import {columnValue} from './model/conditions';
import {rowsOf} from './model/new-data';
import {myFollowsAll} from './new-data';
import {reviewsOf} from './model/approvals';
import {newestFirstBy} from './model/time';
import {filterValue, queryIn, ROW_LIMIT} from './rows';
import {userNames} from './users';
import {myPortfolioIds} from './secured-rows';
import {postsOf} from './discussion';
import {postFeedText} from './model/discussion';

const iso = (d: any): string | null => d ? dayjs(d).toISOString() : null;

/** Sync runs read per portfolio: the feed looks this far back in the sync history. */
const SYNC_RUNS = 50;

type SyncChanges = {programs?: ProgramChange[], newData?: {program: string, line: string}[]} | undefined;

/** The portfolio's recent sync runs, read once per feed. */
async function syncRunsOf(portfolioId: string): Promise<SyncLogRow[]> {
  const mine = await myPortfolioIds().catch(() => new Set<string>());
  if (mine && !mine.has(portfolioId))
    return [];
  return portalDb.syncLogs.query({filter: `portfolio_id = ${filterValue(portfolioId)}`, sort: '!started_on', limit: SYNC_RUNS})
    .catch(() => []);
}

/** Everything one program's feed is made of. [syncRuns] spares a second read of the portfolio's history. */
export async function feedRecordsOf(program: ProgramRow, portfolio: PortfolioRow, syncRuns?: SyncLogRow[]): Promise<FeedRecords> {
  const by = `program_id = "${program.id}"`;
  const [stages, decisions, evidence, runs, artifacts, documents, snapshots, roleChanges, follows, predictions, syncs, observations] = await Promise.all([
    stagesOf(program.stage_model_id),
    portalDb.stageTransitions.query({filter: by, limit: ROW_LIMIT}),
    portalDb.evidences.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.modelRuns.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.criterionArtifacts.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.documents.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.programSnapshots.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    portalDb.compoundRoleChanges.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    // the news told to the current user: another user's follows are theirs alone
    portalDb.followings.query({filter: `${by} and user_id = ${filterValue(DG.User.current().id)}`, limit: ROW_LIMIT}).catch(() => []),
    portalDb.predictions.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
    syncRuns ? Promise.resolve(syncRuns) : syncRunsOf(portfolio.id),
    portalDb.observations.query({filter: by, limit: ROW_LIMIT}).catch(() => []),
  ]);
  const statusIds = [...new Set(artifacts.map((a) => a.criterion_status_id))];
  const statuses = await queryIn(portalDb.criterionStatuses, 'id', statusIds).catch(() => []);
  const criterionIds = [...new Set(statuses.map((s) => s.gate_criterion_id))];
  const criteria = await queryIn(portalDb.gateCriterions, 'id', criterionIds).catch(() => []);
  // compound codes from the program's own links; the compound table is closed to most users
  const compoundIds = [...new Set(roleChanges.map((c) => c.compound_id))];
  const links = compoundIds.length ? await portalDb.programCompounds.query({filter: by, limit: ROW_LIMIT}).catch(() => []) : [];
  const codes = new Map(links.filter((l) => l.compound_code).map((l) => [l.compound_id, l.compound_code!]));
  const rest = compoundIds.filter((id) => !codes.has(id));
  for (const c of await queryIn(portalDb.compounds, 'id', rest).catch(() => []))
    codes.set(c.id, c.code);
  const [programRoles, modelCards, posts] = await Promise.all([
    programRolesOf(program.id).catch(() => []),
    approvedCardsOf([...runs.map((x) => x.model_id), ...predictions.map((x) => x.model_id)]).catch(() => []),
    // the discussion is the platform chat on the program's folder (F-51)
    postsOf(program).catch(() => []),
  ]);
  const runTitle = (id: string) => {
    const r = runs.find((x) => x.id === id);
    return r ? r.title ?? r.run_id ?? r.model_name : 'a run';
  };
  return {
    program: {id: program.id, code: program.code, portfolio: portfolio.code},
    stages: stages.map((s) => ({id: s.id, code: s.code, name: s.name})),
    decisions: decisions.map((t) => ({id: t.id, decision: t.decision, from: t.from_stage_id ?? null, to: t.to_stage_id ?? null,
      by: t.decided_by ?? (t.source === 'approver' ? t.author_id : null), on: iso(t.decided_on ?? t.created_on),
      source: t.source ?? null, rationale: t.rationale ?? null, outcome: t.outcome ?? null})),
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
    roleChanges: roleChanges.map((c) => ({id: c.id, compound: codes.get(c.compound_id) ?? c.compound_id,
      from: c.from_role ?? null, to: c.to_role, by: c.changed_by ?? null, on: iso(c.changed_on ?? c.created_on)})),
    syncs: syncs.map((s) => {
      const ch = s.changes as SyncChanges;
      return {id: s.id, on: iso(s.finished_on ?? s.started_on), by: s.run_by ?? null,
        change: ch?.programs?.find((p) => p.program === program.code) ?? null,
        newData: (ch?.newData ?? []).filter((n) => n.program === program.code).map((n) => n.line)};
    }).filter((s) => s.change || s.newData.length),
    news: follows.filter((f) => f.last_news && f.news_on).map((f) => ({id: f.id, line: f.last_news!, on: iso(f.news_on)})),
    programRoles,
    modelCards,
    observations: observations.map((o) => ({id: o.id, subject: o.subject_kind, code: o.subject_code ?? null, title: o.title,
      by: o.observed_by ?? o.author_id, on: iso(o.observed_on ?? o.created_on),
      retractedBy: o.retracted_by ?? null, retractedOn: o.status === 'retracted' ? iso(o.retracted_on) : null,
      retraction: o.retraction ?? null})),
    discussion: posts.map((p) => ({id: p.id, on: p.on, by: p.by, text: postFeedText(p), question: p.question})),
  };
}

/** The program's roles as they stand, with the names of their groups. */
async function programRolesOf(programId: string): Promise<NonNullable<FeedRecords['programRoles']>> {
  const rows = await portalDb.programRoles.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT});
  const names = new Map<string, string>();
  for (const id of new Set(rows.map((r) => r.group_id))) {
    const g = await grok.dapi.groups.find(id).catch(() => null);
    names.set(id, g ? g.friendlyName || g.name : id);
  }
  return rows.map((r) => ({id: r.id, group: names.get(r.group_id) ?? r.group_id, role: r.role, active: r.active !== false,
    by: r.granted_by ?? r.author_id ?? null, on: iso(r.granted_on ?? r.created_on)}));
}

/** The approved cards of [modelIds]: the model, the version, who approved it (the request's last
 * reviewer, when the approval workflow approved it) and when. */
async function approvedCardsOf(modelIds: (string | null | undefined)[]): Promise<NonNullable<FeedRecords['modelCards']>> {
  const ids = [...new Set(modelIds.filter((x): x is string => !!x))];
  if (!ids.length)
    return [];
  const [models, cards] = await Promise.all([
    queryIn(portalDb.models, 'id', ids),
    queryIn(portalDb.modelCards, 'model_id', ids),
  ]);
  const approved = cards.filter((c) => c.approval === 'approved' && c.approved_on);
  const requests = await queryIn(portalDb.approvals, 'target', approved.map((c) => c.id),
    'kind = "model_card" and status = "approved"').catch(() => []);
  return approved.map((c) => {
    const reviews = reviewsOf(requests.filter((r) => r.target === c.id).sort(newestFirstBy((r) => r.closed_on))[0]?.reviews);
    return {id: c.id, model: models.find((m) => m.id === c.model_id)?.name ?? c.model_id, version: c.model_version,
      by: reviews.length ? reviews[reviews.length - 1].by : null, on: iso(c.approved_on)};
  });
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
  return feedOfFollows(await myFollowsAll(), filter);
}

/** The feed of [follows], as `myFeedOf` reads it. */
export async function feedOfFollows(follows: FollowingRow[], filter: FeedFilter = {}): Promise<{entries: FeedEntry[], programs: number}> {
  const ids = [...new Set(follows.map((f) => f.program_id))];
  if (!ids.length)
    return {entries: [], programs: 0};
  // a followed program may be gone or not visible to the user: then nothing is left to read
  const programs = await queryIn(portalDb.programs, 'id', ids);
  const portfolios = await queryIn(portalDb.portfolios, 'id', programs.map((p) => p.portfolio_id));
  const syncs = new Map<string, Promise<SyncLogRow[]>>();
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
export function feedUserNames(entries: FeedEntry[]): Promise<Map<string, string>> {
  return userNames(entries.map((e) => e.who));
}

/** The feed as a table: the stable contract of programFeed and myFeed. */
export function feedFrame(entries: FeedEntry[], names: Map<string, string>, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'at', entries.map((e) => new Date(e.at))),
    textColumn('portfolio', entries.map((e) => e.portfolio ?? '')),
    textColumn('program', entries.map((e) => e.program)),
    textColumn('kind', entries.map((e) => e.kind)),
    textColumn('what', entries.map((e) => e.text)),
    textColumn('who', entries.map((e) => e.who ? names.get(e.who) ?? e.who : '')),
    textColumn('work unit', entries.map((e) => e.unit ?? '')),
    textColumn('record', entries.map((e) => e.ref ? `${e.ref.kind}:${e.ref.id}` : '')),
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
