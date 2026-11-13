/**
 * NCA Studies URL router — the single `metaUrl` parse + 3-state dispatch
 * + the ONLY writer of `view.path`.
 *
 * @see help/develop/.../apps/routing.md (pathSegments[0]/[1] dispatch idiom)
 *
 * The app entry `ncaApp(path)` (package.ts) delegates here. The path string is
 * the trailing segment(s) after `/apps/Nca` supplied by the `path {meta.url:true}`
 * parameter; we split it into three states:
 *   - `[]`            → Projects landing (the Hub view, AC-A5 empty-path branch)
 *   - `[studyId]`   → that project's cached two-tier overview
 *   - `[pid, aid]`    → that analysis (lazy compute, AD-4)
 *   - `['analyses', label]` → a STANDALONE analysis (nca-workflow-stages F-S2.6); after a move
 *     into a study the same link redirects there (via the entry's `origin`)
 *
 * **AC-A5 (path-aware hub-dedup, both directions).** The dispatcher runs the
 * router FIRST for a non-empty path (→ project/analysis), and the hub-dedup
 * short-circuit fires ONLY for the empty-path landing (returns the already-open
 * hub instead of a duplicate). Not a blanket "router always first" — that would
 * silently break hub-tab consolidation for the landing case (FM-4 / R1 SY-F6).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {NcaAppView} from './nca-app-view';
import {STANDALONE_ROUTE_SEGMENT, NCA_APP_ROOT, getNcaRepository} from '../adapter';
import type {StudyHeader} from '../adapter';
import {renderStudyOverview} from './nca-study-overview';
import {openAnalysis} from './nca-analysis-open';

/** The app's URL root — defined in the adapter (the headless pipeline builds links too) and
 *  re-exported here, where every link builder lives. */
export {NCA_APP_ROOT};

/** Parsed router state — a pure string→state transform (AC-A1). The segments are
 *  human-readable **labels** (the URL is keyed on `label`, not uuid — see the
 *  module header / datagrok-app-tree-pattern.md §8), resolved to a study/analysis
 *  by {@link routeNcaApp}. */
export type RouteState =
  | {readonly kind: 'landing'}
  | {readonly kind: 'overview'; readonly studyLabel: string}
  | {readonly kind: 'analysis'; readonly studyLabel: string; readonly analysisLabel: string}
  | {readonly kind: 'standalone-analysis'; readonly label: string};

/** Decode one path segment, tolerating a value that is already decoded or carries
 *  a stray `%` (a malformed escape throws — fall back to the raw segment). */
function decodeSeg(seg: string): string {
  try {
    return decodeURIComponent(seg);
  } catch {
    return seg;
  }
}

/**
 * Parse the `metaUrl` path into a {@link RouteState} (AC-A1). Pure — no platform
 * reads. Path segments are `encodeURIComponent`-encoded labels (written by
 * {@link updateNcaStudyPath}); they are decoded back to raw labels here.
 * `analysisLabelParam` is the optional query-param fallback (AC-A3 / DG-1): if the
 * path has exactly one segment but an `analysisId` query param is present
 * (`/apps/Nca/<study>?analysisId=<analysis>`), it routes to the analysis state via
 * the SAME state machine. The query value arrives already URL-decoded (from
 * `URLSearchParams`), so it is NOT decoded again. Trailing slashes and empty
 * segments are tolerated.
 */
export function parseStudyPath(
  path: string | null,
  analysisLabelParam?: string | null,
): RouteState {
  const segments = (path ?? '').split('/').filter((s) => s.length > 0).map(decodeSeg);
  if (segments.length === 0)
    return {kind: 'landing'};
  // `analyses/<label>` — a standalone analysis (F-S2.6). `analyses` alone → the hub (its
  // Ungrouped section lists them). The segment is a reserved study label, so this cannot
  // shadow a study.
  if (segments[0] === STANDALONE_ROUTE_SEGMENT)
    return segments.length >= 2 ? {kind: 'standalone-analysis', label: segments[1]} : {kind: 'landing'};
  if (segments.length === 1) {
    const aid = analysisLabelParam ?? null;
    return aid !== null && aid.length > 0 ?
      {kind: 'analysis', studyLabel: segments[0], analysisLabel: aid} :
      {kind: 'overview', studyLabel: segments[0]};
  }
  // Two or more segments: first is the study, second the analysis. Extra
  // segments (none defined in v1) are ignored.
  return {kind: 'analysis', studyLabel: segments[0], analysisLabel: segments[1]};
}

/** The app-relative tail of a study / study-analysis link (`/<study>[/<analysis>]`), keyed on
 *  **labels**, each `encodeURIComponent`-encoded to match the router's per-segment decode (a
 *  literal `/` in a label stays inside ONE segment). Empty study → `''` (the app root). */
function studyTail(studyLabel?: string | null, analysisLabel?: string | null): string {
  if (studyLabel === undefined || studyLabel === null || studyLabel.length === 0)
    return '';
  const study = `/${encodeURIComponent(studyLabel)}`;
  return analysisLabel !== undefined && analysisLabel !== null && analysisLabel.length > 0 ?
    `${study}/${encodeURIComponent(analysisLabel)}` : study;
}

/** The app-relative tail of a STANDALONE analysis link (`/analyses/<label>`, F-S2.6). */
function standaloneTail(label: string): string {
  return `/${STANDALONE_ROUTE_SEGMENT}/${encodeURIComponent(label)}`;
}

/**
 * Write a view's URL. The platform builds the URL as `view.basePath + view.path` — a PLAIN
 * concatenation (probed on local 2026-10-06: `/apps/Nca` + `/apps/Nca/x` → `/apps/Nca/apps/Nca/x`,
 * no prefix stripping), so `view.path` carries ONLY the app-relative tail; writing the absolute
 * `/apps/NCA/…` there was the doubled-prefix bug (`/apps/Nca/apps/NCA/analyses/…`). A view opened
 * inside the app already has the platform's root (`/apps/Nca`, or `/browse/apps/Nca` in Browse) —
 * kept as is. Any other base is replaced by {@link NCA_APP_ROOT}: null on a fresh `View.fromRoot` /
 * detached preview (the URL would lose the `/apps/Nca` prefix), another app's root, or a path
 * deeper than the app root (not observed — the platform sets the app root itself). The canonical
 * root always yields a URL the router resolves; at worst a Browse prefix is dropped.
 *
 * A TABLE view differs (BUG-20, core `table_view.dart` `set path`): its `path` setter resets the
 * custom `basePath` to `''`, so a root written first is lost and the address bar reads
 * `/analyses/<label>`. Every analysis is a table view, so there the root rides inside `path`.
 * But not on every table view: one the app itself opened (Hub → Run) has `basePath` = the app
 * root after the write and prepends it (measured 2026-11-02 on the local stand: `/apps/Nca` +
 * `/apps/Nca/analyses/x` → `/apps/Nca/apps/Nca/analyses/x`; writing `/analyses/x` → one prefix).
 * So the path is re-written as the tail when the root survived the write.
 */
function writeViewPath(view: DG.ViewBase, tail: string): void {
  const base = view.basePath ?? '';
  const isRoot = (b: string): boolean => b.toLowerCase().endsWith(NCA_APP_ROOT.toLowerCase());
  const root = isRoot(base) ? base : NCA_APP_ROOT;
  if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
    view.path = `${root}${tail}`;
    if (isRoot(view.basePath ?? '')) view.path = tail;
    return;
  }
  if (root !== base)
    view.basePath = root;
  view.path = tail;
}

/**
 * THE single writer of `view.path` (AC-A2). Every in-app navigation —
 * landing→study, study→analysis, in-app drill-in — calls this; no other code
 * assigns `view.path` in the projects subsystem. The resulting URL is
 * `<app root>/<study>[/<analysis>]`, so a copy-paste of the browser URL round-trips.
 * Resolution back to a study is {@link routeNcaApp}'s job.
 */
export function updateNcaStudyPath(
  view: DG.ViewBase,
  studyLabel?: string | null,
  analysisLabel?: string | null,
): void {
  writeViewPath(view, studyTail(studyLabel, analysisLabel));
}

/**
 * THE single source of the absolute deep-link path (`/apps/Nca/<study>[/<analysis>]`) that
 * every copy-link button puts after `window.location.origin`. Built from the SAME tail as the
 * `view.path` writer, so a copy-link can never drift from the address bar or from what
 * `routeNcaApp` resolves (the contract-drift that shipped a uuid-keyed analysis copy-link under
 * label routing). Empty study → `/apps/Nca`.
 */
export function ncaDeepLinkPath(studyLabel?: string | null, analysisLabel?: string | null): string {
  return `${NCA_APP_ROOT}${studyTail(studyLabel, analysisLabel)}`;
}

/** Write the URL of a STANDALONE analysis (`<app root>/analyses/<label>`, F-S2.6). The
 *  standalone counterpart of {@link updateNcaStudyPath} — the only other `view.path` writer. */
export function updateNcaAnalysisPath(view: DG.ViewBase, label: string): void {
  writeViewPath(view, standaloneTail(label));
}

/** The absolute deep-link path of a standalone analysis (`/apps/Nca/analyses/<label>`) — shared
 *  by the hub / tree Copy-link actions; same tail as the path writer, so a copied link always
 *  matches the address bar and what the router resolves. */
export function ncaStandaloneLinkPath(label: string): string {
  return `${NCA_APP_ROOT}${standaloneTail(label)}`;
}

/** The Hub landing view, with AC-A5 hub-dedup: return the already-open hub instead
 *  of a duplicate. Shared by the empty-path branch and the stale-label fallback. */
async function landingView(): Promise<DG.ViewBase> {
  const existing = NcaAppView.findHubView();
  if (existing) {
    grok.shell.v = existing;
    return existing;
  }
  const app = new NcaAppView();
  await app.render();
  return app.view;
}

/**
 * Dispatch an `ncaApp(path)` invocation (AC-A1 / AC-A5). Empty path → hub-dedup
 * landing; non-empty → resolve the label segment(s) to a study/analysis and route
 * there. Labels are unique by construction (create/rename gate), so the lookup is
 * unambiguous; a stale label (study renamed/deleted since the link was made) falls
 * back to the landing hub, and a stale analysis label degrades to the study
 * overview — the accepted rename trade-off of label-keyed URLs.
 */
export async function routeNcaApp(path: string | null): Promise<DG.ViewBase> {
  // Query-param fallback (DG-1): read `analysisId` from the live URL so a
  // depth-1-path + query-param deep-link still resolves to the analysis state.
  //
  // Two signals, each authoritative for one thing:
  //   - the `path` PARAM carries INTENT: empty = open the app / Hub (landing),
  //     non-empty = a deep-link. But it is PRE-DECODED by the platform (a label's
  //     `%2F` arrives as a literal `/`, mis-splitting "PK/TK study"), so it is lossy.
  //   - `window.location.pathname` PRESERVES percent-encoding (lossless) but can be
  //     STALE — after opening a study it still points at that study, so a fresh
  //     app-open (Hub intent, empty param) would wrongly re-resolve to it.
  // So: empty param → Hub (ignore the stale URL); non-empty param → deep-link, and
  // THEN prefer the encoding-preserving window.location tail (fixes %2F), falling
  // back to the param for non-browser contexts. Segments decode in `parseStudyPath`.
  let analysisLabelParam: string | null = null;
  let locTail = '';
  try {
    analysisLabelParam = new URLSearchParams(window.location.search).get('analysisId');
    const m = window.location.pathname.match(/\/apps\/[^/]+\/(.*)$/i);
    locTail = m && m[1] ? m[1] : '';
  } catch {/* non-browser context — primary path scheme only */}

  const rawPath = path ?? '';
  const effectivePath = rawPath.length > 0 ? (locTail.length > 0 ? locTail : rawPath) : '';
  const state = parseStudyPath(effectivePath, analysisLabelParam);

  if (state.kind === 'landing')
    return landingView();

  // AC-A5 non-empty branch: router routes FIRST (no hub-dedup interception).
  const target = resolveDeepLink(state, await getNcaRepository().listStudies());
  if (target.kind === 'overview')
    return renderStudyOverview(target.studyId);
  if (target.kind === 'analysis')
    return openAnalysis(target.analysisId);
  if (target.notice !== undefined)
    grok.shell.info(target.notice);
  return landingView();
}

/** Where a deep link lands, by id — the view-free half of {@link routeNcaApp}. A `landing` target
 *  may carry the toast that says why the link did not resolve. */
export type RouteTarget =
  | {readonly kind: 'landing'; readonly notice?: string}
  | {readonly kind: 'overview'; readonly studyId: string}
  | {readonly kind: 'analysis'; readonly analysisId: string};

/**
 * Resolve a parsed deep link among the caller's visible studies (pure — the test seam for the
 * ambiguous branches). A study link: no study of that label → landing (stale / renamed / deleted —
 * the rename trade-off); several → landing with a toast (SD-26: sharing can bring a second study of
 * the same name into view — never a guess); a stale analysis label → the study overview.
 *
 * An `analyses/<label>` link (F-S2.6): a standalone analysis with that label opens. A miss means it
 * was moved into a study (D5) or deleted: the visible studies are scanned for analyses moved from
 * that label (`movedFromLabel`, SD-22) — exactly one → open it there (the copied link survives the
 * move, AC-S2.6.1); several → landing with a toast; none → landing.
 */
export function resolveDeepLink(state: RouteState, studies: readonly StudyHeader[]): RouteTarget {
  if (state.kind === 'landing')
    return {kind: 'landing'};
  if (state.kind === 'standalone-analysis') {
    const label = state.label;
    const own = studies.find((s) => s.kind === 'standalone' && s.analyses.length > 0 && s.analyses[0].label === label);
    if (own !== undefined)
      return {kind: 'analysis', analysisId: own.analyses[0].id};
    const matches = studies.filter((s) => s.kind === 'study')
      .flatMap((s) => s.analyses.filter((a) => a.movedFromLabel === label));
    if (matches.length === 1)
      return {kind: 'analysis', analysisId: matches[0].id};
    return matches.length > 1 ?
      {kind: 'landing', notice: `'${label}' was moved into ${matches.length} studies — open it from the hub.`} :
      {kind: 'landing'};
  }
  const hits = studies.filter((s) => s.kind === 'study' && s.label === state.studyLabel);
  if (hits.length === 0)
    return {kind: 'landing'};
  if (hits.length > 1)
    return {kind: 'landing', notice: `More than one study is named '${state.studyLabel}'. Open it from the Hub.`};
  const study = hits[0];
  if (state.kind === 'overview')
    return {kind: 'overview', studyId: study.id};
  const analysis = study.analyses.find((a) => a.label === state.analysisLabel);
  return analysis === undefined ? {kind: 'overview', studyId: study.id} : {kind: 'analysis', analysisId: analysis.id};
}
