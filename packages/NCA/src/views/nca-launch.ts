/**
 * NCA launch helpers — pure functions that compose the analysis flow.
 *
 * Three exports, one responsibility each:
 * - `loadSampleByName(file)` — fetch a bundled fixture CSV from the package
 *   AppData path and parse to DG.DataFrame. Source-of-truth for the file path
 *   convention (`System:AppData/Nca/<file>`); see BUG-01-hotfix commit 6822fd9
 *   for the rationale on Pascal-case `Nca`.
 * - `resolveAnalysisName(requested, taken)` — verbatim if a name is supplied
 *   (callers pass the table name); else first free name in `'New NCA'`, `'New NCA (2)'`, ...
 *   per spec §4 decision 4.
 * - `launchAnalysis(df, requestedName)` — orchestrator: resolve name → add
 *   table view → mount mapping pane. Returns the new TableView.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {mountMappingPane} from './nca-mapping-pane';
import type {MountMappingPaneOptions, MappingPaneHandle} from './nca-mapping-pane';
import type {SourcePointer} from '../adapter';
import {normalizeNumericNaMarkers} from '../adapter/na-normalize';
import {reparseFloatColumnsToFloat64} from '../adapter/csv-float64';
import {DEMO_TABLE_FILES, DEMO_ID_TAG, PROVENANCE_TAG} from './nca-demo-tables';
import {keepBrowseInFront} from './shell-panels';

/** AppData segment for the published package; case-sensitive (BUG-01 hotfix). */
const APP_DATA_PREFIX = 'System:AppData/Nca';

/**
 * Load a bundled fixture CSV into a fresh `DG.DataFrame`.
 *
 * Throws (rejects) with a descriptive `Error` if the file is unreachable or
 * parse fails; callers should surface via `grok.shell.error`. The DataFrame's
 * `name` property is set to the filename stem (e.g. `'theoph'`) for cosmetic
 * use; the *view*'s name is set independently by `launchAnalysis`.
 */
export async function loadSampleByName(file: string): Promise<DG.DataFrame> {
  const path = `${APP_DATA_PREFIX}/${file}`;
  const csv = await grok.dapi.files.readAsText(path);
  const df = DG.DataFrame.fromCsv(csv);
  // BUG-01: recover Float64 precision lost by DG's Float32 CSV parse, while we
  // still hold the source text (bundled-sample path). See adapter/csv-float64.ts.
  reparseFloatColumnsToFloat64(df, csv);
  // Filename STEM — directory stripped as well as extension. `file` may be a
  // nested path (`demo/indomethacin_iv.csv`), and this line's own contract says
  // "the filename stem (e.g. 'theoph')", which a raw extension-strip would break
  // by yielding 'demo/indomethacin_iv'. Behaviour-neutral for flat callers.
  const stem = file.replace(/\.csv$/i, '');
  df.name = stem.slice(stem.lastIndexOf('/') + 1);
  return df;
}

/** Read a frame tag, guarded against all three "unset" spellings. `DataFrame.
 *  getTag` is *documented* `string | null`, but the Column-side twin returns
 *  `undefined` for an unset tag and that asymmetry has burned this codebase
 *  once (a literal "undefined undefined" grid header). The guard costs nothing. */
function readFrameTag(df: DG.DataFrame, tag: string): string | null {
  const v = df.getTag(tag) as string | null | undefined;
  return (v === null || v === undefined || v === '') ? null : v;
}

/**
 * Ensure every {@link DEMO_TABLE_FILES} entry is in the platform's table list;
 * return the tables in that order.
 *
 * **Idempotent, and rename-proof.** A table already open for a demo entry is
 * REUSED, never re-added, so repeat ★ clicks cannot grow duplicates. The key is
 * {@link DEMO_ID_TAG} rather than `df.name`: the display name is user-mutable,
 * so a renamed demo table used to reappear as a duplicate on the next click.
 * The `df.name` fallback below keeps a table opened BEFORE this change from
 * duplicating.
 *
 * **Standard AppData path (user-directed 2026-07-17).** Loading routes through
 * {@link loadSampleByName} → `System:AppData/Nca/<file>` rather than a second
 * loader, which also keeps the BUG-01 Float64 re-parse on the demo path: a demo
 * table must carry the same numeric precision as any other bundled sample, or
 * the first numbers a new user ever sees are the wrong ones.
 *
 * **Semantic typing is explicit, and must be.** `grok.shell.addTable` performs
 * NO implicit semantic-type detection — verified empirically against a live
 * server (2026-07-28): the `Smiles` column's `semType` is unset both before and
 * after the `addTable` call. Without the `detectSemanticTypes` call below the
 * structure column stays plain text, and `nca-add-to-study.ts:206` — which reads
 * the LIVE `semType` once and persists it as `isStructure`, never re-deriving it
 * — would latch `false` permanently, so the demo's real structures would never
 * render in a study grid. This is the same call the Admetica demo loader makes
 * (`admetica-app.ts:35,121`) for the same reason.
 *
 * Three branches, and only one of them can prove the frame is already typed:
 *
 * - **Reused via {@link DEMO_ID_TAG}** — already typed: this code loaded and
 *   typed it earlier this session, and only this code writes that tag. No
 *   detection needed.
 * - **Adopted via the `df.name` fallback** — typed state UNKNOWABLE: the branch
 *   exists precisely for a frame this code did *not* load. Detect, and backfill
 *   the tag.
 * - **Freshly loaded** — not typed. Detect.
 */
export async function ensureDemoTables(): Promise<DG.DataFrame[]> {
  const tables: DG.DataFrame[] = [];
  for (const entry of DEMO_TABLE_FILES) {
    // Branch 1 — reused by stable id. Survives a user rename.
    const tagged = grok.shell.tables.find(
      (t: DG.DataFrame) => readFrameTag(t, DEMO_ID_TAG) === entry.file);
    if (tagged !== undefined) {
      tables.push(tagged);
      continue;
    }

    // Branch 2 — adopted by display name: a demo table opened before the tag
    // existed. Stamp the tags FIRST so it graduates to branch 1 on the next
    // click; without the backfill it would re-enter this branch every time and
    // re-run the whole-frame scan below, which a table-count assertion cannot
    // see (the count stays correct while the work repeats).
    const named = grok.shell.tables.find((t: DG.DataFrame) => t.name === entry.name);
    if (named !== undefined) {
      named.setTag(DEMO_ID_TAG, entry.file);
      named.setTag(PROVENANCE_TAG, entry.provenance);
      // eslint-disable-next-line no-await-in-loop
      await grok.data.detectSemanticTypes(named);
      tables.push(named);
      continue;
    }

    // Branch 3 — fresh load. Sequential by design: two small files, and each
    // iteration's already-open check must see the previous iteration's addTable.
    // eslint-disable-next-line no-await-in-loop
    const df = await loadSampleByName(entry.file);
    df.name = entry.name;
    df.setTag(DEMO_ID_TAG, entry.file);
    df.setTag(PROVENANCE_TAG, entry.provenance);
    // eslint-disable-next-line no-await-in-loop
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTable(df);
    tables.push(df);
  }
  return tables;
}

/**
 * Pick a view name. A non-empty request wins verbatim (no dedup — platform
 * handles duplicate view labels gracefully). Callers pass the TABLE name (the Hub
 * form and the study "+ add analysis" dialog have no name input since
 * hub-ungrouped-polish, 2026-10-06) or a sample's view name. Empty / whitespace /
 * null → first-free in the auto-increment series `'New NCA'`, `'New NCA (2)'`, ...
 *
 * Algorithm matches spec §4 decision 4 and §8.2.
 */
export function resolveAnalysisName(
  requested: string | null,
  taken: ReadonlySet<string>,
): string {
  const trimmed = requested?.trim() ?? '';
  if (trimmed.length > 0)
    return trimmed;
  const base = 'New NCA';
  if (!taken.has(base))
    return base;
  for (let i = 2; ; i++) {
    const candidate = `${base} (${i})`;
    if (!taken.has(candidate))
      return candidate;
  }
}

/**
 * Mount the analysis surface (mapping pane + recipe replay) onto an
 * ALREADY-CREATED table view — the shared seam for BOTH the click-open path
 * (`launchAnalysis` → here) and the cold-deep-link hydration path
 * (`nca-analysis-open.ts` → here) (AD-4 / AC-D1). Deliberately does NOT create a
 * view or touch the toolbox — that is the caller's job; calling
 * `grok.shell.addTableView` here would spawn an orphaned double-view inside the
 * router (architect Finding 1).
 *
 * The NA-marker normalization (`normalizeNumericNaMarkers`) stays at each CALL
 * SITE — it must run before the table view's grid binds so the numeric pickers
 * surface (the documented ordering invariant); the helper is the shared
 * mapping-mount + recipe seam, the single code path both entry points route
 * through so byte-identical state is a guarantee, not coincidence.
 *
 * @returns the {@link MappingPaneHandle} so the hydration path can replay a
 *   recipe (`opts.recipe`) and `triggerCompute()` without a user click.
 */
export function mountAnalysisIntoView(
  view: DG.TableView,
  df: DG.DataFrame,
  opts: MountMappingPaneOptions = {},
): MappingPaneHandle {
  return mountMappingPane(view, df, opts);
}

/**
 * Open `df` as a table view, name it, and mount the mapping pane on the left.
 *
 * @returns the created `DG.TableView` so callers can attach view-side state.
 */
export function launchAnalysis(
  df: DG.DataFrame,
  requestedName: string | null,
  source?: SourcePointer,
  targetStudyId?: string,
): DG.TableView {
  // Retype NONMEM-style numeric columns carrying `.`/empty NA markers from
  // STRING to FLOAT in place, so DV/AMT/RATE surface in the mapping-pane's
  // numeric pickers (Ed #15). Runs before the grid + mapping pane are built so
  // both reflect the normalized types. See adapter/na-normalize.ts.
  normalizeNumericNaMarkers(df);
  const taken = new Set<string>(
    Array.from(grok.shell.tableViews, (v: DG.TableView) => v.name),
  );
  const name = resolveAnalysisName(requestedName, taken);
  // The toolbox stays closed: the mapping pane is the analysis's left-hand surface, and the side
  // panel keeps the Browse tree for switching between analyses (shell-panels.ts).
  keepBrowseInFront();
  const view = grok.shell.addTableView(df);
  view.name = name;
  mountAnalysisIntoView(view, df, {
    ...(source !== undefined ? {source} : {}),
    ...(targetStudyId !== undefined ? {targetStudyId} : {}),
  });
  return view;
}
