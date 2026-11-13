/**
 * Hub "Ungrouped" section — standalone analyses, not (yet) in a study.
 * Modelled on the Studies table (`nca-studies-list.ts`):
 * same `ui.h3` header row, same row padding (`td`), same blue action icons.
 *
 * Columns: Name (link — opens the analysis where it stopped) / Updated / Stage / actions:
 * Move to study… · Copy link · Rename · Delete (user-directed 2026-09-30). The Stage badge is
 * `formatStageBadge` of the cached stage with the engine-staleness overlay — the ONE wording (AD-9),
 * in plain text. Move to study is offered only when the analysis has results (`stage ≠ setup`, P7),
 * and warns first when profiles are still to review (D14).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {formatStageBadge, STANDALONE_ROUTE_SEGMENT, uniqueLabel, getNcaRepository, isStale} from '../adapter';
import type {StudyHeader, AnalysisHeader, StageCache} from '../adapter';
import {
  td, blueIcon, HUB_BLUE, HUB_DISABLED_GREY, shortDate, beginHostRender, hubTable, CANNOT_UNDO,
} from './nca-studies-list';
import {openAnalysis} from './nca-analysis-open';
import {moveStandaloneIntoStudy} from './nca-add-to-study';
import {ncaStandaloneLinkPath} from './nca-studies-router';
import {richTip} from './rich-tooltip';

/** Tooltip of the disabled Move to study icon: why it is off AND how to turn it on. */
const MOVE_DISABLED_TIP =
  'Compute the analysis first: open it and click Confirm. Only computed analyses can be moved to a study.';

/** A standalone container and its one analysis. */
type Standalone = {container: StudyHeader; entry: AnalysisHeader};

/** The caller's standalone analyses (one per container). */
async function listStandalone(): Promise<Standalone[]> {
  return (await getNcaRepository().listStudies())
    .filter((s) => s.kind === 'standalone' && s.analyses.length > 0)
    .map((s) => ({container: s, entry: s.analyses[0]}));
}

/** Render the Ungrouped list into `host`, replacing its contents. Re-callable. Returns how many
 *  analyses were listed — the hub hides the whole section (title included) when there are none
 *  (user-directed 2026-09-29). A listing error still renders + counts 1, so the failure stays
 *  visible. `onCount` receives the same number on EVERY render, including the list's own
 *  re-renders after a delete / move. Overlapping renders of one host are safe: only the newest
 *  writes the DOM and reports a count (`beginHostRender`). */
export async function renderAnalysesListInto(host: HTMLElement, onCount?: (n: number) => void): Promise<number> {
  const isCurrent = beginHostRender(host);
  let items: Standalone[];
  try {
    items = await listStandalone();
  } catch (e) {
    if (!isCurrent()) return 0;
    ui.empty(host);
    host.appendChild(ui.divText(`Could not list analyses: ${e instanceof Error ? e.message : String(e)}`));
    onCount?.(1);
    return 1;
  }
  // Clear only AFTER the read and only if still the newest render (BUG-14).
  if (!isCurrent()) return 0;
  ui.empty(host);
  onCount?.(items.length);
  if (items.length === 0) return 0;
  const refresh = (): void => {void renderAnalysesListInto(host, onCount);};
  host.appendChild(hubTable(['Name', 'Updated', 'Stage', ''], items.map((it) => analysisRow(it, refresh))));
  return items.length;
}

/** The Stage badge element: the ONE wording, plain text (no colour — user-directed 2026-09-30). */
export function stageBadge(stageCache: StageCache | null, stale: boolean): HTMLElement {
  return ui.divText(formatStageBadge(stageCache, stale));
}

/** Can this standalone analysis join a study (P7 / AC-S1.7.2: it has results)? */
export function canMove(entry: AnalysisHeader): boolean {
  return entry.stageCache !== null && entry.stageCache.stage !== 'setup' && entry.hasCache;
}

function analysisRow(it: Standalone, refresh: () => void): HTMLTableRowElement {
  const {container, entry} = it;
  const tr = document.createElement('tr');
  const MAX = 24;
  const display = entry.label.length > MAX ? `${entry.label.slice(0, MAX)}…` : entry.label;
  const link = ui.link(display, () => {void openAnalysis(entry.id);});
  link.style.color = HUB_BLUE;
  ui.tooltip.bind(link, richTip({title: entry.label, summary: 'Open the analysis where it stopped'}));
  // Move to study — FIRST of the icons, as `layer-plus`, and ALWAYS shown (user-directed 2026-09-30).
  // Before results exist it is shown DISABLED — the platform disabled grey, a not-allowed cursor, a
  // no-op click — with a tooltip that says why and how to enable it (user-directed 2026-10-06).
  const movable = canMove(entry);
  const moveIcon = blueIcon('layer-plus', () => {
    if (movable) void moveToStudyFlow(it, refresh);
  }, movable ? 'Move to study…' : MOVE_DISABLED_TIP);
  if (!movable) {
    moveIcon.style.color = HUB_DISABLED_GREY;
    moveIcon.style.cursor = 'not-allowed';
  }
  const actionRow = ui.divH([
    moveIcon,
    blueIcon('share', () => copyAnalysisLink(entry.label), 'Copy a shareable link to this analysis'),
    blueIcon('pen', () => renameAnalysisFlow(it, refresh), 'Rename analysis'),
    blueIcon('trash-alt', () => confirmDeleteAnalysis(it, refresh), 'Delete analysis'),
  ]);
  actionRow.style.gap = '12px';
  tr.appendChild(td(link));
  tr.appendChild(td(shortDate(entry.updatedAt || container.updatedAt)));
  tr.appendChild(td(stageBadge(entry.stageCache, isStale(entry))));
  tr.appendChild(td(actionRow));
  return tr;
}

/** Copy the standalone deep link (`/apps/Nca/analyses/<label>`) — the SAME builder the router
 *  resolves, so a copied link survives a later move into a study (the moved-from redirect). */
export function copyAnalysisLink(label: string): void {
  const url = `${window.location.origin}${ncaStandaloneLinkPath(label)}`;
  void navigator.clipboard.writeText(url).then(
    () => grok.shell.info('Link copied'),
    () => grok.shell.warning(`Could not copy the link: ${url}`));
}

/** Rename a standalone analysis (its container's label follows — one transaction). Unique among
 *  the caller's standalone analyses so `analyses/<label>` stays unambiguous. */
export function renameAnalysisFlow(it: Standalone, onDone: () => void): void {
  const input = ui.input.string('Analysis name', {value: it.entry.label});
  ui.dialog('Rename analysis').add(input.root).onOK(() => {
    void (async () => {
      try {
        const taken = new Set((await listStandalone()).filter((x) => x.container.id !== it.container.id)
          .map((x) => x.entry.label));
        const label = uniqueLabel((input.value ?? '').trim() || it.entry.label, taken);
        await getNcaRepository().renameAnalysis(it.entry.id, label);
        onDone();
      } catch (e) {
        grok.shell.error(`Could not rename: ${e instanceof Error ? e.message : String(e)}`);
      }
    })();
  }).show();
}

/** Final delete behind a blocking confirm (user decision 2026-10-21: no restore). */
export function confirmDeleteAnalysis(it: Standalone, onDone: () => void): void {
  const name = document.createElement('b');
  name.textContent = it.entry.label;
  const body = document.createElement('div');
  body.style.maxWidth = '420px';
  body.style.whiteSpace = 'normal';
  body.append('Delete "', name, `"? ${CANNOT_UNDO}`);
  ui.dialog('Delete analysis').add(body).onOK(() => {
    void (async () => {
      try {
        await getNcaRepository().deleteStudy(it.container.id);
        grok.shell.info(`Deleted analysis '${it.entry.label}'.`);
        onDone();
      } catch (e) {
        grok.shell.error(`Could not delete: ${e instanceof Error ? e.message : String(e)}`);
      }
    })();
  }).show();
}

/**
 * D14 soft gate (F-S6.5): adding an analysis with profiles still to review warns first —
 * `Review first` (open it on Review λz) or `Add anyway` (the study row then carries
 * `N not reviewed`). Resolves true when the caller should proceed.
 */
export function confirmNotReviewed(pending: number, onReviewFirst: () => void): Promise<boolean> {
  if (pending <= 0) return Promise.resolve(true);
  return new Promise((resolve) => {
    let decided = false;
    const dlg = ui.dialog('Profiles not reviewed')
      .add(ui.divText(`${pending} ${pending === 1 ? 'profile is' : 'profiles are'} not reviewed. The study ` +
        `grid will show "${pending} not reviewed" for this analysis until they are.`,
      {style: {maxWidth: '420px', whiteSpace: 'normal'}}))
      .addButton('Review first', () => {decided = true; dlg.close(); onReviewFirst(); resolve(false);})
      .onOK(() => {decided = true; resolve(true);});
    dlg.onClose.subscribe(() => {if (!decided) resolve(false);});
    // The OK button reads "Add anyway".
    dlg.show();
    const ok = dlg.root.querySelector('button[name="button-OK"]') as HTMLElement | null;
    if (ok !== null) ok.textContent = 'Add anyway';
  });
}

/** Hub "Move to study…" (D5 / F-S2.5): pick a study (or a new one), pass the D14 warning, then the
 *  ONE move implementation. Only studies the caller can edit are offered (SD-13). */
export async function moveToStudyFlow(it: Standalone, onDone: () => void): Promise<void> {
  const {container, entry} = it;
  const proceed = await confirmNotReviewed(entry.stageCache?.stage === 'review' ? entry.stageCache.pending : 0,
    () => {void openAnalysis(entry.id);});
  if (!proceed) return;
  const all = (await getNcaRepository().listStudies()).filter((s) => s.kind === 'study');
  const editable = all.filter((s) => s.access.canEdit);
  const NEW = '(new study)';
  const byLabel = new Map(editable.map((s) => [s.label, s.id]));
  const choice = ui.input.choice<string>('Study', {value: NEW, items: [NEW, ...byLabel.keys()]});
  const newName = ui.input.string('New study name', {value: 'New study'});
  const sync = (): void => {newName.root.style.display = choice.value === NEW ? '' : 'none';};
  choice.onChanged.subscribe(sync);
  sync();
  ui.dialog('Move to study').add(choice.root).add(newName.root).onOK(() => {
    void (async () => {
      const moved = await moveStandaloneIntoStudy(container.id, choice.value === NEW ?
        {kind: 'new', label: uniqueLabel((newName.value ?? '').trim() || 'New study',
          new Set([...all.map((s) => s.label), STANDALONE_ROUTE_SEGMENT]))} :
        {kind: 'existing', studyId: byLabel.get(choice.value as string)!});
      if (moved !== null) onDone();
    })();
  }).show();
}
