/**
 * Hub "Studies" list (Slice C): one row per study (label, updated, analyses summary) + per-row
 * link / share / rename / delete. Delete is final — there is no Deleted section or restore.
 *
 * The server filters what the caller may see; per-row actions follow its `~can_*`.
 * Reuse: the Hub card/section idiom (`nca-hub-*` classes from
 * `nca-app-view.ts`) + `ui.table` rows.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  getNcaRepository, uniqueLabel, STANDALONE_ROUTE_SEGMENT, formatStudyAnalysesSummary,
} from '../adapter';
import type {StudyHeader} from '../adapter';
import {renderStudyOverview} from './nca-study-overview';
import {ncaDeepLinkPath} from './nca-studies-router';
import {richTip} from './rich-tooltip';
import {openShareDialog} from './nca-share-dialog';

/** The labels already in use by the caller's visible studies, optionally excluding one study id
 *  (the one being renamed). `analyses` is the standalone route segment — a reserved study label. */
async function takenLabels(excludeId?: string): Promise<Set<string>> {
  const studies = (await getNcaRepository().listStudies()).filter((s) => s.kind === 'study');
  const taken = new Set(studies.filter((s) => s.id !== excludeId).map((s) => s.label));
  taken.add(STANDALONE_ROUTE_SEGMENT);
  return taken;
}

/** The Hub accent blue for icons + links (user-directed 2026-06-29: `--blue-1`,
 *  literal fallback `#2083D5`), so actions read as distinct from plain text. */
export const HUB_BLUE = 'var(--blue-1, #2083D5)';

/** The platform's standard disabled-control colour (`.d4-disabled` → `--grey-3`, e.g. compute-utils
 *  `pipeline-view.css`): a disabled Hub icon is grey, not a faded {@link HUB_BLUE}
 *  (user-directed 2026-10-06, hub-ungrouped-polish). */
export const HUB_DISABLED_GREY = 'var(--grey-3, #B8BAC0)';

/** The user-facing name of the standalone-analyses section (analyses not in a study). The Hub title,
 *  the Browse-tree group and every message that says where a saved analysis lands read THIS, so they
 *  cannot drift (user-directed 2026-10-06, hub-ungrouped-polish). Display only: identifiers and the
 *  `analyses/<label>` URL segment keep their names. */
export const UNGROUPED_SECTION = 'Ungrouped';

/** Tooltip of the {@link UNGROUPED_SECTION} title: what the section holds and how to leave it. */
export const UNGROUPED_TOOLTIP =
  'Saved analyses not grouped into a study yet. Use the Move to study icon to add one to a study.';

/** The tail of every delete confirmation: delete is final (user decision 2026-10-21). */
export const CANNOT_UNDO = 'This cannot be undone.';

/** A Hub icon in {@link HUB_BLUE}. */
export function blueIcon(name: string, handler: () => void, tooltip: string): HTMLElement {
  const el = ui.iconFA(name, handler, tooltip);
  el.style.color = HUB_BLUE;
  el.classList.add('nca-hub-icon');
  return el;
}

/** A Hub icon shown disabled — the platform grey, a not-allowed cursor and a no-op click, with a
 *  tooltip that says why (the Move-to-study precedent, hub-ungrouped-polish). */
export function disabledIcon(name: string, tooltip: string): HTMLElement {
  const el = blueIcon(name, () => {}, tooltip);
  el.style.color = HUB_DISABLED_GREY;
  el.style.cursor = 'not-allowed';
  return el;
}

/** Tooltip of an action the caller's access does not allow (SD-13). */
export const NO_EDIT_TIP = 'View only — ask the study owner for Edit access.';

/** Per-host generation of the Hub's ASYNC list renders (studies, ungrouped analyses, deleted). */
const hostRenderGen = new WeakMap<HTMLElement, number>();

/**
 * Start an async render of `host`. The returned predicate stays true only while this is the NEWEST
 * render of that host. A list render reads storage (await) and must clear + fill its host only
 * afterwards, and only if still current: the Hub refreshes on every `NCA_STUDIES_CHANGED_EVENT`
 * and a list also re-renders itself after a delete / rename, so renders overlap — and a
 * clear-then-await-then-append render let BOTH clear before either appended, leaving the table
 * twice (user QA 2026-10-06). The newest render wins; a superseded one writes nothing.
 */
export function beginHostRender(host: HTMLElement): () => boolean {
  const gen = (hostRenderGen.get(host) ?? 0) + 1;
  hostRenderGen.set(host, gen);
  return () => hostRenderGen.get(host) === gen;
}

/** Render the studies list into `host`, replacing its contents. Re-callable (change event / after
 *  create / rename / delete) — re-reads storage each time. Overlapping calls on one host are safe:
 *  only the newest writes the DOM ({@link beginHostRender}). */
export async function renderStudiesListInto(host: HTMLElement): Promise<void> {
  const isCurrent = beginHostRender(host);
  let studies: StudyHeader[];
  try {
    studies = (await getNcaRepository().listStudies()).filter((s) => s.kind === 'study');
  } catch (e) {
    if (!isCurrent()) return;
    ui.empty(host);
    host.appendChild(ui.divText(`Could not list studies: ${e instanceof Error ? e.message : String(e)}`,
      {classes: 'text-muted-foreground'}));
    return;
  }
  // Clear only AFTER the read and only if still the newest render (see beginHostRender).
  if (!isCurrent()) return;
  ui.empty(host);
  if (studies.length === 0) {
    const empty = ui.divText('No studies yet.');
    empty.style.marginLeft = '2px';
    empty.style.marginTop = '12px';
    empty.style.color = 'var(--grey-4, #9497A0)';
    host.appendChild(empty);
    return;
  }
  host.appendChild(studiesTable(studies, host));
}

/** A `td` with the standard cell padding. Shared with the study-overview table so
 *  the row padding stays a single knob. */
export function td(child: HTMLElement | string): HTMLTableCellElement {
  const cell = document.createElement('td');
  cell.style.padding = '8px 28px 2px 0'; // bottom 2 (tight row-to-row gap; header→rows unchanged)
  cell.style.verticalAlign = 'middle';
  cell.append(typeof child === 'string' ? document.createTextNode(child) : child);
  return cell;
}

/** A Hub table with `ui.h3` headers (the Studies / Ungrouped idiom). */
export function hubTable(headers: string[], rows: HTMLTableRowElement[]): HTMLElement {
  const table = document.createElement('table');
  table.className = 'nca-studies-table';
  table.style.marginTop = '12px';
  table.style.width = 'max-content'; // pack columns to content (parity with the study view)
  const thead = document.createElement('thead');
  const htr = document.createElement('tr');
  for (const h of headers) {
    const th = document.createElement('th');
    th.style.textAlign = 'left';
    th.style.padding = '0 28px 2px 0';
    const h3 = ui.h3(h);
    h3.style.margin = '0'; // h3 default block margins otherwise inflate the header→rows gap
    th.appendChild(h3);
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);
  const tbody = document.createElement('tbody');
  for (const r of rows) tbody.appendChild(r);
  table.appendChild(tbody);
  return table;
}

function studiesTable(studies: StudyHeader[], host: HTMLElement): HTMLElement {
  const refresh = (): void => {void renderStudiesListInto(host);};
  return hubTable(['Name', 'Updated', 'Analyses', ''], studies.map((s) => studyRow(s, refresh)));
}

function studyRow(p: StudyHeader, refresh: () => void): HTMLTableRowElement {
  const tr = document.createElement('tr');
  // Name = a clickable link that opens the study overview in a new view. The
  // displayed text is capped at 40 chars (…); the full name lives in the tooltip. 40, not 20
  // (user-directed 2026-11-04): 20 cut every demo study label (`Demo · rat dose escalation`).
  const MAX = 40;
  const display = p.label.length > MAX ? `${p.label.slice(0, MAX)}…` : p.label;
  const link = ui.link(display, () => {void openStudyOverview(p.id);});
  link.style.color = HUB_BLUE;
  ui.tooltip.bind(link, richTip({title: p.label, summary: 'Open study in a new view'}));
  tr.appendChild(td(link));
  tr.appendChild(td(shortDate(p.updatedAt)));
  // D19 / F-S6.3: `3 analyses · 1 in review` — the ONE summary wording (workflow-stage.ts).
  tr.appendChild(td(formatStudyAnalysesSummary(p.analyses.map((a) => a.stageCache),
    p.analyses.reduce((n, a) => n + a.pipelineFlagCount, 0))));
  tr.appendChild(td(studyActions(p, refresh)));
  return tr;
}

/** Copy link · Share · Rename · Delete — each per the caller's access (SD-13). Shared with the
 *  Browse tree's context menu through {@link renameFlow} / {@link confirmDelete}. */
function studyActions(p: StudyHeader, refresh: () => void): HTMLElement {
  const icons: HTMLElement[] = [
    blueIcon('share', () => copyStudyLink(p.label), 'Copy a shareable link to this study'),
    p.access.canShare ?
      blueIcon('user-plus', () => {void openShareDialog(p, refresh);}, 'Share this study with a group (View)') :
      disabledIcon('user-plus', 'Only the study owner can share it.'),
    p.access.canEdit ?
      blueIcon('pen', () => {void renameFlow(p, refresh);}, 'Rename study') :
      disabledIcon('pen', NO_EDIT_TIP),
    p.access.canDelete ?
      blueIcon('trash-alt', () => confirmDelete(p, refresh), 'Delete study') :
      disabledIcon('trash-alt', NO_EDIT_TIP),
  ];
  const actions = ui.divH(icons);
  actions.style.gap = '12px'; // +4px over the default 8px (user-directed 2026-06-29)
  return actions;
}

/** Navigate to a study's overview view. `mount: 'view'` (default) adds it as a permanent tab;
 *  `mount: 'preview'` mounts it in the platform preview slot (auto-replaces any prior preview) —
 *  the browse tree's single-click peek. The overview computes nothing (cache-only). */
export async function openStudyOverview(studyId: string, mount: 'view' | 'preview' = 'view'): Promise<DG.ViewBase> {
  const view = await renderStudyOverview(studyId);
  if (mount === 'preview') {
    // Drop the current preview slot occupant before showing ours (DiffStudio releaseSharedPreview:
    // grok.shell.preview = null, NEVER view.close() — closing the app-landing view closes the app).
    grok.shell.preview = null;
    grok.shell.addPreview(view as DG.ViewBase);
  } else
    grok.shell.addView(view as DG.ViewBase);
  return view as DG.ViewBase;
}

/** Create-study flow: prompt for a label, create a creator-only empty study (the platform grants
 *  the author everything — `promotion: eager`), refresh the list. */
export async function createStudyFlow(host: HTMLElement): Promise<void> {
  const labelInput = ui.input.string('Study name', {value: 'New study'});
  ui.dialog('New NCA study')
    .add(labelInput.root)
    .onOK(() => {
      void (async () => {
        try {
          const raw = (labelInput.value ?? '').trim() || 'New study';
          const s = await getNcaRepository().createStudy(uniqueLabel(raw, await takenLabels()));
          await renderStudiesListInto(host);
          grok.shell.info(`Created study '${s.label}'.`);
        } catch (e) {
          grok.shell.error(`Could not create the study: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}

/** Rename = label-only; touches no id. Re-applies the `label-<n>` uniqueness gate (excluding the
 *  study's own current label) so a rename can't collide with a sibling and break label-keyed
 *  routing (§6). */
export async function renameFlow(project: {id: string; label: string}, onDone: () => void): Promise<void> {
  const labelInput = ui.input.string('Study name', {value: project.label});
  ui.dialog('Rename study')
    .add(labelInput.root)
    .onOK(() => {
      void (async () => {
        try {
          const raw = (labelInput.value ?? '').trim() || project.label;
          await getNcaRepository().renameStudy(project.id, uniqueLabel(raw, await takenLabels(project.id)));
          onDone();
        } catch (e) {
          grok.shell.error(`Could not rename the study: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}

/** Final delete behind a blocking confirm (user decision 2026-10-21: no Deleted section, no
 *  restore): "Delete "X" and its N analyses? This cannot be undone." */
export function confirmDelete(project: {id: string; label: string; analysisCount: number}, onDone: () => void): void {
  const n = project.analysisCount;
  const name = document.createElement('b');
  name.textContent = project.label;
  const body = document.createElement('div');
  body.style.maxWidth = '420px';
  body.style.whiteSpace = 'normal';
  body.append('Delete "', name, n === 0 ? '"?' : `" and its ${n} ${n === 1 ? 'analysis' : 'analyses'}?`,
    ` ${CANNOT_UNDO}`);
  ui.dialog('Delete study')
    .add(body)
    .onOK(() => {
      void (async () => {
        try {
          await getNcaRepository().deleteStudy(project.id);
          grok.shell.info(`Deleted study '${project.label}'.`);
          onDone();
        } catch (e) {
          grok.shell.error(`Could not delete the study: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}

/** Trim an ISO timestamp to the date portion for the terse list meta. Shared with
 *  the study-overview table. */
export function shortDate(iso: string): string {
  const t = iso.indexOf('T');
  return t > 0 ? iso.slice(0, t) : iso;
}

/** Copy a shareable deep-link to the study to the clipboard. The URL mirrors the
 *  study's address-bar URL (`/apps/Nca/<study label>`), so pasting it opens the study
 *  overview directly; the label is `encodeURIComponent`-encoded to match the
 *  router's decode. Falls back to surfacing the URL in a toast when the clipboard
 *  API is unavailable (non-secure context). */
export function copyStudyLink(studyLabel: string): void {
  const url = `${window.location.origin}${ncaDeepLinkPath(studyLabel)}`;
  const clip = window.navigator?.clipboard;
  if (clip?.writeText) {
    clip.writeText(url).then(
      () => grok.shell.info('Study link copied to clipboard.'),
      () => grok.shell.info(url));
  } else
    grok.shell.info(url);
}
