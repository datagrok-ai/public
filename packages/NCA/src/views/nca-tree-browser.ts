/**
 * NCA Browse-tree builder — populates the app's tree node with one **group per
 * study** and one **leaf per analysis** (two-level hierarchy). A **single-click
 * previews** the study/analysis in the platform preview slot (which auto-replaces
 * the prior preview); a **double-click** or the context-menu **Run** opens the full
 * detail view; **↑/↓/←/→** navigate the tree (platform built-in + a focus-guard so
 * a mounted preview doesn't steal focus). The tree rebuilds reactively on study
 * mutation.
 *
 * Hybrid of the two reference apps:
 *   - **Hit Design** mechanic — `getOrCreateGroup(studyLabel)` nests analyses
 *     under their study; the platform pairs this builder to the app via the
 *     `//meta.app` binding on `ncaAppTreeBrowser` (package.ts).
 *   - **MPO Profiles** mechanics — load from a manager (`NcaRepository.listStudies`),
 *     open a custom detail view (a study group reuses the Hub's `openStudyOverview`;
 *     an analysis leaf uses the cold `openAnalysisByDeepLink`), both taking a
 *     `mount: 'view' | 'preview'` so the same builders serve the single-click preview
 *     and the double-click/Run full open, and stay live via a module-scoped resubscribe
 *     to `NCA_STUDIES_CHANGED_EVENT` (Hit Design builds once; we don't).
 *
 * Captions use the study/analysis `label`. Labels are unique by construction (the
 * create/rename `label-<n>` gate); `legacyDisambiguate` is a defensive caption
 * fallback for pre-gate `study.json` files only, so `getOrCreateGroup` can never
 * merge two distinct studies.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {Subscription} from 'rxjs';

import {NCA_STUDIES_CHANGED_EVENT, uniqueLabel, getNcaRepository} from '../adapter';
import type {StudyHeader, AnalysisHeader} from '../adapter';
import {openAnalysis} from './nca-analysis-open';
import {NCA_RESULTS_PRESENTED_EVENT} from './shell-panels';
import {
  openStudyOverview, renameFlow, confirmDelete, copyStudyLink, UNGROUPED_SECTION,
} from './nca-studies-list';
import {renameAnalysisFlow, confirmDeleteAnalysis, copyAnalysisLink} from './nca-study-overview';
import {
  copyAnalysisLink as copyStandaloneLink, renameAnalysisFlow as renameStandaloneFlow,
  confirmDeleteAnalysis as confirmDeleteStandalone, moveToStudyFlow, canMove,
} from './nca-analyses-list';

/** Module-scoped so each tree (re)build replaces — not stacks — the change-event
 *  subscription (MPO precedent: unsubscribe-then-resubscribe avoids leaks across
 *  re-entries when the node is collapsed/expanded repeatedly). */
let ncaTreeSub: Subscription | null = null;

/** Single root-level selection subscription that opens the debounced PREVIEW. It
 *  REPLACES the old per-node `onSelected` full-opens so the platform's built-in
 *  ↑/↓/←/→ tree navigation is not broken by a focus-stealing preview fired on every
 *  node (DiffStudio `setupTreeItemSelection`, app.ts). Unsubscribe-then-resubscribe
 *  per (re)build, exactly like `ncaTreeSub`. */
let ncaTreeSelSub: Subscription | null = null;
/** One-shot focus restore for the next view change (the preview mount). */
let ncaTreeFocusSub: Subscription | null = null;
/** Teardown for the short post-mount focus-guard window (see `guardTreeFocus`). */
let ncaTreeFocusGuardCleanup: (() => void) | null = null;
/** Pending debounced-preview timer id; a following double-click cancels it. */
let pendingPreviewTimer: number | null = null;
/** The analysis view the tree last PREVIEWED, and whether the analyst is still walking the tree
 *  (a selection made it true; a mousedown outside the tree ends it). An analysis preview computes
 *  after it mounts — on a large table well past `guardTreeFocus`'s 1.5 s — and the docked results
 *  grid then takes keyboard focus, so ↑/↓ stopped walking the tree (user-reported 2026-10-12).
 *  `NCA_RESULTS_PRESENTED_EVENT` for THIS preview hands focus back while the walk is still on. */
let lastAnalysisPreview: DG.ViewBase | null = null;
let treeWalkActive = false;
let ncaResultsPresentedSub: Subscription | null = null;
/** The Browse tree's root element, as of the latest (re)build — the walk's "inside" test. */
let ncaTreeRoot: HTMLElement | null = null;

/** Debounce before a selection opens its preview. Long enough that (a) a double-click
 *  cancels it and opens the full view without racing a preview, and (b) an arrow-key
 *  sweep does not preview every node it passes over. Matches DiffStudio's DBL_CLICK_DELAY. */
const PREVIEW_DEBOUNCE_MS = 250;

/** What a tree node points at — read back by the centralized selection handler to open
 *  the matching preview. Stashed on `node.value` (an object, per the `getOrCreateGroup` /
 *  `item` value param). A `study` node previews its overview; an `analysis` node previews
 *  its (computing) analysis view. */
export type NcaNodeMeta =
  | {kind: 'study'; studyId: string}
  | {kind: 'analysis'; analysisId: string};

/** A selected node's value as OUR meta, or null. Tests the kind and the id, never mere presence:
 *  the ⚠ marker items carry null, but the APP node carries the platform's own value (the app
 *  function), which the handler once read as an analysis and opened `undefined` (QA 2026-10-21). */
export function asNcaNodeMeta(value: unknown): NcaNodeMeta | null {
  if (value === null || typeof value !== 'object') return null;
  const v = value as Record<string, unknown>;
  if (v.kind === 'study' && typeof v.studyId === 'string') return {kind: 'study', studyId: v.studyId};
  if (v.kind === 'analysis' && typeof v.analysisId === 'string') return {kind: 'analysis', analysisId: v.analysisId};
  return null;
}

/** Cancel any preview scheduled by the centralized selection handler. */
function cancelPendingPreview(): void {
  if (pendingPreviewTimer !== null) {
    clearTimeout(pendingPreviewTimer);
    pendingPreviewTimer = null;
  }
}

/** Drop whatever is in the platform's single preview slot. `grok.shell.preview = null` is the
 *  graceful slot-clear DiffStudio uses (`releaseSharedPreview`, app.ts:631) — it removes the
 *  transient preview WITHOUT `view.close()`. That distinction is load-bearing: `.close()` on
 *  the app-landing view closes the whole app (user QA 2026-09-21 — an earlier `.close()`
 *  attempt made the app vanish on click and cascaded into a bogus "study unavailable"),
 *  whereas clearing the slot removes only the transient preview and leaves permanent views
 *  (including the app) alone. The `mount: 'preview'` openers clear the slot right before their
 *  own `addPreview`, so a new peek REPLACES the old preview (incl. the app-landing preview);
 *  this helper is for the FULL-open path (double-click / Run), which drops the transient
 *  preview before adding the permanent view (mirrors DiffStudio's Run: release → `addView`). */
function clearNcaPreviewSlot(): void {
  try {
    grok.shell.preview = null;
  } catch {/* no active preview slot */}
}

/** Double-click a node's caption → cancel any pending single-click preview and open the
 *  FULL view (`addView`). This is where the old single-click-opens behavior now lives. */
function wireOpenOnDoubleClick(node: DG.TreeViewNode, open: () => Promise<unknown>): void {
  node.captionLabel.addEventListener('dblclick', (ev) => {
    ev.preventDefault();
    ev.stopPropagation();
    cancelPendingPreview();
    clearNcaPreviewSlot(); // drop the transient preview before the full view opens
    void open();
  });
}

/** Capture the browse-panel scroll position and arm a one-shot focus restore for the next
 *  view change (the preview mount), so focus returns to the tree root and the platform's
 *  arrow-key navigation keeps targeting the tree. Ported from DiffStudio
 *  `armTreeFocusRestore`/`restoreTreeFocus` (the Tutorials Demo-app `_guardTreeFocus`
 *  pattern; TODO.md:1448 — nav itself is platform-built-in, NOT a custom keydown handler). */
function armTreeFocusRestore(treeNode: DG.TreeViewGroup): void {
  const treeRoot = treeNode.rootNode.root;
  const panelRoot = treeRoot.parentElement;
  if (panelRoot === null)
    return;
  const scrollY = panelRoot.scrollTop;
  ncaTreeFocusSub?.unsubscribe();
  ncaTreeFocusSub = grok.events.onCurrentViewChanged.subscribe(() => {
    ncaTreeFocusSub?.unsubscribe();
    ncaTreeFocusSub = null;
    panelRoot.scrollTo(0, scrollY);
    treeRoot.focus();
    guardTreeFocus(treeNode);
  });
}

/** For ~1.5s after a preview mounts, redirect focus back to the tree if a grid/viewer
 *  steals it (otherwise ↑/↓ would target the preview). Cleared on the first user
 *  mousedown or after the timeout. Ported from DiffStudio `guardTreeFocus`. */
function guardTreeFocus(treeNode: DG.TreeViewGroup): void {
  ncaTreeFocusGuardCleanup?.();
  const treeRoot = treeNode.rootNode.root;
  const onFocusIn = (e: FocusEvent): void => {
    if (!treeRoot.contains(e.target as Node))
      treeRoot.focus();
  };
  const clear = (): void => {
    document.removeEventListener('focusin', onFocusIn, true);
    document.removeEventListener('mousedown', clear, true);
    clearTimeout(timer);
    ncaTreeFocusGuardCleanup = null;
  };
  const timer = window.setTimeout(clear, 1500);
  ncaTreeFocusGuardCleanup = clear;
  document.addEventListener('focusin', onFocusIn, true);
  document.addEventListener('mousedown', clear, true);
}

/** A small inline spinner shown while the study list loads (Hit Design idiom). */
function makeLoader(): HTMLElement {
  const d = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  d.innerHTML = '<div class="grok-loader"><div></div><div></div><div></div><div></div></div>';
  return d;
}

/** Map each visible study id → a unique caption. Labels are unique among one user's studies
 *  (create/rename gate), but sharing can bring a second study of the same name into view (SD-26):
 *  it gets a `label-<n>` caption, keeping `getOrCreateGroup` collision-free (§6). */
function disambiguate(studies: StudyHeader[]): Map<string, string> {
  const out = new Map<string, string>();
  const seen = new Set<string>();
  for (const s of studies) {
    const caption = uniqueLabel(s.label, seen);
    seen.add(caption);
    out.set(s.id, caption);
  }
  return out;
}

/** The reused Hub/overview action handlers take an `onDone`/`refresh` callback to
 *  re-render their own table; the tree instead auto-rebuilds on
 *  `NCA_STUDIES_CHANGED_EVENT` (fired by the coordinator's save/delete inside those
 *  handlers), so here the callback is a no-op — no manual re-render needed. */
const treeAutoRefreshes = (): void => {};

/** Right-click a study group → Run / Copy Link / (sep) / Rename / Delete (reuses the
 *  Hub's study handlers; MPO precedent, datagrok-app-tree-pattern.md §8). "Run" opens
 *  the full overview (`addView`) — the same as a double-click (nca-tree-node-interactions
 *  D3); single-click only previews. */
function wireStudyActions(node: DG.TreeViewNode, study: StudyHeader): void {
  node.captionLabel.addEventListener('contextmenu', (ev) => {
    ev.preventDefault();
    ev.stopPropagation();
    const menu = DG.Menu.popup()
      .item('Run', () => {
        cancelPendingPreview(); // a single-click may have armed a preview timer; cancel it (B3 race)
        clearNcaPreviewSlot();
        void openStudyOverview(study.id, 'view');
      })
      .item('Copy Link', () => copyStudyLink(study.label));
    // SD-13: write actions only where the caller's access allows them.
    if (study.access.canEdit || study.access.canDelete) menu.separator();
    if (study.access.canEdit) menu.item('Rename', () => {void renameFlow(study, treeAutoRefreshes);});
    if (study.access.canDelete) menu.item('Delete', () => confirmDelete(study, treeAutoRefreshes));
    menu.show({x: ev.clientX, y: ev.clientY, causedBy: ev});
  });
}

/** Right-click an analysis leaf → Run / Copy Link / (sep) / Rename / Delete (reuses the
 *  study overview's analysis handlers). "Run" opens (and computes) the full analysis view
 *  (`addView`) — the same as a double-click; single-click only previews. */
function wireAnalysisActions(node: DG.TreeViewNode, study: StudyHeader, analysis: AnalysisHeader): void {
  node.captionLabel.addEventListener('contextmenu', (ev) => {
    ev.preventDefault();
    ev.stopPropagation();
    const menu = DG.Menu.popup()
      .item('Run', () => {
        cancelPendingPreview(); // a single-click may have armed a preview timer; cancel it (B3 race)
        clearNcaPreviewSlot();
        void openAnalysis(analysis.id, 'view');
      })
      .item('Copy Link', () => copyAnalysisLink(study.label, analysis.label));
    if (study.access.canEdit) {
      menu.separator()
        .item('Rename', () => renameAnalysisFlow(study.id, analysis, treeAutoRefreshes))
        .item('Delete', () => confirmDeleteAnalysis(analysis, treeAutoRefreshes));
    }
    menu.show({x: ev.clientX, y: ev.clientY, causedBy: ev});
  });
}

/** Right-click a standalone analysis → Run / Copy Link / (sep) / Rename / Delete / Move to study…
 *  — the SAME handlers as the hub Ungrouped section (tree/hub parity, datagrok-app-tree-pattern.md).
 *  Move is offered only when the analysis has results (stage ≠ setup, P7). */
function wireStandaloneActions(node: DG.TreeViewNode, container: StudyHeader): void {
  const entry = container.analyses[0];
  const it = {container, entry};
  node.captionLabel.addEventListener('contextmenu', (ev) => {
    ev.preventDefault();
    ev.stopPropagation();
    const menu = DG.Menu.popup()
      .item('Run', () => {
        cancelPendingPreview();
        clearNcaPreviewSlot();
        void openAnalysis(entry.id, 'view');
      })
      .item('Copy Link', () => copyStandaloneLink(entry.label))
      .separator()
      .item('Rename', () => renameStandaloneFlow(it, treeAutoRefreshes))
      .item('Delete', () => confirmDeleteStandalone(it, treeAutoRefreshes));
    if (canMove(entry))
      menu.item('Move to study…', () => {void moveToStudyFlow(it, treeAutoRefreshes);});
    menu.show({x: ev.clientX, y: ev.clientY, causedBy: ev});
  });
}

/**
 * Build (and keep live) the NCA Browse-tree under `treeNode`. Called by the
 * platform when the app's tree node is expanded (`ncaAppTreeBrowser`).
 */
export async function buildNcaTree(treeNode: DG.TreeViewGroup): Promise<void> {
  // Rebuilds overlap: a Move fires TWO change events (save into the study, delete from Ungrouped),
  // and every event starts a refresh. So a refresh reads storage FIRST and rebuilds the tree
  // synchronously afterwards, only if no newer refresh started meanwhile. The old
  // clear → await → append order let two refreshes both clear, then both append — duplicated
  // leaves inside the reused study groups and a second Ungrouped group (the same race as the Hub
  // lists, user QA 2026-10-06) — and the second one read the expand state off the cleared tree.
  let buildGen = 0;
  const refresh = async (): Promise<void> => {
    const gen = ++buildGen;
    // A loader only on the first build (empty tree); a rebuild keeps the current tree on screen
    // until the new one is ready. The NEWEST refresh's clear below always removes it, so a
    // superseded refresh leaves it in place (no empty-tree flash while the newer one reads).
    if (treeNode.children.length === 0) treeNode.item(makeLoader());

    let all: StudyHeader[] | null = null;
    let listError: unknown = null;
    try {
      all = await getNcaRepository().listStudies();
    } catch (e) {
      listError = e;
    }
    if (gen !== buildGen) return; // superseded — the newer refresh rebuilds (and drops the loader)

    // Preserve the user's expand/collapse state across a rebuild. A rebuild fires
    // on a structural change (create / rename / delete / add-analysis); without
    // this capture the tree would snap back to a fixed layout and collapse
    // whatever the user had open. Keyed by caption (the group `text`). Captured
    // here — after the reads, from the tree as it stands — and the tree is then
    // cleared and rebuilt with no await in between.
    const wasExpanded = new Set<string>();
    let hadGroups = false;
    for (const child of treeNode.children) {
      if (child instanceof DG.TreeViewGroup) {
        hadGroups = true;
        if (child.expanded) wasExpanded.add(child.text);
      }
    }

    treeNode.items.forEach((it) => it.remove()); // clear existing children incl. the loader (MPO)

    if (all === null) {
      treeNode.item(`⚠ Could not list studies: ${listError instanceof Error ? listError.message : String(listError)}`);
      return;
    }
    const studies = all.filter((s) => s.kind === 'study');
    // nca-workflow-stages F-S6.2 (D4): the standalone analyses, from the same read.
    const standalone = all.filter((s) => s.kind === 'standalone' && s.analyses.length > 0);

    const captionFor = disambiguate(studies);
    let first = true;
    for (const study of studies) {
      // Intermediate level: one group per study (Hit Design getOrCreateGroup).
      // Restore prior expansion; on the very first build, expand the first study.
      const caption = captionFor.get(study.id)!;
      const expanded = hadGroups ? wasExpanded.has(caption) : first;
      // Stash the study id on the node so the centralized selection handler (below) can
      // open its preview. Single-click PREVIEWs; a double-click / "Run" opens the full view.
      const studyMeta: NcaNodeMeta = {kind: 'study', studyId: study.id};
      const group = treeNode.getOrCreateGroup(caption, studyMeta, expanded);
      first = false;
      wireOpenOnDoubleClick(group, () => openStudyOverview(study.id, 'view'));
      wireStudyActions(group, study);

      // Leaf level: one item per analysis. The open reads the analysis fresh BEFORE mounting — a
      // tree leaf can be clicked without first selecting its study group, and the list captured
      // at tree-build time may be stale by then. This is the same self-sufficient path the URL
      // deep-link uses (AC-D1: one shared mount).
      for (const analysis of study.analyses) {
        const analysisMeta: NcaNodeMeta = {kind: 'analysis', analysisId: analysis.id};
        const leaf = group.item(analysis.label, analysisMeta);
        wireOpenOnDoubleClick(leaf, () => openAnalysis(analysis.id, 'view'));
        wireAnalysisActions(leaf, study, analysis);
      }
    }
    // nca-workflow-stages F-S6.2 (D4): standalone analyses as their own group, after the studies,
    // captioned with the Hub section's name (UNGROUPED_SECTION — one name for one entity,
    // user-directed 2026-10-06). `group()` (not getOrCreateGroup) so a study that happens to carry
    // the same caption can never merge into it. The expand-state lookup MUST use the same caption,
    // or the group collapses on every refresh.
    if (standalone.length > 0) {
      const aGroup = treeNode.group(UNGROUPED_SECTION, null,
        hadGroups ? wasExpanded.has(UNGROUPED_SECTION) : false);
      for (const container of standalone) {
        const entry = container.analyses[0];
        const meta: NcaNodeMeta = {kind: 'analysis', analysisId: entry.id};
        const leaf = aGroup.item(entry.label, meta);
        wireOpenOnDoubleClick(leaf, () => openAnalysis(entry.id, 'view'));
        wireStandaloneActions(leaf, container);
      }
    }
  };

  await refresh();

  // Stay live: rebuild whenever a study is saved / deleted (Hit Design lacks this).
  ncaTreeSub?.unsubscribe();
  ncaTreeSub = grok.events.onCustomEvent(NCA_STUDIES_CHANGED_EVENT).subscribe(() => {void refresh();});

  // Single-click / arrow-key PREVIEW. One subscription on the tree root (not per-node
  // `onSelected`) so the platform's built-in ↑/↓/←/→ navigation keeps working — a preview
  // fired on every node would otherwise steal focus off the tree (the focus theft that made
  // arrow keys feel dead; TODO.md:1447/1448). Debounced so a double-click cancels it
  // (→ full open) and an arrow sweep doesn't preview every node. Unsubscribe-then-resubscribe
  // per (re)build, like `ncaTreeSub` (DiffStudio `setupTreeItemSelection`).
  ncaTreeSelSub?.unsubscribe();
  ncaTreeSelSub = treeNode.rootNode.onSelectedNodeChanged.subscribe((node) => {
    // The tree outlives a single build; self-detach once its DOM is gone.
    if (!treeNode.root.isConnected) {
      ncaTreeSelSub?.unsubscribe();
      ncaTreeSelSub = null;
      return;
    }
    cancelPendingPreview();
    // Ignore selections outside our subtree (other apps share the one browse tree).
    if (!node || !treeNode.root.contains(node.root))
      return;
    // Keep the tree focused after the preview mounts, so arrow keys keep targeting it.
    armTreeFocusRestore(treeNode);
    treeWalkActive = true;
    const meta = asNcaNodeMeta(node.value);
    if (meta === null)
      return;
    // The mount='preview' openers each clear the preview slot (grok.shell.preview = null)
    // right before their own addPreview, so a new peek REPLACES the previous preview — a prior
    // peek OR the app-landing preview (DiffStudio releaseSharedPreview → addPreview). No
    // per-view tracking / .close() here: .close() on the app-landing view closes the app.
    pendingPreviewTimer = window.setTimeout(() => {
      pendingPreviewTimer = null;
      if (meta.kind === 'study') {
        void openStudyOverview(meta.studyId, 'preview');
        return;
      }
      const opening = openAnalysis(meta.analysisId, 'preview');
      // Resolves before the preview's compute (deferred ≥ 50 ms in the open path) can present,
      // so the identity check below sees this view. A failed open just leaves no preview.
      void opening.then((v) => {lastAnalysisPreview = v;}).catch(() => {lastAnalysisPreview = null;});
    }, PREVIEW_DEBOUNCE_MS);
  });

  // The preview's results arrive after its compute; give the keyboard back to the tree then, as
  // long as the analyst is still walking it. Deferred one tick so it lands after the analysis
  // view's own post-present side-panel reset (nca-mapping-pane.ts, keepBrowseInFront).
  ncaTreeRoot = treeNode.rootNode.root;
  armTreeWalkEnd();
  ncaResultsPresentedSub?.unsubscribe();
  ncaResultsPresentedSub = grok.events.onCustomEvent(NCA_RESULTS_PRESENTED_EVENT).subscribe((v) => {
    if (!treeWalkActive || v !== lastAnalysisPreview || !treeNode.root.isConnected) return;
    setTimeout(() => {
      treeNode.rootNode.root.focus();
      guardTreeFocus(treeNode);
    }, 0);
  });
}

/** A mousedown outside the Browse tree ends the keyboard walk (the analyst has moved to the view).
 *  Registered ONCE, on the first tree build (not at module import); it reads the latest
 *  `ncaTreeRoot`, so a rebuild needs no re-registration. */
let treeWalkEndArmed = false;
function armTreeWalkEnd(): void {
  if (treeWalkEndArmed) return;
  treeWalkEndArmed = true;
  document.addEventListener('mousedown', (e) => {
    if (ncaTreeRoot === null || !(e.target instanceof Node) || !ncaTreeRoot.contains(e.target))
      treeWalkActive = false;
  }, true);
}
