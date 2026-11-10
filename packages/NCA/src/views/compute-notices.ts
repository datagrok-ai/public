/**
 * Compute notices — per-Compute collector + persistent dismissible banners.
 *
 * Replaces the flood of transient `grok.shell.info/warning` balloons the analysis
 * decorators fire on Compute (they stacked and auto-dismissed → unreadable) with
 * persistent banners overlaid on top of the panel each notice belongs to.
 * No analytical change — the same notices, better surfaced (honest uncertainty).
 *
 * Channel: each decorator calls {@link collectNotice} (instead of `grok.shell.*`)
 * during the Compute pass; {@link flushNotices} renders the collected notices at
 * the end of the pass. Lifecycle (B4): {@link resetNotices} at the Compute entry
 * tears down the previous run's banners and clears the collection, so a notice
 * that recurs reappears even if the user closed it last time, and a notice that
 * no longer applies disappears.
 *
 * Anchoring (B2/B3): a notice names its *logical* home panel ({@link NoticeAnchor});
 * the consumer supplies an `anchorRoot` resolver that maps the anchor to the home
 * panel's grid/viewer root element (held in `view.temp`), or to the always-present
 * results grid root as the fallback when that panel's tab is not docked.
 *
 * Placement (bug grid-banner-header-overlap): the banner is NOT overlaid on the
 * grid — an absolute top strip lands on the grid's own canvas-drawn column header
 * and covers it. Instead each dockable viewer is wrapped at dock time by
 * {@link attachBannerHost} into a flex column `[bannerHost, viewer@flex:1]`; the
 * banner renders as a thin flow strip in `bannerHost` ABOVE the grid, so the grid
 * (header included) simply starts below it. Still per-panel (travels with the tab,
 * zero global-viewport cost) — the grid loses only the strip's height, never
 * scrolls. `flush`/`refresh` map the resolved viewer root to its registered host
 * internally, so callers keep passing the viewer root.
 *
 * Reuse: models the `errorPanel` / `.nca-mapping-pane-sparse-hint`
 * neutral strip (the conservative default: neutral grey, small text) and the
 * `view.temp`-keyed lifecycle the dock helpers already use. Severity is conveyed by a
 * leading icon, NOT a background hue (audit-checklist C-05/C-27: Info/Warning is a
 * categorical label → neutral; user-approved design decision 2026-06-23).
 */

import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';

type NoticeSeverity = 'info' | 'warning';

/** Logical home panel a notice belongs to. The consumer's `anchorRoot` resolver
 *  maps this to a DOM root (the panel's grid/viewer root in `view.temp`); when the
 *  panel's tab is not docked it falls back to the results grid (B3). */
export type NoticeAnchor =
  | 'byGroup'
  | 'bioavailability'
  | 'doseProp'
  | 'accumulation'
  | 'exposureMargin'
  | 'sexComparison'
  | 'sparse'
  | 'composite'
  | 'meanSd'
  /** The per-subject results grid itself (e.g. the review-decisions-reset notice). */
  | 'results';

export interface ComputeNotice {
  /** Stable identity for the notice (debug / future de-dup); not user-visible. */
  key: string;
  /** Which panel this notice is "about" — drives banner placement. */
  anchor: NoticeAnchor;
  severity: NoticeSeverity;
  /** CT-19 message: fact + consequence + action. Rendered verbatim. */
  text: string;
}

/** Resolves a notice's logical anchor to the DOM root it should overlay. Returns
 *  the home panel's grid/viewer root when docked, else the results grid root
 *  (B3 fallback), else `null` (nothing docked — nothing to overlay). */
export type AnchorRootResolver = (anchor: NoticeAnchor) => HTMLElement | null;

const COLLECTOR_TAG = 'ncaNoticeCollector';
const BANNER_ELS_TAG = 'ncaNoticeBannerEls';

/** Maps a dockable viewer's root → the flow strip that hosts its banners (created
 *  by {@link attachBannerHost}). A WeakMap so a torn-down viewer's entry is GC'd
 *  with it. `flush`/`refresh` translate a resolved viewer root to its host through
 *  {@link hostOrSelf}, so banners land in the strip ABOVE the grid, not over it. */
const bannerHosts = new WeakMap<HTMLElement, HTMLElement>();
/** viewerRoot → the split element returned by {@link attachBannerHost} (docked in
 *  place of the bare root). Kept only for the idempotent re-attach guard. */
const bannerWraps = new WeakMap<HTMLElement, HTMLElement>();
/** viewerRoot → its collapsed-badge overlay: a small absolute container pinned to the
 *  grid's top-left corner INSIDE the viewer root (so it paints over the grid canvas,
 *  which a badge stranded in the sibling split-pane cannot — bug grid-banner-header-
 *  overlap). Collapsed rows are reparented into it and STACK vertically, so several
 *  simultaneously-dismissed notices on one panel stay individually visible + clickable
 *  (the pre-split overlay stacked them the same way; a per-row `top:0` overlay would
 *  pile them pixel-identically and hide all but the last). */
const collapsedOverlays = new WeakMap<HTMLElement, HTMLElement>();
/** banner container → the rows it created. Teardown removes these explicitly because a
 *  collapsed row is reparented OUT of its container (into the collapsed overlay), so
 *  removing the container alone would strand it on the grid. */
const bannerRows = new WeakMap<HTMLElement, HTMLElement[]>();

/** Remove a banner container and every row it created, wherever each row currently
 *  lives (a collapsed row is reparented into the viewer root; `remove()` on an already
 *  detached node is a safe no-op). */
function removeBanner(container: HTMLElement): void {
  const rows = bannerRows.get(container);
  if (rows != null) for (const r of rows) r.remove();
  container.remove();
}

/**
 * Wrap a freshly-created (not-yet-docked) grid/chart viewer root in a DG **vertical
 * split** `ui.splitV([bannerPane, viewerPane])` and register `viewerRoot → bannerHost`.
 * Dock the RETURNED split in place of the bare viewer root; the viewer fills the panel
 * minus a thin top strip that its banners occupy (bug grid-banner-header-overlap: the
 * previous absolute overlay covered the grid's own column header). Idempotent — a
 * second call for an already-wrapped root returns its existing split.
 *
 * Why `ui.splitV` and NOT a `ui.divV`/`ui.box` flex column: a DG grid inside a plain
 * flex wrapper does NOT inherit the dock's sizing — the grid's WIDTH collapses to its
 * intrinsic ~400px (measured: `flex:1` sizes only the column's main axis, and DG pins
 * the grid root's width so `align-items:stretch` can't widen it). `ui.splitV` builds a
 * real DG `d4` split-layout node, which DOES size the viewer pane to full width + the
 * remaining height (measured: grid 1145px wide, banner strip 30-49px tall). The banner
 * pane is `flex:0 0 auto` so it hugs its content — one row, several stacked rows, or a
 * collapsed badge — and shrinks to 0 height when there are no notices.
 *
 * `leadingControl` (optional): a persistent control the caller wants pinned to the TOP
 * of the strip, ABOVE any banner rows. It is placed once, as the host's first child, and
 * is NOT a tracked banner container, so {@link resetNotices}/{@link flushNotices} — which
 * only add/remove banner containers — leave it in place across every Compute; banner
 * containers append after it, so notices stack directly below the control. A strip
 * carrying one therefore no longer shrinks to 0 (it hugs the control's height).
 *
 * As of 2026-09-20 NO nca-studio caller passes `leadingControl` — the results grid's
 * `Compact` toggle used it (user-directed 2026-08-28) until that toggle moved into the
 * ribbon Settings dialog (nca-settings-gear-dialog), so every dock site now calls the
 * 1-arg form. The parameter is retained as generic API (exercised directly by
 * compute-notices-tests.ts); a future cleanup may drop it if no consumer returns.
 */
export function attachBannerHost(
  viewerRoot: HTMLElement, leadingControl?: HTMLElement,
): HTMLElement {
  const existingWrap = bannerWraps.get(viewerRoot);
  if (existingWrap != null) return existingWrap;

  const host = ui.div([], 'nca-notice-host');
  if (leadingControl != null)
    host.appendChild(ui.div([leadingControl], 'nca-results-control-slot'));
  const bannerPane = ui.box(host);
  bannerPane.style.flex = '0 0 auto'; // content-height strip; 0 when no notices
  const viewerPane = ui.box(viewerRoot); // DG box → the viewer fills its pane (full W+H)
  const split = ui.splitV([bannerPane, viewerPane], null, false);
  split.style.height = '100%';
  split.style.width = '100%';

  // Collapsed-badge overlay: an empty (zero-size) corner container INSIDE the viewer
  // root now, populated on dismiss. Living inside the grid root, it paints over the
  // canvas; sized to its content, it never blocks grid clicks while empty.
  const collapsedOverlay = ui.divV([], 'nca-notice-collapsed-overlay');
  viewerRoot.appendChild(collapsedOverlay);

  bannerHosts.set(viewerRoot, host);
  bannerWraps.set(viewerRoot, split);
  collapsedOverlays.set(viewerRoot, collapsedOverlay);
  return split;
}

/** A resolved viewer root → its banner-host strip when wrapped, else the root
 *  itself (tests and any un-wrapped fallback append directly to the given root). */
function hostOrSelf(root: HTMLElement): HTMLElement {
  return bannerHosts.get(root) ?? root;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function tempBagOf(view: DG.TableView): Record<string, any> {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return view.temp as Record<string, any>;
}

/**
 * Start a fresh notice collection for this Compute and tear down the previous
 * run's banners (B4). Call at the Compute entry, before the decorators run. The
 * dismissed state is implicitly per-Compute: closing a banner only collapses it
 * to its icon (the row stays); the next Compute rebuilds every still-relevant
 * banner expanded from scratch.
 */
export function resetNotices(view: DG.TableView): void {
  const bag = tempBagOf(view);
  // `!= null`, not `!== undefined`: `view.temp` returns dart-side null for a
  // never-set key (code-quality-guardrails Datagrok API footgun / BUG-06).
  const prior = bag[BANNER_ELS_TAG] as HTMLElement[] | null;
  if (prior != null)
    for (const el of prior) removeBanner(el); // also clears rows reparented to a grid
  bag[BANNER_ELS_TAG] = [];
  bag[COLLECTOR_TAG] = [];
}

/**
 * Run `fn`, RETURN the notices it contributed, and restore the current Compute's
 * collector to its pre-`fn` state (AC-7.2). Channel: an interactive refresh re-invokes
 * the Compute-time decorators to update their DATA (grid frame swaps) AND wants their
 * freshly-evaluated notices so it can live-replace each panel's banner via
 * {@link refreshNotices} — without those notices leaking into the last Compute's
 * collector (`collectNotice` is a plain push with no dedup, so leaving them there would
 * stack duplicates on the next flush). The caller decides where the returned notices
 * render; this function only isolates them. Assumes `fn` APPENDS to the collector
 * (which `collectNotice` does), so the tail past the pre-`fn` length is exactly `fn`'s
 * contribution.
 */
export function collectNoticesFrom(view: DG.TableView, fn: () => void): ComputeNotice[] {
  const bag = tempBagOf(view);
  const start = ((bag[COLLECTOR_TAG] as ComputeNotice[] | null) ?? []).length;
  fn();
  const live = (bag[COLLECTOR_TAG] as ComputeNotice[] | null) ?? [];
  const added = live.slice(start);
  bag[COLLECTOR_TAG] = live.slice(0, start); // restore — drop fn's additions
  return added;
}

/** Contribute one notice to the current Compute's collection (channel: replaces a
 *  `grok.shell.info/warning` call). No-op-safe before {@link resetNotices}. */
export function collectNotice(view: DG.TableView, notice: ComputeNotice): void {
  const bag = tempBagOf(view);
  const collector = (bag[COLLECTOR_TAG] as ComputeNotice[] | null) ?? [];
  collector.push(notice);
  bag[COLLECTOR_TAG] = collector;
}

/**
 * Render the collected notices as persistent banners over their home panels.
 * Call once at the end of the Compute pass, after all decorators have run and
 * docked their grids (so the anchor roots exist).
 *
 * Notices that resolve to the SAME root (e.g. several "not applicable" notices
 * falling back to the results grid, or By-group's mixed-basis + mixed-dose-norm
 * pair) stack as rows in one overlay container on that root.
 */
export function flushNotices(view: DG.TableView, anchorRoot: AnchorRootResolver): void {
  const bag = tempBagOf(view);
  const notices = (bag[COLLECTOR_TAG] as ComputeNotice[] | null) ?? [];
  if (notices.length === 0) return;

  // Group by the resolved VIEWER root so multiple notices on one panel share a
  // single top-strip container (stacked rows), and B3 fallbacks pile onto the
  // results grid together. Keyed on the viewer root (not the host) so we still hold
  // it at build time — for the collapsed-badge overlay — without a reverse lookup.
  const byViewer = new Map<HTMLElement, ComputeNotice[]>();
  for (const n of notices) {
    const resolved = anchorRoot(n.anchor);
    if (resolved === null) continue; // nothing docked to anchor — drop silently
    const list = byViewer.get(resolved) ?? [];
    list.push(n);
    byViewer.set(resolved, list);
  }

  const created = (bag[BANNER_ELS_TAG] as HTMLElement[] | null) ?? [];
  for (const [viewerRoot, list] of byViewer) {
    // Render into the strip ABOVE the grid (the host), not the grid root itself.
    const container = buildBannerContainer(hostOrSelf(viewerRoot), list,
      collapsedOverlays.get(viewerRoot) ?? null);
    created.push(container);
  }
  bag[BANNER_ELS_TAG] = created;
  bag[COLLECTOR_TAG] = []; // consumed — guard against a stray double-flush
}

/**
 * Replace the banner over ONE specific root in place, without touching banners on
 * other roots. For the By-group panel, whose mixed-basis / mixed-dose-norm notice
 * depends on the live Group-by selection and is re-evaluated on every Group-by
 * change (no recompute) — so it cannot route through the per-Compute collector.
 * The By-group grid uniquely owns its root (the aggregate is only built when a
 * grouping is active), so removing every banner parented to that root and
 * rebuilding is safe and never disturbs the fallback-to-results banners.
 */
export function refreshNotices(
  view: DG.TableView, notices: ComputeNotice[], viewerRoot: HTMLElement | null,
): void {
  if (viewerRoot === null) return;
  // Same viewer-root → banner-host translation as flushNotices (the By-group caller
  // passes the aggregate grid's viewer root; banners belong in the strip above it).
  const root = hostOrSelf(viewerRoot);
  const bag = tempBagOf(view);
  const tracked = (bag[BANNER_ELS_TAG] as HTMLElement[] | null) ?? [];
  const kept: HTMLElement[] = [];
  for (const el of tracked) {
    // The container stays parented to `root` (the host) even when its rows reparent
    // out on collapse, so this per-root filter is intact; removeBanner clears rows too.
    if (el.parentElement === root) removeBanner(el);
    else kept.push(el);
  }
  if (notices.length > 0)
    kept.push(buildBannerContainer(root, notices, collapsedOverlays.get(viewerRoot) ?? null));
  bag[BANNER_ELS_TAG] = kept;
}

/** Build the top-strip container (one per host) holding one dismissible row per
 *  notice, and append it as a flow child of `root` (the banner-host strip that
 *  {@link attachBannerHost} placed above the grid, or — in tests — the given root).
 *  `collapsedOverlay` is the grid-corner container a dismissed row reparents into
 *  (null in tests / un-wrapped fallback, where a row just collapses in place). */
function buildBannerContainer(
  root: HTMLElement, notices: ComputeNotice[], collapsedOverlay: HTMLElement | null,
): HTMLElement {
  const container = ui.divV([], 'nca-notice-banner');
  const rows: HTMLElement[] = [];
  for (const n of notices) {
    const r = buildNoticeRow(n, collapsedOverlay);
    rows.push(r);
    container.appendChild(r);
  }
  bannerRows.set(container, rows);
  root.appendChild(container);
  return container;
}

/**
 * One notice row in its expanded state: leading severity icon + message text +
 * dismiss control. Dismissing does NOT remove the notice — it collapses the row
 * to its icon alone (the message moves to the icon's tooltip), and clicking the
 * icon re-expands it. So a warning can be tucked into the panel corner but is
 * never lost from view (honest uncertainty); the per-Compute reset
 * is what actually clears banners that no longer apply.
 *
 * Collapse relocates the badge (bug grid-banner-header-overlap): the expanded row
 * lives in the banner strip ABOVE the grid (reserves layout so the header stays
 * visible), but on dismiss it is REPARENTED into `collapsedOverlay` — the grid-corner
 * container that {@link attachBannerHost} placed INSIDE the viewer root. A badge left
 * behind in the sibling split-pane is painted over / clipped by the grid canvas and
 * disappears; inside the viewer root (as the pre-split banner did) it paints over the
 * grid's top-left corner, and the now-empty strip collapses to zero so the grid
 * reclaims the space. Badges STACK vertically in that container, so two notices
 * dismissed on the same panel both stay visible + clickable. Re-expand puts the row
 * back at its original strip spot via a placeholder. `collapsedOverlay` is null in
 * tests / un-wrapped fallback, where the row just collapses in place.
 */
function buildNoticeRow(notice: ComputeNotice, collapsedOverlay: HTMLElement | null = null): HTMLElement {
  const icon = ui.iconFA(notice.severity === 'warning' ? 'exclamation-triangle' : 'info-circle');
  icon.classList.add('nca-notice-icon');

  const text = ui.divText(notice.text, 'nca-notice-text');

  const row = ui.divH([icon, text], 'nca-notice-row');

  // Marks the row's home in the strip so re-expand restores its original order.
  let placeholder: Comment | null = null;

  const collapse = (): void => {
    if (row.classList.contains('nca-notice-collapsed')) return;
    row.classList.add('nca-notice-collapsed');
    if (collapsedOverlay === null) return;
    placeholder = document.createComment('nca-notice');
    row.parentElement?.insertBefore(placeholder, row);
    collapsedOverlay.appendChild(row); // stacks under any already-collapsed badge
  };

  const expand = (): void => {
    if (!row.classList.contains('nca-notice-collapsed')) return;
    row.classList.remove('nca-notice-collapsed');
    if (placeholder?.parentElement != null) {
      placeholder.parentElement.insertBefore(row, placeholder);
      placeholder.remove();
      placeholder = null;
    }
  };

  const close = ui.iconFA('times', collapse,
    'Dismiss — collapses to the icon; click the icon to reopen');
  close.classList.add('nca-notice-close');
  row.appendChild(close);

  // The icon doubles as the restore affordance: it re-expands the row on click,
  // and while collapsed it carries the full message as a tooltip so the text
  // stays reachable on hover. The tooltip is gated on the collapsed state (lazy
  // getter returns null when expanded) — an expanded row already shows its text,
  // so a duplicate icon tooltip there would be noise. In the expanded state the
  // click is a harmless no-op (the class isn't set), so binding it is safe.
  ui.tooltip.bind(icon, () => row.classList.contains('nca-notice-collapsed') ? notice.text : null);
  icon.addEventListener('click', expand);
  return row;
}
