/**
 * Compute notices — collector + banner lifecycle category.
 *
 * @see src/views/compute-notices.ts (implementation under test)
 *
 * Coverage (DOM-level; runs in the grok-test browser harness where `ui` + DOM
 * are available):
 * - B1: collectNotice + flushNotices renders one overlay container per resolved
 *   root, one row per notice, attached to that root.
 * - B3: notices whose anchors resolve to the SAME root (the results-grid
 *   fallback) stack as rows in one container — none dropped.
 * - B4: resetNotices tears down the prior run's banners and clears the
 *   collection, so a recurring notice reappears on the next flush even after it
 *   was dismissed.
 * - Dismiss/restore: the × control collapses its row to the icon (the message
 *   moves to the icon's tooltip) without removing the row or container; clicking
 *   the icon re-expands it — the warning is tucked away, never lost.
 * - refreshNotices (By-group live path): replaces the banner over one root in
 *   place without disturbing banners on other roots; an empty list clears it.
 *
 * A `view` here is just an object with a `temp` bag (the module casts
 * `view.temp` to a plain record), so a stub suffices; roots are real detached
 * DOM elements attached to `document.body` for the duration of each test.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import type * as DG from 'datagrok-api/dg';

import {
  resetNotices, collectNotice, flushNotices, refreshNotices, attachBannerHost,
  collectNoticesFrom,
} from '../views/compute-notices';
import type {ComputeNotice, NoticeAnchor} from '../views/compute-notices';

/** A stub TableView carrying only the `temp` bag the notice module reads. */
function stubView(): DG.TableView {
  return {temp: {} as Record<string, unknown>} as unknown as DG.TableView;
}

function makeRoot(): HTMLElement {
  const el = document.createElement('div');
  document.body.appendChild(el);
  return el;
}

function notice(key: string, anchor: NoticeAnchor, text: string,
  severity: ComputeNotice['severity'] = 'info'): ComputeNotice {
  return {key, anchor, severity, text};
}

function rows(container: HTMLElement): number {
  return container.querySelectorAll('.nca-notice-row').length;
}

category('Compute notices: banner lifecycle', () => {
  test('B1: collect + flush builds one container per root, one row per notice', async () => {
    const view = stubView();
    const fpctRoot = makeRoot();
    try {
      resetNotices(view);
      collectNotice(view, notice('A', 'bioavailability', 'first'));
      collectNotice(view, notice('B', 'bioavailability', 'second', 'warning'));
      // Both anchors resolve to the same (F%) root → one stacked container.
      flushNotices(view, () => fpctRoot);

      const containers = fpctRoot.querySelectorAll('.nca-notice-banner');
      expect(containers.length, 1, 'one container on the F% root');
      expect(rows(fpctRoot), 2, 'two stacked notice rows');
      // Severity is conveyed by the leading icon, not a hue.
      expect(fpctRoot.querySelectorAll('.nca-notice-icon').length, 2, 'each row has an icon');
      expect(fpctRoot.querySelectorAll('.nca-notice-close').length, 2, 'each row has a × control');
      expect(fpctRoot.textContent?.includes('first') ?? false, true, 'first message present');
      expect(fpctRoot.textContent?.includes('second') ?? false, true, 'second message present');
    } finally {
      fpctRoot.remove();
    }
  });

  test('B3: notices for undocked panels fall onto one shared (results) root', async () => {
    const view = stubView();
    const resultsRoot = makeRoot();
    try {
      resetNotices(view);
      // F% N/A + dose-prop underpowered: neither tab docked → both resolve to results.
      collectNotice(view, notice('F', 'bioavailability', 'F% not applicable'));
      collectNotice(view, notice('D', 'doseProp', 'dose-prop underpowered', 'warning'));
      flushNotices(view, () => resultsRoot);

      expect(resultsRoot.querySelectorAll('.nca-notice-banner').length, 1, 'single container');
      expect(rows(resultsRoot), 2, 'both notices stacked on the results root (none dropped)');
    } finally {
      resultsRoot.remove();
    }
  });

  test('AC-7.2: collectNoticesFrom returns fn-collected notices and restores the collector', async () => {
    const view = stubView();
    const root = makeRoot();
    try {
      resetNotices(view);
      collectNotice(view, notice('A', 'bioavailability', 'compute-time notice'));
      // A dependent-aggregate refresh re-invokes decorators that re-collect notices;
      // collectNoticesFrom hands those back for live per-panel refresh WITHOUT leaking
      // them into the Compute's collector.
      const captured = collectNoticesFrom(view, () => {
        collectNotice(view, notice('B', 'bioavailability', 'refresh notice 1', 'warning'));
        collectNotice(view, notice('C', 'doseProp', 'refresh notice 2', 'warning'));
      });
      expect(captured.length, 2, 'both refresh notices are returned to the caller');
      expect(captured.map((n) => n.key).join(','), 'B,C', 'in collection order');
      // The Compute collector is restored → a flush renders only A (B/C were isolated).
      flushNotices(view, () => root);
      expect(rows(root), 1, 'only the Compute-time notice A survives the collector');
      expect(root.textContent?.includes('compute-time notice') ?? false, true, 'A present');
      expect(root.textContent?.includes('refresh notice') ?? false, false, 'B/C not in the Compute set');
    } finally {
      root.remove();
    }
  });

  test('flush with nothing collected is a no-op', async () => {
    const view = stubView();
    const root = makeRoot();
    try {
      resetNotices(view);
      flushNotices(view, () => root);
      expect(root.childElementCount, 0, 'no banner when nothing was collected');
    } finally {
      root.remove();
    }
  });

  test('B4: reset tears down prior banners and a dismissed notice reappears next flush', async () => {
    const view = stubView();
    const root = makeRoot();
    try {
      resetNotices(view);
      collectNotice(view, notice('R', 'sparse', 'recurring note', 'warning'));
      flushNotices(view, () => root);
      expect(rows(root), 1, 'banner shown on first run');

      // User dismisses it → the row collapses to its icon (kept, not removed) so
      // the warning stays reachable; the per-Compute reset is what clears it.
      (root.querySelector('.nca-notice-close') as HTMLElement).click();
      expect(rows(root), 1, 'row kept (collapsed), not removed on dismiss');
      expect(root.querySelectorAll('.nca-notice-collapsed').length, 1, 'row is collapsed to its icon');
      expect(root.querySelectorAll('.nca-notice-banner').length, 1, 'container kept so the icon remains');

      // Next Compute: reset (clears the collection + any leftover DOM), recollect, reflush.
      resetNotices(view);
      collectNotice(view, notice('R', 'sparse', 'recurring note', 'warning'));
      flushNotices(view, () => root);
      expect(rows(root), 1, 'recurring notice reappears regardless of prior dismissal');
      expect(root.querySelectorAll('.nca-notice-collapsed').length, 0, 'reappears expanded, not collapsed');
    } finally {
      root.remove();
    }
  });

  test('dismiss collapses to the icon; clicking the icon restores the message', async () => {
    const view = stubView();
    const root = makeRoot();
    try {
      resetNotices(view);
      collectNotice(view, notice('C', 'sparse', 'collapse me', 'warning'));
      flushNotices(view, () => root);
      expect(root.querySelectorAll('.nca-notice-collapsed').length, 0, 'starts expanded');

      // Dismiss → collapse: the row stays (row count unchanged) but gains the
      // collapsed class, and the icon stays inside it as the restore handle.
      // (An icon-count assertion would be tautological — collapse only toggles a
      // class, so the icon is present whether or not dismiss did anything; assert
      // the icon is INSIDE the collapsed row instead.)
      (root.querySelector('.nca-notice-close') as HTMLElement).click();
      expect(rows(root), 1, 'row kept on dismiss (collapsed, not removed)');
      expect(root.querySelectorAll('.nca-notice-collapsed').length, 1, 'collapsed on dismiss');
      expect(root.querySelectorAll('.nca-notice-collapsed .nca-notice-icon').length, 1,
        'icon stays inside the collapsed row as the restore handle');

      // Click the icon → restore the full row.
      (root.querySelector('.nca-notice-icon') as HTMLElement).click();
      expect(root.querySelectorAll('.nca-notice-collapsed').length, 0, 'icon click re-expands the row');
    } finally {
      root.remove();
    }
  });

  test('B4: reset drops a notice that no longer applies', async () => {
    const view = stubView();
    const root = makeRoot();
    try {
      resetNotices(view);
      collectNotice(view, notice('X', 'sparse', 'transient'));
      flushNotices(view, () => root);
      expect(rows(root), 1, 'shown when it applied');

      // Next Compute: notice no longer collected → no banner after reset+flush.
      resetNotices(view);
      flushNotices(view, () => root);
      expect(root.childElementCount, 0, 'stale notice gone');
    } finally {
      root.remove();
    }
  });

  test('attachBannerHost wraps the viewer so the banner strip sits ABOVE the grid (header not covered)', async () => {
    const view = stubView();
    // A stand-in for a freshly-created, not-yet-docked grid viewer root.
    const viewerRoot = document.createElement('div');
    const wrap = attachBannerHost(viewerRoot);
    document.body.appendChild(wrap);
    try {
      // Structure: a DG ui.splitV whose ABOVE pane holds the banner host and whose
      // below pane holds the viewer (bug grid-banner-header-overlap: the old absolute
      // overlay landed on the grid's own column header and covered it). splitV — not a
      // ui.divV flex column — because only a DG split-layout sizes the viewer pane to
      // full width (a plain flex wrapper collapsed the grid to ~400px; verified live).
      // Assert the structural contract robustly to splitV's internal pane wrappers:
      // the split contains both the host and the viewer, and the host precedes the
      // viewer in DOM order (i.e. the banner strip sits ABOVE the grid).
      const host = wrap.querySelector('.nca-notice-host') as HTMLElement;
      expect(host !== null, true, 'the split contains a banner host strip');
      expect(wrap.contains(viewerRoot), true, 'the split contains the viewer');
      expect(host.contains(viewerRoot), false, 'the viewer is NOT inside the banner strip');
      // eslint-disable-next-line no-bitwise
      const hostBeforeViewer = (host.compareDocumentPosition(viewerRoot) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0;
      expect(hostBeforeViewer, true, 'banner strip precedes the grid in DOM order (renders above it)');

      // A notice resolving to the viewer root must land in the host strip (above the
      // grid), NOT inside the viewer root itself (which would re-cover the header).
      resetNotices(view);
      collectNotice(view, notice('H', 'doseProp', 'header must stay visible', 'warning'));
      flushNotices(view, () => viewerRoot);
      expect(host.querySelectorAll('.nca-notice-banner').length, 1, 'banner rendered in the host strip');
      expect(viewerRoot.querySelectorAll('.nca-notice-banner').length, 0,
        'nothing appended over the grid viewer root — the header is never covered');

      // Idempotent: a second call returns the same split (no double-wrap).
      expect(attachBannerHost(viewerRoot), wrap, 'second call returns the existing split');
    } finally {
      wrap.remove();
    }
  });

  test('attachBannerHost leadingControl: persistent slot above the banners, survives reset/flush', async () => {
    // The results grid passes its Compact toggle as attachBannerHost's leadingControl so
    // the presentation switch sits directly above the grid. Contract: the control lands
    // in a `.nca-results-control-slot` that is the host's FIRST child (above any banner
    // rows), and — unlike a banner container — it is NOT torn down by resetNotices /
    // rebuilt by flushNotices, so the strip stops collapsing to 0 and the toggle persists
    // across every Compute. Other panels pass no control and keep the collapse-to-0 host.
    const view = stubView();
    const viewerRoot = document.createElement('div');
    const control = document.createElement('button'); // stand-in for the toggle root
    control.classList.add('nca-test-compact');
    const wrap = attachBannerHost(viewerRoot, control);
    document.body.appendChild(wrap);
    try {
      const host = wrap.querySelector('.nca-notice-host') as HTMLElement;
      const slot = host.querySelector('.nca-results-control-slot') as HTMLElement;
      expect(slot !== null, true, 'the leading control is wrapped in a control slot');
      expect(slot.contains(control), true, 'the slot holds the passed control');
      expect(host.firstElementChild, slot, 'the slot is the host\'s FIRST child (above banners)');

      // A banner now flushes onto the SAME host — it must append AFTER the slot.
      resetNotices(view);
      collectNotice(view, notice('N', 'doseProp', 'a banner', 'warning'));
      flushNotices(view, () => viewerRoot);
      expect(host.querySelectorAll('.nca-notice-banner').length, 1, 'banner rendered in the host');
      expect(host.firstElementChild, slot, 'the control slot still precedes the banner');
      expect(host.contains(control), true, 'the control survives a flush');

      // reset tears the banner down; the slot (not a tracked banner container) stays.
      resetNotices(view);
      expect(host.querySelectorAll('.nca-notice-banner').length, 0, 'banner cleared on reset');
      expect(host.contains(control), true, 'the control survives reset — the strip stays non-empty');
      expect(host.firstElementChild, slot, 'the slot is still the first (and now only) child');

      // A bare host (no control — every other panel) has no slot.
      const bareRoot = document.createElement('div');
      const bareWrap = attachBannerHost(bareRoot);
      const bareHost = bareWrap.querySelector('.nca-notice-host') as HTMLElement;
      expect(bareHost.querySelector('.nca-results-control-slot'), null,
        'a panel that passes no control has no control slot (collapse-to-0 preserved)');

      return 'leadingControl slot persists above banners across reset/flush; bare host has none';
    } finally {
      wrap.remove();
    }
  });

  test('collapse reparents the badge INTO the grid (corner overlay); expand restores it; reset clears it', async () => {
    // Regression: a collapsed badge left in the sibling banner pane is painted over /
    // clipped by the grid canvas and vanishes. It must move into the viewer root — the
    // pre-split overlay home — so it floats over the grid's corner. (An earlier fix was
    // verified by geometry alone and shipped an invisible badge; this locks the DOM home.)
    const view = stubView();
    const viewerRoot = document.createElement('div');
    const wrap = attachBannerHost(viewerRoot);
    document.body.appendChild(wrap);
    try {
      const host = wrap.querySelector('.nca-notice-host') as HTMLElement;
      resetNotices(view);
      collectNotice(view, notice('C', 'doseProp', 'collapse me', 'warning'));
      flushNotices(view, () => viewerRoot);
      expect(host.querySelectorAll('.nca-notice-row').length, 1, 'expanded row starts in the host strip');
      expect(viewerRoot.querySelectorAll('.nca-notice-row').length, 0, 'nothing in the grid while expanded');

      // Dismiss → the row is REPARENTED into the viewer root as a corner overlay, and
      // the host strip no longer holds it (so the strip collapses and the grid reclaims).
      (wrap.querySelector('.nca-notice-close') as HTMLElement).click();
      expect(viewerRoot.querySelectorAll('.nca-notice-row.nca-notice-collapsed').length, 1,
        'collapsed badge moved into the grid viewer root');
      expect(host.querySelectorAll('.nca-notice-row').length, 0, 'the strip no longer holds the row');

      // Click the icon → the row returns to its original spot in the host strip.
      (viewerRoot.querySelector('.nca-notice-icon') as HTMLElement).click();
      expect(host.querySelectorAll('.nca-notice-row').length, 1, 'expand restores the row to the strip');
      expect(viewerRoot.querySelectorAll('.nca-notice-row').length, 0, 'grid overlay cleared on expand');

      // Collapse again, then reset → the reparented badge must NOT be stranded on the grid.
      (wrap.querySelector('.nca-notice-close') as HTMLElement).click();
      expect(viewerRoot.querySelectorAll('.nca-notice-row').length, 1, 'collapsed again into the grid');
      resetNotices(view);
      expect(viewerRoot.querySelectorAll('.nca-notice-row').length, 0, 'reset clears the reparented badge (no leak)');
      expect(host.querySelectorAll('.nca-notice-row').length, 0, 'reset clears the strip too');
    } finally {
      wrap.remove();
    }
  });

  test('two notices dismissed on one panel STACK as distinct badges (neither is lost)', async () => {
    // Regression (decision-audit FAIL): reparenting each collapsed row to top:0/left:0
    // piled ≥2 badges pixel-identically, hiding all but the last — a warning "lost from
    // view" (violates CLAUDE.md §4). They must reparent into ONE overlay container and
    // stack. meanSd (fragmented + cross-dose-pool) and accumulation (needs-sex + warnings)
    // routinely emit two notices on one anchor, so this is a real case, not a corner.
    const view = stubView();
    const viewerRoot = document.createElement('div');
    const wrap = attachBannerHost(viewerRoot);
    document.body.appendChild(wrap);
    try {
      const overlay = viewerRoot.querySelector('.nca-notice-collapsed-overlay') as HTMLElement;
      expect(overlay !== null, true, 'attachBannerHost created a collapsed-badge overlay in the grid');

      resetNotices(view);
      collectNotice(view, notice('A', 'meanSd', 'first warning', 'warning'));
      collectNotice(view, notice('B', 'meanSd', 'second warning', 'info'));
      flushNotices(view, () => viewerRoot);
      expect(rows(wrap.querySelector('.nca-notice-host') as HTMLElement), 2, 'both rows share the strip');

      // Dismiss BOTH. Each must land in the single overlay as its own child — two
      // distinct, individually-present badges (stacked), not one occluding the other.
      const closes = wrap.querySelectorAll('.nca-notice-close');
      (closes[0] as HTMLElement).click();
      (closes[1] as HTMLElement).click();
      expect(overlay.querySelectorAll('.nca-notice-row.nca-notice-collapsed').length, 2,
        'both dismissed badges are present in the overlay (neither lost)');
      expect(overlay.children.length, 2, 'they are distinct stacked siblings, not piled into one');
      // Each badge keeps its own restore icon — both remain reachable.
      expect(overlay.querySelectorAll('.nca-notice-icon').length, 2, 'each badge has its own restore icon');

      // Restore the first → it returns to the strip; the second badge stays put.
      (overlay.querySelector('.nca-notice-icon') as HTMLElement).click();
      expect(rows(wrap.querySelector('.nca-notice-host') as HTMLElement), 1, 'first row restored to the strip');
      expect(overlay.querySelectorAll('.nca-notice-row').length, 1, 'second badge still in the overlay');
    } finally {
      wrap.remove();
    }
  });

  test('refreshNotices replaces one root in place, leaving other roots untouched', async () => {
    const view = stubView();
    const aggRoot = makeRoot();
    const resultsRoot = makeRoot();
    try {
      resetNotices(view);
      // A flush banner on the results root (e.g. an F% fallback).
      collectNotice(view, notice('F', 'bioavailability', 'results-side note'));
      flushNotices(view, () => resultsRoot);
      // By-group live banner on its own (aggregate) root.
      refreshNotices(view, [notice('M', 'byGroup', 'mixed basis', 'warning')], aggRoot);
      expect(rows(aggRoot), 1, 'by-group banner present');
      expect(rows(resultsRoot), 1, 'results banner present');

      // Group-by changes → re-evaluate: a different mixed-basis message replaces it.
      refreshNotices(view, [notice('M2', 'byGroup', 'mixed dose norm', 'warning')], aggRoot);
      expect(rows(aggRoot), 1, 'still exactly one by-group row (replaced in place)');
      expect(aggRoot.textContent?.includes('mixed dose norm') ?? false, true, 'new message shown');
      expect(aggRoot.textContent?.includes('mixed basis') ?? false, false, 'old message gone');
      expect(rows(resultsRoot), 1, 'results banner untouched by the by-group refresh');

      // Grouping turned off → clear the by-group banner only.
      refreshNotices(view, [], aggRoot);
      expect(aggRoot.childElementCount, 0, 'by-group banner cleared');
      expect(rows(resultsRoot), 1, 'results banner still untouched');
    } finally {
      aggRoot.remove();
      resultsRoot.remove();
    }
  });
});
