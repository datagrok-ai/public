/**
 * Platform verification for the post-Compute dock layout (user-directed
 * 2026-08-13): tables in one group, charts in a second group split off its top.
 *
 * **Why this file exists.** The layout rests on platform behaviour that cannot be
 * read off the typings — and the first draft of it was reverted precisely because
 * that behaviour was about to be assumed. The charts group has NO permanent
 * member (every chart is conditionally closed: `Mean ± SD` on an empty profile,
 * `Combined` on a cleared mapping, `Facet` on a subject-set change, `Composite`
 * when the design stops being sparse), so it must be addressed through whichever
 * of its tabs is still docked. Getting "still docked" wrong opens a SECOND charts
 * row instead of raising — a layout the user never asked for, reported as
 * success.
 *
 * The PROBE below is the load-bearing test: it establishes, by measurement, which
 * of three plausible liveness signals actually tells the truth. Two of the three
 * are wrong, and the intuitive one (`isConnected`) is wrong in the direction that
 * silently breaks the layout.
 *
 * Measured, not inferred (rule 16): every check reads DOM / dock state back from
 * a real `TableView` on the live host.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  isDockedPanel, buildScrollingFacetGrid, FACET_TILE_MIN_PX, scrollFacetTileIntoView,
  buildFacetGrid, getMaxFacetCols, buildFacetDataFrames, applyFacetLambdaZLine,
  activateMeanSdTab, MEAN_SD_CHART_NODE_TAG,
} from '../views/nca-mapping-pane';
import {buildFacetVerdictColors, buildFacetLambdaZParams} from '../views/results-view';
import {STATUS_ACCEPTED, STATUS_FLAGGED, STATUS_NOT_ESTIMATED} from '../views/status-palette';
import type {FacetLambdaZParams} from '../views/results-view';
import {attachBannerHost} from '../views/compute-notices';

/** A throwaway host view plus a labelled panel factory. */
function host(): {view: DG.TableView; panel: (t: string) => HTMLElement; close: () => void} {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', [1, 2, 3]),
  ]);
  const view = grok.shell.addTableView(df);
  return {
    view,
    panel: (t: string) => ui.divText(t),
    close: () => {
      try {
        view.close();
      } catch {
        // already closed — teardown is best-effort
      }
    },
  };
}

category('Dock layout (2026-08-13)', () => {
  test('the results grid FILLs onto view.grid as a TAB, not a split', async () => {
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      expect(mainNode !== null, true, 'the main grid node resolves before any dock');

      const el = h.panel('NCA results');
      const node = h.view.dockManager.dock(el, DG.DOCK_TYPE.FILL, mainNode, 'NCA results');
      expect(node !== null && node !== undefined, true, 'FILL returned a dock node');

      // The measurable signature of a TAB GROUP, established by the PROBE below:
      // only the ACTIVE tab's content stays attached to the document. The panel
      // just docked is active, so the grid it filled onto must now be detached.
      // A SPLIT would leave both on screen — which is exactly the old 'down'
      // layout this change replaces.
      expect(el.isConnected, true, 'the newly docked results panel is the active tab');
      expect(h.view.grid.root.isConnected, false,
        'the source-data grid went to the background — i.e. FILL produced a tab ' +
        'group, not a split (a split would keep both attached)');
      expect(isDockedPanel(h.view, h.view.grid.root), true,
        'the backgrounded grid is still DOCKED — detached from the DOM is not ' +
        'the same as closed, which is the whole point of isDockedPanel');

      return 'FILL onto view.grid produced a tab group';
    } finally {
      h.close();
    }
  });

  test('PROBE: which liveness signal tells the truth about a docked panel', async () => {
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      const elA = h.panel('A');
      const a = h.view.dockManager.dock(elA, DG.DOCK_TYPE.TOP, mainNode, 'A', 0.5);

      // --- one panel, unambiguously alive and active ---
      expect(elA.isConnected, true, 'a lone docked panel is attached');
      expect(isDockedPanel(h.view, elA), true,
        'findNode resolves the element that was passed to dock()');
      expect((h.view.dockManager.findNode(a.container.containerElement) ?? null) !== null, false,
        'findNode does NOT resolve container.containerElement — a healthy panel ' +
        'reads as absent through it, so it must never be used as the lookup key');

      // --- a second tab joins the SAME group; A is now the background tab ---
      const elB = h.panel('B');
      h.view.dockManager.dock(elB, DG.DOCK_TYPE.FILL, a, 'B');

      // THE decisive measurement. If isConnected meant "alive", A would read true
      // here. It reads false — it means "is the active tab" — so a liveness check
      // built on it prunes every background chart, and once it prunes them all,
      // dockChartTab re-splits and opens a second charts row.
      expect(elA.isConnected, false,
        'an INACTIVE tab is detached from the document — isConnected is an ' +
        'active-tab test, not a liveness test');
      expect(elB.isConnected, true, 'the newly docked tab is the active one');
      expect(isDockedPanel(h.view, elA), true,
        'the background tab is still docked — the signal dockChartTab relies on');

      // --- and it flips honestly when the panel is really gone ---
      h.view.dockManager.close(a);
      expect(isDockedPanel(h.view, elA), false, 'a CLOSED panel reads as not docked');

      return 'isConnected = active-tab; findNode(dockedElement) = liveness';
    } finally {
      h.close();
    }
  });

  test('the control-strip wrapper pushed the grid past its dock pane; the bare host does not', async () => {
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      const frame = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.INT, 'k',
          Array.from({length: 60}, (_, i) => i)),
      ]);

      // The SHIPPED wrapper, built from the real pieces it was built from — a real
      // grid inside a real `attachBannerHost` (a `ui.splitV`, not a plain div),
      // inside the `divV([strip, body])` column. A stand-in div does NOT reproduce
      // this: a plain flex item just shrinks.
      const wrappedGrid = DG.Viewer.grid(frame);
      const strip = ui.divH([ui.divText('Full')]);
      strip.style.flex = '0 0 auto';
      strip.style.padding = '2px 8px';
      const wrapBody = attachBannerHost(wrappedGrid.root);
      const wrapped = ui.divV([strip, wrapBody]);
      wrapped.style.height = '100%';
      wrapped.style.width = '100%';
      wrapped.style.minHeight = '0';
      h.view.dockManager.dock(wrapped, DG.DOCK_TYPE.TOP, mainNode, 'wrapped', 0.5);

      // The SHIPPED fix: the bare banner host every sibling panel docks.
      const bareGrid = DG.Viewer.grid(frame.clone());
      const bare = attachBannerHost(bareGrid.root);
      h.view.dockManager.dock(bare, DG.DOCK_TYPE.FILL, mainNode, 'bare');

      // Let the dock manager lay both out before measuring.
      await new Promise((r) => window.setTimeout(r, 300));

      // How far the grid's bottom edge sits BELOW its own dock pane's bottom —
      // the on-screen defect, measured directly rather than inferred from a
      // scrollHeight proxy.
      const escape = (panel: HTMLElement, gridRoot: HTMLElement): number => {
        const pane = panel.getBoundingClientRect();
        const g = gridRoot.getBoundingClientRect();
        return Math.round(g.bottom - pane.bottom);
      };
      const wrappedEscape = escape(wrapped, wrappedGrid.root);
      const bareEscape = escape(bare, bareGrid.root);
      const stripH = Math.round(strip.getBoundingClientRect().height);

      // Measured on the real host: wrappedEscape=20, bareEscape=0, strip=20.
      //
      // The escape is EXACTLY the strip's height, which is what makes deleting the
      // strip the fix rather than a workaround. Note the shape that produces it:
      // `attachBannerHost` is a `ui.splitV` carrying its own `height: 100%`, so it
      // does NOT resolve the overage by flex-shrink the way a plain child would —
      // a first attempt at this test used stand-in divs, measured 0px, and would
      // have "disproved" a real defect. Keep the real components here.
      expect(wrappedEscape > 0, true,
        `the wrapped panel pushes the grid past its pane — the defect ` +
        `(escape=${wrappedEscape}px)`);
      expect(Math.abs(wrappedEscape - stripH) <= 2, true,
        `the escape equals the control strip's height, pinning the root cause ` +
        `(escape=${wrappedEscape}px, strip=${stripH}px)`);
      expect(bareEscape <= 0, true,
        `the bare banner host — what every panel now docks — fits its pane ` +
        `(escape=${bareEscape}px)`);

      return `wrappedEscape=${wrappedEscape}px bareEscape=${bareEscape}px strip=${stripH}px`;
    } finally {
      h.close();
    }
  });

  test('a background tab can be brought back to the front (results must not stay buried)', async () => {
    // `NCA results` docks FIRST — it establishes the tables group every companion
    // table then FILLs into, and the charts row splits off. So the LAST decorator
    // to dock (Dose-Normalized on a dose-ranging study) ends up the active tab and
    // buries the main result. Docking results last is not an option; re-activating
    // it is. This pins the API that does so.
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      const elA = h.panel('NCA results');
      const a = h.view.dockManager.dock(elA, DG.DOCK_TYPE.FILL, mainNode, 'NCA results');
      const elB = h.panel('Dose-Normalized');
      h.view.dockManager.dock(elB, DG.DOCK_TYPE.FILL, a, 'Dose-Normalized');

      // Established by the PROBE above: exactly the ACTIVE tab is attached.
      expect(elA.isConnected, false, 'the later tab buried the results tab');
      expect(elB.isConnected, true, 'the later tab is the active one');

      // The call under test.
      a.parent.container.setActiveChild(a.container);

      expect(elA.isConnected, true,
        'setActiveChild(resultsNode.container) on the PARENT container brings the ' +
        'results tab back to the front');
      expect(elB.isConnected, false, 'and the companion table goes to the background');
      expect(isDockedPanel(h.view, elB), true,
        'backgrounded, not closed — the companion table is still docked');

      return 'results tab re-activation verified';
    } finally {
      h.close();
    }
  });

  test('activateMeanSdTab fronts the tab stored under MEAN_SD_CHART_NODE_TAG (Group-by → Mean±SD)', async () => {
    // The Group-by handler calls activateMeanSdTab to bring the Mean±SD chart forward
    // when a non-empty grouping is chosen. It reads the chart's dock node from
    // view.temp[MEAN_SD_CHART_NODE_TAG] (set by ensureMeanSdChart) and fronts it — the
    // same setActiveChild mechanism the results-tab test above pins, here reached
    // through the tag the chart dock stores. This covers the new helper end to end.
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      // Mean±SD docks FIRST and is buried by a later chart — exactly the post-Compute
      // shape (charts group: Mean±SD, then Combined/Facet FILL on top).
      const elMsd = h.panel('Mean ± SD');
      const msd = h.view.dockManager.dock(elMsd, DG.DOCK_TYPE.TOP, mainNode, 'Mean ± SD', 0.5);
      const elOther = h.panel('Combined');
      h.view.dockManager.dock(elOther, DG.DOCK_TYPE.FILL, msd, 'Combined');
      expect(elMsd.isConnected, false, 'the later chart buried Mean ± SD');

      // Wire the node under the tag the way ensureMeanSdChart does, then front it.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (h.view.temp as Record<string, any>)[MEAN_SD_CHART_NODE_TAG] = msd;
      activateMeanSdTab(h.view);
      expect(elMsd.isConnected, true, 'activateMeanSdTab brought Mean ± SD to the front');
      expect(elOther.isConnected, false, 'the sibling chart went to the background');
      expect(isDockedPanel(h.view, elOther), true, 'backgrounded, not closed');

      // No node under the tag → cosmetic no-op, never a throw (pre-Compute / tab closed).
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (h.view.temp as Record<string, any>)[MEAN_SD_CHART_NODE_TAG] = null;
      activateMeanSdTab(h.view); // must not throw

      return 'activateMeanSdTab fronts the tagged chart tab; null node is a no-op';
    } finally {
      h.close();
    }
  });

  test('the facet grid past one viewport keeps tile height and SCROLLS', async () => {
    // The facet grid no longer refuses above 20 subjects; it scrolls. That claim
    // is only worth making if the tile really keeps its 5×4 height instead of
    // being squeezed — which is exactly what `ui.splitV` would have done and why
    // this layout is a CSS grid instead. Measured on the real host.
    const h = host();
    const plots: DG.Viewer[] = [];
    try {
      const N = 60; // three viewports' worth
      for (let i = 0; i < N; i++) {
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [0, 1, 2, 4, 8]),
          DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', [10, 8, 6, 3, 1]),
        ]);
        // The same layout-critical options `buildFacetGrid` gives its plots —
        // same viewer type (scatterPlot) and, in particular, `autoLayout: false`,
        // a plausible reason a tile would refuse to shrink that a bare default
        // viewer would not reproduce.
        plots.push(DG.Viewer.scatterPlot(df, {
          xColumnName: 't', yColumnName: 'c',
          autoLayout: false, showXAxis: true, showYAxis: true,
          showXSelector: false, showYSelector: false,
          markerType: DG.MARKER_TYPE.CIRCLE, yAxisType: 'linear' as const,
          linesOrderColumnName: 't', linesWidth: 2,
          title: `S${i}`, showTitle: true,
          description: `S${i}`, descriptionPosition: 'Bottom' as const,
          descriptionVisibilityMode: 'Always' as const,
        }));
      }
      const root = buildScrollingFacetGrid(plots);
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      h.view.dockManager.dock(root, DG.DOCK_TYPE.FILL, mainNode, 'Facet');
      // Let the dock lay out and the ResizeObserver fire.
      await new Promise((r) => window.setTimeout(r, 500));


      expect(root.children.length, N, 'every plot is in the grid — none dropped');
      const tracks = getComputedStyle(root).gridTemplateColumns.split(' ');
      const cols = tracks.length;
      expect(cols, 5, `five columns side by side; got ${cols}`);

      // HORIZONTAL fit. Reported defect: the 4th and 5th columns were sliced off
      // the right edge with no way to reach them — `overflow-x: hidden` CLIPS.
      // Five 1fr tracks can only overflow if something inside refuses to shrink,
      // so measure the track width against the pane rather than trusting `1fr`.
      const itemW = (root.children[0] as HTMLElement).getBoundingClientRect().width;
      expect(root.scrollWidth <= root.clientWidth + 1, true,
        `the grid must fit its pane horizontally: scrollWidth=${root.scrollWidth} ` +
        `clientWidth=${root.clientWidth} tracks=[${tracks.join('|')}] item0=${itemW}`);

      const viewportH = root.clientHeight;
      const rowPx = parseFloat(getComputedStyle(root).gridAutoRows);
      expect(viewportH > 0, true, 'the docked grid has a laid-out height');
      // THE claim, in two halves. A squeezed layout would put 12 rows in the pane
      // (rowPx ≈ h/12); a blindly-quartered one would put four at any pane height,
      // which on a 686px pane meant 171px tiles with their axes stripped off.
      expect(rowPx >= FACET_TILE_MIN_PX, true,
        `a tile is never crushed below the height the chart needs for its axes: ` +
        `got ${rowPx}px, floor ${FACET_TILE_MIN_PX}px (pane ${viewportH}px)`);
      const quarter = Math.floor(viewportH / 4);
      expect(rowPx, Math.max(FACET_TILE_MIN_PX, quarter),
        `a tile fills a quarter of the pane when that clears the floor: pane ` +
        `${viewportH}px → quarter ${quarter}px, floor ${FACET_TILE_MIN_PX}px`);
      // And the surplus is reachable rather than clipped away.
      expect(root.scrollHeight > root.clientHeight, true,
        `60 plots in a 20-plot viewport must overflow into scroll: ` +
        `scrollHeight=${root.scrollHeight} clientHeight=${root.clientHeight}`);

      return `N=${N} cols=${cols} rowPx=${rowPx} viewport=${viewportH} ` +
        `scrollHeight=${root.scrollHeight}`;
    } finally {
      for (const p of plots) {
        try {
          p.close();
        } catch {/* best-effort */}
      }
      h.close();
    }
  });

  test('selecting a results row scrolls its facet tile into view, and only that container', async () => {
    // Dimming the unselected tiles is a half-answer once the grid scrolls: the
    // highlighted tile is usually below the fold, so a click lit up something the
    // scientist could not see. Two claims, and the second is the one that makes
    // this safe in a docked app — the page must not move.
    const h = host();
    const plots: DG.Viewer[] = [];
    try {
      const N = 60;
      for (let i = 0; i < N; i++) {
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [0, 1, 2, 4, 8]),
          DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', [10, 8, 6, 3, 1]),
        ]);
        const p = DG.Viewer.scatterPlot(df, {xColumnName: 't', yColumnName: 'c'});
        p.root.dataset.ncaFacetSubject = `S${i}`;
        plots.push(p);
      }
      const root = buildScrollingFacetGrid(plots);
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;
      h.view.dockManager.dock(root, DG.DOCK_TYPE.FILL, mainNode, 'Facet');
      await new Promise((r) => window.setTimeout(r, 500));

      expect(root.scrollHeight > root.clientHeight, true, 'the grid genuinely scrolls');
      const pageBefore = window.scrollY;

      // A tile far below the fold — row 10 of 12.
      const far = plots[50].root;
      scrollFacetTileIntoView(far);
      await new Promise((r) => window.setTimeout(r, 100));
      const cRect = root.getBoundingClientRect();
      const fRect = far.getBoundingClientRect();
      expect(root.scrollTop > 0, true,
        `the container scrolled to reach a far tile; scrollTop=${root.scrollTop}`);
      expect(fRect.top >= cRect.top - 1 && fRect.bottom <= cRect.bottom + 1, true,
        `the tile is now inside the viewport: tile ${Math.round(fRect.top)}..` +
        `${Math.round(fRect.bottom)} vs container ${Math.round(cRect.top)}..` +
        `${Math.round(cRect.bottom)}`);
      expect(window.scrollY, pageBefore,
        'ONLY the facet container moved — scrollIntoView would have walked every ' +
        'scrollable ancestor and dragged the page with it');

      // A tile already on screen must not jump.
      const settled = root.scrollTop;
      const visible = plots[51].root; // same row as 50
      scrollFacetTileIntoView(visible);
      expect(root.scrollTop, settled, 'an already-visible tile causes no movement');

      return `scrollTop=${root.scrollTop} of ${root.scrollHeight - root.clientHeight}`;
    } finally {
      for (const p of plots) {
        try {
          p.close();
        } catch {/* best-effort */}
      }
      h.close();
    }
  });

  test('a chart docked after the group anchor closed joins the survivors, not a new row', async () => {
    const h = host();
    try {
      const mainNode = h.view.dockManager.findNode(h.view.grid.root) ?? null;

      // Build the charts group the way `dockChartTab` does: the first member
      // splits, the rest FILL onto a surviving member.
      const elA = h.panel('Mean ± SD');
      const a = h.view.dockManager.dock(elA, DG.DOCK_TYPE.TOP, mainNode, 'Mean ± SD', 0.5);
      const elB = h.panel('Combined');
      h.view.dockManager.dock(elB, DG.DOCK_TYPE.FILL, a, 'Combined');

      // Close the ANCHOR while a sibling still lives — the Group-by → "(none)"
      // shaped scenario that made the first draft of this layout unsound.
      h.view.dockManager.close(a);
      expect(isDockedPanel(h.view, elA), false, 'the closed anchor is pruned');
      expect(isDockedPanel(h.view, elB), true, 'the surviving member is still docked');

      // Re-pick exactly as dockChartTab does: first remembered panel that the
      // dock manager still resolves.
      const survivors = [elA, elB].filter((e) => isDockedPanel(h.view, e));
      expect(survivors.length, 1, 'exactly one member survives');
      const elC = h.panel('Facet');
      h.view.dockManager.dock(
        elC, DG.DOCK_TYPE.FILL, h.view.dockManager.findNode(survivors[0])!, 'Facet');

      // Same group ⇒ exactly ONE of them is attached (tab groups show one tab at
      // a time — established by the PROBE). Two separate groups would each show
      // their own active tab, so BOTH would be attached. That is the defect this
      // asserts against: a second charts row opening beside the first.
      const attached = [elB, elC].filter((e) => e.isConnected).length;
      expect(attached, 1,
        'the new chart joined the SURVIVING group — two attached panels would ' +
        'mean a second charts row was opened instead');
      expect(isDockedPanel(h.view, elB), true, 'the survivor is still docked');
      expect(isDockedPanel(h.view, elC), true, 'the new chart is docked');

      return 'group integrity survives anchor death';
    } finally {
      h.close();
    }
  });

  test('buildFacetGrid: options + categorical colour + trimmed perimeter (facet-scatterplot-lines)', async () => {
    // Asserts the REAL buildFacetGrid output, not a hand-copied mirror (CLAUDE.md
    // rule 16). Each subject renders as a scatterPlot whose points are JOINED into a
    // line by `linesOrderColumnName`, and whose per-tile hue rides on a CATEGORICAL
    // `color` column keyed by a `subj '<id>'` token — a bare numeric subject made DG
    // infer a NUMERIC column, silently disabling colour-coding (the all-tiles-blue bug).
    const mk = (subject: string): {subject: string; df: DG.DataFrame} => ({
      subject,
      df: DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [0, 1, 2, 4, 8]),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', [10, 8, 6, 3, 1]),
      ]),
    });
    // Seven subjects → split layout (≤ one viewport); getMaxFacetCols(7) = 4 columns,
    // so row 1 (tiles 4–6) is a PARTIAL row — the case that exercises the trim math.
    const subjects = ['1', '2', '3', '4', '5', '6', '7'];
    const {plots} = buildFacetGrid(subjects.map(mk), 't', 'c');
    expect(plots.length, 7, 'one tile per subject');

    // (1) Line + colour options on every tile.
    for (let i = 0; i < plots.length; i++) {
      const look = plots[i].getOptions().look;
      expect(look.linesOrderColumnName, 't', `tile ${i}: points ordered into a line by time`);
      expect(look.colorColumnName, 'color', `tile ${i}: colour bound to the colour column`);
      expect(look.yColumnName, 'c', `tile ${i}: y bound to the conc column (singular)`);
      expect(look.legendVisibility, 'Never', `tile ${i}: colour legend hidden`);
    }

    // (2) The colour column is CATEGORICAL (the load-bearing fix) and its value is the
    // `subj '<id>'` token — a numeric-inferred column is the bug this guards against.
    for (let i = 0; i < plots.length; i++) {
      const col = plots[i].dataFrame.col('color');
      expect(col !== null, true, `tile ${i}: has the colour column`);
      expect(col!.type, DG.COLUMN_TYPE.STRING, `tile ${i}: colour column is STRING, not numeric`);
      expect(col!.get(0), `subj '${subjects[i]}'`, `tile ${i}: category token`);
      expect(col!.colors.getType(), DG.COLOR_CODING_TYPE.CATEGORICAL,
        `tile ${i}: categorical colour-coding is set`);
    }

    // (3) Perimeter trim (tagFacetEdgeTiles, cols = 4 for n = 7 → lastRow = 1). The
    // internal separators are now COMPLETE across the partial last row: a tile keeps
    // its RIGHT border unless it is the rightmost column, and its BOTTOM border unless
    // it is the last row — so the empty row-1/col-3 slot is framed by tile 3's bottom
    // and tile 6's right, and only the true outer edges are trimmed.
    const cl = (i: number): DOMTokenList => plots[i].root.classList;
    expect(cl(3).contains('nca-facet-last-col') && !cl(3).contains('nca-facet-last-row'), true,
      'tile 3 (row-0 rightmost column) trims the RIGHT edge only — its bottom separator spans the partial row');
    expect(cl(4).contains('nca-facet-last-row') && !cl(4).contains('nca-facet-last-col'), true,
      'tile 4 (last row, not last column) trims the BOTTOM edge only');
    expect(cl(6).contains('nca-facet-last-row') && !cl(6).contains('nca-facet-last-col'), true,
      'tile 6 (last filled tile, col 2 of 4) trims BOTTOM only — its right separator frames the empty corner');
    expect(!cl(0).contains('nca-facet-last-col') && !cl(0).contains('nca-facet-last-row'), true,
      'tile 0 (interior) keeps both separators');

    return 'buildFacetGrid: line/colour options + categorical colour + trimmed perimeter';
  });
});

/**
 * Unit coverage for the facet λz-overlay helpers (spike `facet-lambdaz-line-points`,
 * VR-10). Backfilled at the review gate on user direction — the spike shipped these
 * as throwaway-visual helpers with no dedicated tests. Pure DataFrame→value logic,
 * so each test builds a minimal frame and asserts the returned value/column.
 */
category('Facet λz overlay (facet-lambdaz-line-points)', () => {
  test('getMaxFacetCols: per-N split column table (B5)', async () => {
    // User-directed table 1→1, 2→2, 3→3, 4→2, 5→3, 6→3, 7→4, 8→4, 9→5, 10→5.
    const expected = [1, 2, 3, 2, 3, 3, 4, 4, 5, 5];
    for (let n = 1; n <= 10; n++)
      expect(getMaxFacetCols(n), expected[n - 1], `N=${n} → ${expected[n - 1]} columns`);
    return 'getMaxFacetCols table';
  });

  test('applyFacetLambdaZLine: display anchor + slope, semi-log formula (B1)', async () => {
    // Points laid EXACTLY on log10(C) = 3 − 1·t over the window, so the mean anchor
    // recovers a = 3 with slopeLog = λz/ln10 = 1 (λz = ln10). t=1 is OUTSIDE [2,8] and
    // must be ignored; it sits off the line so its inclusion would move the anchor.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Time', [1, 2, 4, 8]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'log(Conc)', [99, 1, -1, -5]),
    ]);
    const lz: FacetLambdaZParams = {lambdaZ: Math.LN10, intercept: 0, tFirst: 2, tLast: 8};
    applyFacetLambdaZLine(df, lz, 'Conc', 'Time');
    const items = df.meta.formulaLines.items;
    expect(items.length, 1, 'one λz line');
    expect(items[0].formula, '${log(Conc)} = 3 - 1 * ${Time}',
      'anchor a=3 (mean over window points), slope λz/ln10=1, log column on the left');
    expect(items[0].min, 2, 'line bounded to tFirst');
    expect(items[0].max, 8, 'line bounded to tLast');

    // No line when the fit failed (λz ≤ 0) …
    applyFacetLambdaZLine(df, {lambdaZ: 0, intercept: 0, tFirst: 2, tLast: 8}, 'Conc', 'Time');
    expect(df.meta.formulaLines.items.length, 0, 'no line when λz ≤ 0');
    // … or the window holds no plotted point.
    applyFacetLambdaZLine(df, {lambdaZ: 1, intercept: 0, tFirst: 100, tLast: 200}, 'Conc', 'Time');
    expect(df.meta.formulaLines.items.length, 0, 'no line when the window has no point');
    return 'applyFacetLambdaZLine anchor + guards';
  });

  test('buildFacetDataFrames: _lz_role = in-window AND measurable (B2)', async () => {
    // S1 has a t=4 BLQ (conc=0) INSIDE the window — must be 'excluded', not mislabelled
    // 'included'; t=1 is outside the window; S2 has no λz → every point 'excluded'.
    const main = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1', 'S1', 'S2']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Time', [1, 2, 4, 8, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Conc', [100, 50, 0, 10, 25]),
    ]);
    const lzMap = new Map<string, FacetLambdaZParams>([
      ['S1', {lambdaZ: 0.5, intercept: 0, tFirst: 2, tLast: 8}],
      // S2 intentionally absent → NaN window → all excluded.
    ]);
    const frames = buildFacetDataFrames(main, 'Subject', 'Time', 'Conc', lzMap);
    const roleOf = (subject: string): string[] => {
      const entry = frames.find((f) => f.subject === subject)!;
      const col = entry.df.col('_lz_role')!;
      return Array.from({length: entry.df.rowCount}, (_, i) => col.getString(i));
    };
    expect(roleOf('S1').join(','), 'excluded,included,excluded,included',
      't=1 out-of-window; t=2/t=8 in-window+measurable; t=4 in-window BLQ (conc=0) → excluded');
    expect(roleOf('S2').join(','), 'excluded', 'no λz fit → excluded');
    return 'buildFacetDataFrames _lz_role';
  });

  test('buildFacetLambdaZParams: per-subject params, first-seen-wins', async () => {
    // Values chosen to be EXACT in the frame's float storage (dyadic fractions), so the
    // read-back compares cleanly — the point is the keying, not float precision.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SUBJECT', ['A', 'B', 'A']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDA_Z', [0.5, 0.25, 9.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDAZ_INTERCEPT', [1.25, 2.5, 9.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDAZ_TFIRST', [2, 4, 16]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDAZ_TLAST', [8, 16, 32]),
    ]);
    const m = buildFacetLambdaZParams(df);
    const a = m.get('A')!;
    expect(a.lambdaZ, 0.5, 'A λz (FIRST row, not the duplicate)');
    expect(a.intercept, 1.25, 'A intercept');
    expect(a.tFirst, 2, 'A tFirst');
    expect(a.tLast, 8, 'A tLast');
    expect(m.get('B')!.lambdaZ, 0.25, 'B λz');
    return 'buildFacetLambdaZParams';
  });

  test('buildFacetVerdictColors: verdict → tile hue (B3)', async () => {
    // acc: all criteria met; flag: span violated only; ne: all n/a (NaN) → Not estimated.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SUBJECT', ['acc', 'flag', 'ne']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'SPAN_RATIO', [10, 0.1, NaN]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDAZ_ADJ_R2', [0.999, 0.999, NaN]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LAMBDAZ_NPOINTS', [10, 10, NaN]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'PCT_EXTRAP', [1, 1, NaN]),
    ]);
    // The verdict hues are the Radix "teal-success" triad from status-palette.ts (the SINGLE
    // source; user-directed 2026-09-09): Accepted → teal, Flagged → amber, Not estimated → red.
    // Assert against the exported constants, never a literal, so the test cannot drift.
    const m = buildFacetVerdictColors(df);
    expect(m.get('acc'), STATUS_ACCEPTED, 'Accepted → Radix teal-11');
    expect(m.get('flag'), STATUS_FLAGGED, 'Flagged → Radix amber-11');
    expect(m.get('ne'), STATUS_NOT_ESTIMATED, 'Not estimated → Radix red-11');
    return 'buildFacetVerdictColors';
  });
});
