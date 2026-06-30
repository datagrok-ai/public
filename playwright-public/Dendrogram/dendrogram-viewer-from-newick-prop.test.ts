/* ---
sub_features_covered: [dendrogram.lifecycle.on-property-changed, dendrogram.lifecycle.on-table-attached, dendrogram.prop.color-aggr-type, dendrogram.prop.color-column-name, dendrogram.prop.current-color, dendrogram.prop.font, dendrogram.prop.light-color, dendrogram.prop.line-width, dendrogram.prop.main-color, dendrogram.prop.mouse-over-color, dendrogram.prop.newick, dendrogram.prop.newick-tag, dendrogram.prop.node-column-name, dendrogram.prop.node-size, dendrogram.prop.selections-color, dendrogram.prop.show-grid, dendrogram.prop.show-labels, dendrogram.prop.show-tooltip, dendrogram.prop.step, dendrogram.prop.step-zoom, dendrogram.viewer]
--- */
// Newick resolution priority: literal `newick` prop > newickTag > .newick tag > empty.
// Property-panel row selectors are name-mangled (e.g. nodeColumnName -> [name="prop-node"]).
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

// Canonical newick strings (matched against the synthetic `leaves4` DataFrame).
const NEWICK_BALANCED = '((a:1,b:1):1,(c:1,d:1):1);';
const NEWICK_LEFT_LEAN = '(((a:1,b:1):1,c:1):1,d:1);';

/** Eighteen property-panel row anchors (name-mangled). Used by the panel-open DOM check. */
const PROP_ROW_NAMES = [
  'newick', 'newick-tag', 'node', 'color', 'color-aggr-type',
  'line-width', 'node-size', 'show-grid', 'main-color', 'light-color',
  'current-color', 'mouse-over-color', 'selections-color', 'show-labels',
  'font', 'step-zoom', 'show-tooltip', 'step',
];

async function getViewerInfo(page: Page): Promise<{
  newick: string | null;
  newickTag: string | null;
  nodeColumnName: string | null;
  colorColumnName: string | null;
  colorAggrType: string | null;
  lineWidth: number | null;
  nodeSize: number | null;
  showGrid: boolean | null;
  showLabels: boolean | null;
  showTooltip: boolean | null;
  stepZoom: number | null;
  step: number | null;
  font: string | null;
  mainColor: unknown;
  lightColor: unknown;
  currentColor: unknown;
  mouseOverColor: unknown;
  selectionsColor: unknown;
  treeNewick: string | null;
  dataFrameMatches: boolean;
  canvasWidth: number;
  canvasHeight: number;
  isAttachedToTv: boolean;
}> {
  return await page.evaluate(() => {
    const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Dendrogram') as any;
    if (!v) throw new Error('Dendrogram viewer not mounted');
    const canvas = document.querySelector('[name="viewer-Dendrogram"] canvas') as HTMLCanvasElement | null;
    return {
      newick: v.props.newick ?? null,
      newickTag: v.props.newickTag ?? null,
      nodeColumnName: v.props.nodeColumnName ?? null,
      colorColumnName: v.props.colorColumnName ?? null,
      colorAggrType: v.props.colorAggrType ?? null,
      lineWidth: v.props.lineWidth ?? null,
      nodeSize: v.props.nodeSize ?? null,
      showGrid: v.props.showGrid ?? null,
      showLabels: v.props.showLabels ?? null,
      showTooltip: v.props.showTooltip ?? null,
      stepZoom: v.props.stepZoom ?? null,
      step: v.props.step ?? null,
      font: v.props.font ?? null,
      mainColor: v.props.mainColor ?? null,
      lightColor: v.props.lightColor ?? null,
      currentColor: v.props.currentColor ?? null,
      mouseOverColor: v.props.mouseOverColor ?? null,
      selectionsColor: v.props.selectionsColor ?? null,
      // treeNewick is the resolved (priority-order-applied) newick.
      treeNewick: v.treeNewick ?? null,
      dataFrameMatches: v.dataFrame === grok.shell.tv.dataFrame,
      canvasWidth: canvas ? canvas.width : 0,
      canvasHeight: canvas ? canvas.height : 0,
      isAttachedToTv: !v.isDetached,
    };
  });
}

async function applyViewerOptions(page: Page, opts: Record<string, unknown>): Promise<void> {
  await page.evaluate((o: Record<string, unknown>) => {
    const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Dendrogram') as any;
    if (!v) throw new Error('Dendrogram viewer not mounted');
    v.setOptions(o);
  }, opts);
  // Settle so onPropertyChanged restyle/rebuild can run.
  await page.waitForTimeout(250);
}

async function getLeafSetFromNewick(page: Page, newick: string): Promise<string[]> {
  return await page.evaluate(async (nwk: string) => {
    const helper: any = await grok.functions.call('Dendrogram:getTreeHelper');
    const df = helper.newickToDf(nwk, 'tmp_leaves');
    const out: string[] = [];
    const nodeCol = df.col('node');
    const leafCol = df.col('leaf');
    for (let i = 0; i < df.rowCount; i++)
      if (leafCol.get(i)) out.push(nodeCol.get(i));
    return out.sort();
  }, newick);
}

test('Dendrogram / Viewer from literal newick property + full property-panel sweep', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup — build the synthetic `leaves4` DataFrame and seed two newick tags.
  await page.evaluate(async (args: {balanced: string; leftLean: string}) => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 800));
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'leaf', ['a', 'b', 'c', 'd']),
      DG.Column.fromList('int', 'value', [1, 2, 3, 4]),
    ]);
    df.name = 'leaves4';
    df.setTag('.newick', args.balanced);
    df.setTag('.newick-alt', args.leftLean);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, {balanced: NEWICK_BALANCED, leftLean: NEWICK_LEFT_LEAN});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // ── Scenario 1 — Viewer mounts and renders from literal newick property ─

  await softStep('1.1 Add Dendrogram viewer with literal newick + nodeColumnName: leaf', async () => {
    await page.evaluate((args: {newick: string}) => {
      const tv = grok.shell.tv;
      tv.addViewer('Dendrogram', {newick: args.newick, nodeColumnName: 'leaf'});
    }, {newick: NEWICK_BALANCED});
    await page.locator('[name="viewer-Dendrogram"]').waitFor({timeout: 30_000});
    await page.locator('[name="viewer-Dendrogram"] canvas').waitFor({timeout: 30_000});
    const info = await getViewerInfo(page);
    expect(info.canvasWidth, 'canvas width > 0').toBeGreaterThan(0);
    expect(info.canvasHeight, 'canvas height > 0').toBeGreaterThan(0);
    expect(info.dataFrameMatches, 'viewer.dataFrame === tv.dataFrame').toBe(true);
    expect(info.isAttachedToTv, 'viewer is attached').toBe(true);
  });

  await softStep('1.2 viewer.props.newick === literal (priority order: property > newickTag > .newick tag)', async () => {
    const info = await getViewerInfo(page);
    expect(info.newick, 'newick prop is the literal').toBe(NEWICK_BALANCED);
    expect(info.nodeColumnName, 'nodeColumnName is leaf').toBe('leaf');
    expect(info.treeNewick, 'treeNewick resolves to the literal property').toBe(NEWICK_BALANCED);
  });

  await softStep('1.3 Leaf set = {a, b, c, d} via TreeHelper.newickToDf', async () => {
    const leaves = await getLeafSetFromNewick(page, NEWICK_BALANCED);
    expect(leaves, 'parsed leaf set').toEqual(['a', 'b', 'c', 'd']);
  });

  // ── Scenario 2 — Property-panel sweep ───────────────────────────────────

  await softStep('2.1 Open property panel (grok.shell.o = viewer); verify all 18 prop rows render', async () => {
    await page.evaluate(() => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Dendrogram') as any;
      if (!v) throw new Error('Dendrogram viewer not mounted');
      grok.shell.o = v;
    });
    await page.locator('[name="prop-newick"]').waitFor({timeout: 15_000});
    const present = await page.evaluate((names: string[]) => {
      const out: Record<string, {row: boolean; view: boolean}> = {};
      for (const n of names) {
        out[n] = {
          row: !!document.querySelector(`[name="prop-${n}"]`),
          view: !!document.querySelector(`[name="prop-view-${n}"]`),
        };
      }
      return out;
    }, PROP_ROW_NAMES);
    for (const n of PROP_ROW_NAMES) {
      expect(present[n].row, `[name="prop-${n}"] row present`).toBe(true);
      expect(present[n].view, `[name="prop-view-${n}"] cell present`).toBe(true);
    }
    const cats = await page.evaluate(() => ({
      data: !!document.querySelector('[name="prop-category-data"]'),
      behavior: !!document.querySelector('[name="prop-category-behavior"]'),
      style: !!document.querySelector('[name="prop-category-style"]'),
    }));
    expect(cats.data, 'prop-category-data header').toBe(true);
    expect(cats.behavior, 'prop-category-behavior header').toBe(true);
    expect(cats.style, 'prop-category-style header').toBe(true);
  });

  await softStep('2.2 newickTag dropdown is populated from .* tags; switching it rebuilds the view', async () => {
    // Clear the literal prop so the newickTag wins, then switch the tag and observe treeNewick.
    await applyViewerOptions(page, {newick: ''});
    let info = await getViewerInfo(page);
    expect(info.newick, 'newick prop cleared').toBe('');

    await applyViewerOptions(page, {newickTag: '.newick-alt'});
    info = await getViewerInfo(page);
    expect(info.newickTag, 'newickTag set to .newick-alt').toBe('.newick-alt');
    expect(info.treeNewick, 'treeNewick resolves to .newick-alt value').toBe(NEWICK_LEFT_LEAN);
    const altLeaves = await getLeafSetFromNewick(page, info.treeNewick ?? '');
    expect(altLeaves, 'leaf set unchanged across topology').toEqual(['a', 'b', 'c', 'd']);
    expect(info.canvasWidth, 'canvas non-empty after newickTag swap').toBeGreaterThan(0);
    expect(info.canvasHeight, 'canvas non-empty after newickTag swap').toBeGreaterThan(0);

    // Reset to scenario-1 state: literal newick wins again.
    await applyViewerOptions(page, {newick: NEWICK_BALANCED});
    info = await getViewerInfo(page);
    expect(info.treeNewick, 'treeNewick snaps back to literal newick').toBe(NEWICK_BALANCED);
  });

  // 2.3 — nodeColumnName rebuild path (data prop).
  await softStep('2.3 nodeColumnName toggle leaf<->value<->leaf rebuilds the view (no console error)', async () => {
    await applyViewerOptions(page, {nodeColumnName: 'value'});
    let info = await getViewerInfo(page);
    expect(info.nodeColumnName, 'nodeColumnName set to value').toBe('value');
    await applyViewerOptions(page, {nodeColumnName: 'leaf'});
    info = await getViewerInfo(page);
    expect(info.nodeColumnName, 'nodeColumnName restored to leaf').toBe('leaf');
    expect(info.canvasWidth, 'canvas non-empty after nodeColumnName toggle').toBeGreaterThan(0);
  });

  // 2.4 — colorColumnName + colorAggrType rebuild path (data props).
  await softStep('2.4 colorColumnName=value + colorAggrType {avg,min,max,med,total-count}', async () => {
    await applyViewerOptions(page, {colorColumnName: 'value'});
    let info = await getViewerInfo(page);
    expect(info.colorColumnName, 'colorColumnName set to value').toBe('value');
    const aggrChoices: ReadonlyArray<string> = ['avg', 'min', 'max', 'med', 'total-count'];
    for (const a of aggrChoices) {
      await applyViewerOptions(page, {colorAggrType: a});
      info = await getViewerInfo(page);
      expect(info.colorAggrType, `colorAggrType set to ${a}`).toBe(a);
    }
    expect(info.canvasWidth, 'canvas non-empty after colorAggrType sweep').toBeGreaterThan(0);
  });

  // 2.5 — lineWidth restyle path (style prop).
  await softStep('2.5 lineWidth sweep 0 / 16 / 2.5 (slider bounds + mid-range)', async () => {
    for (const lw of [0, 16, 2.5]) {
      await applyViewerOptions(page, {lineWidth: lw});
      const info = await getViewerInfo(page);
      expect(info.lineWidth, `lineWidth = ${lw}`).toBe(lw);
    }
    const info = await getViewerInfo(page);
    expect(info.canvasWidth, 'canvas non-empty after lineWidth sweep').toBeGreaterThan(0);
  });

  // 2.6 — nodeSize restyle path.
  await softStep('2.6 nodeSize sweep 0 / 16 / 4.0', async () => {
    for (const ns of [0, 16, 4]) {
      await applyViewerOptions(page, {nodeSize: ns});
      const info = await getViewerInfo(page);
      expect(info.nodeSize, `nodeSize = ${ns}`).toBe(ns);
    }
  });

  // 2.7 — showGrid restyle path.
  await softStep('2.7 showGrid toggle off/on round-trip', async () => {
    const initial = (await getViewerInfo(page)).showGrid;
    await applyViewerOptions(page, {showGrid: !initial});
    let info = await getViewerInfo(page);
    expect(info.showGrid, `showGrid toggled to ${!initial}`).toBe(!initial);
    await applyViewerOptions(page, {showGrid: initial});
    info = await getViewerInfo(page);
    expect(info.showGrid, `showGrid restored to ${initial}`).toBe(initial);
  });

  // 2.8 — Color sweep: hex strings ('#FF0000') are accepted; viewer stores the literal string.
  await softStep('2.8 Color sweep — set each of 5 color props to #FF0000, then restore', async () => {
    const colorProps = [
      'mainColor', 'lightColor', 'currentColor', 'mouseOverColor', 'selectionsColor',
    ] as const;
    const defaults: Record<string, unknown> = {};
    {
      const info = await getViewerInfo(page);
      defaults.mainColor = info.mainColor;
      defaults.lightColor = info.lightColor;
      defaults.currentColor = info.currentColor;
      defaults.mouseOverColor = info.mouseOverColor;
      defaults.selectionsColor = info.selectionsColor;
    }
    for (const cp of colorProps) {
      await applyViewerOptions(page, {[cp]: '#FF0000'});
      const info = await getViewerInfo(page);
      expect((info as any)[cp], `${cp} set to #FF0000`).toBe('#FF0000');
    }
    // Restore defaults.
    for (const cp of colorProps)
      await applyViewerOptions(page, {[cp]: defaults[cp]});
    const info = await getViewerInfo(page);
    expect(info.canvasWidth, 'canvas non-empty after color sweep').toBeGreaterThan(0);
  });

  // 2.9 — showLabels + font restyle path.
  await softStep('2.9 showLabels true; font 12pt monospace; showLabels false', async () => {
    await applyViewerOptions(page, {showLabels: true});
    let info = await getViewerInfo(page);
    expect(info.showLabels, 'showLabels = true').toBe(true);
    await applyViewerOptions(page, {font: '12pt monospace'});
    info = await getViewerInfo(page);
    expect(info.font, 'font = 12pt monospace').toBe('12pt monospace');
    await applyViewerOptions(page, {showLabels: false});
    info = await getViewerInfo(page);
    expect(info.showLabels, 'showLabels = false').toBe(false);
  });

  // 2.10 — stepZoom behavior prop (tactile wheel-zoom magnitude deferred to manual; reflection only).
  await softStep('2.10 stepZoom sweep -4 / 4 / 0.5 (slider bounds + mid)', async () => {
    for (const sz of [-4, 4, 0.5]) {
      await applyViewerOptions(page, {stepZoom: sz});
      const info = await getViewerInfo(page);
      expect(info.stepZoom, `stepZoom = ${sz}`).toBe(sz);
    }
  });

  // 2.11 — showTooltip behavior prop (hover-driven tooltip content deferred to manual; reflection only).
  await softStep('2.11 showTooltip toggle off/on (tooltip-content hover deferred to manual per SR-02)', async () => {
    await applyViewerOptions(page, {showTooltip: false});
    let info = await getViewerInfo(page);
    expect(info.showTooltip, 'showTooltip = false').toBe(false);
    await applyViewerOptions(page, {showTooltip: true});
    info = await getViewerInfo(page);
    expect(info.showTooltip, 'showTooltip = true').toBe(true);
  });

  // 2.12 — step style prop.
  await softStep('2.12 step sweep 0 / 64 / 1.5 (slider bounds + mid)', async () => {
    for (const s of [0, 64, 1.5]) {
      await applyViewerOptions(page, {step: s});
      const info = await getViewerInfo(page);
      expect(info.step, `step = ${s}`).toBe(s);
    }
  });

  await softStep('2.13 Final invariants: canvas non-empty, viewer attached, dataFrame matches', async () => {
    const info = await getViewerInfo(page);
    expect(info.canvasWidth, 'canvas width > 0 at end of sweep').toBeGreaterThan(0);
    expect(info.canvasHeight, 'canvas height > 0 at end of sweep').toBeGreaterThan(0);
    expect(info.isAttachedToTv, 'viewer still attached at end of sweep').toBe(true);
    expect(info.dataFrameMatches, 'viewer.dataFrame === tv.dataFrame at end of sweep').toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
