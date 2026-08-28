/* ---
realizes: []
--- */

import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

// The server lane of the mixed bar-chart scenario: the Data panel step is the only one whose
// subject IS server state — a layout saved through dapi.layouts, closed, and reloaded, plus the
// rebind onto the spgi-100 table. Everything else about the bar chart is client behaviour and
// lives in bar-chart-spec.ts on the local lane.

test('Bar chart tests — data panel', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await openDatagrok(page);

  await v.installEventWaits(page);

  await softStep('Data panel', async () => {
    await v.closeAllAndWait(page);

    const result = await page.evaluate(async () => {
      const w = window as any;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const df2 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const views = Array.from(grok.shell.views).filter((view: any) => view.type === 'TableView');
      const demogView = views.find((view: any) => view.dataFrame.name !== 'SPGI') as any;
      if (demogView)
        await w.__settled('grok.events.onCurrentViewChanged', () => { grok.shell.v = demogView; }, 500);

      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      await w.__settled('grok.events.onViewerAdded', () => icon.click(), 1000);

      const bc = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      const r: any[] = [];

      for (const src of ['Filtered', 'Selected', 'All']) {
        await w.__settled('viewer:Bar chart.onViewerRendered', () => {
          bc.props.rowSource = src;
        }, 200);
        r.push(bc.props.rowSource);
      }
      bc.props.rowSource = 'All';

      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      await w.__settled('viewer:Bar chart.onViewerRendered', () => { bc.dataFrame = spgi; }, 500);
      r.push(bc.dataFrame.name);

      await w.__settled('viewer:Bar chart.onViewerRendered', () => {
        bc.props.filter = '${CAST Idea ID} < 634835';
      }, 500);
      r.push(bc.props.filter);

      await w.__settled('viewer:Bar chart.onViewerRendered', () => {
        bc.props.colorColumnName = 'Chemical Space Y';
      }, 300);
      r.push(bc.props.colorColumnName);

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;

      await w.__settled('grok.events.onViewerClosed', () => bc.close(), 500);

      const saved = await grok.dapi.layouts.find(layoutId);
      await w.__settled('grok.events.onViewLayoutApplied', () => grok.shell.tv.loadLayout(saved), 3000);

      const bc2 = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
      r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result.slice(0, 3)).toEqual(['Filtered', 'Selected', 'All']);
    expect(result[3]).toBe('SPGI');
    expect(result[4]).toBe('${CAST Idea ID} < 634835');
    expect(result[5]).toBe('Chemical Space Y');
    expect(result[6]).toBe('Chemical Space Y');
    expect(result[7]).toBe('${CAST Idea ID} < 634835');
  });

  await softStep('No page errors', async () => {
    expect(pageErrors).toEqual([]);
    expect(consoleErrors).toEqual([]);
  });

  v.finishSpec();
});
