/* ---
realizes: []
--- */

import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// The server lane of the mixed histogram scenario: the layout save/close/reload round-trip
// is the only step whose subject IS server state. The rest stays in histogram-spec.ts.

test('Histogram tests — layout persistence', async ({page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const viewerAlive = () => page.evaluate(() =>
    !!grok.shell.tv.viewers.find(v => v.type === 'Histogram')
    && !!document.querySelector('[name="viewer-Histogram"]'));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  await softStep('Layout persistence', async () => {
    await v.setViewerProps(page, 'Histogram', [
      {set: {valueColumnName: 'WEIGHT', bins: 15, splitColumnName: 'RACE', splitStack: true}, wait: 500},
    ]);
    const result = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const tv = grok.shell.tv;

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000)); 

      h.close();
      for (let i = 0; i < 20 && Array.from(tv.viewers).some((x: any) => x.type === 'Histogram'); i++)
        await new Promise(res => setTimeout(res, 25));

      const saved = await grok.dapi.layouts.find(layoutId);
      let layoutApplied = false;
      const sub = grok.events.onViewLayoutApplied.subscribe(() => { layoutApplied = true; });
      const t0 = Date.now();
      tv.loadLayout(saved);
      while (Date.now() - t0 < 3000) {
        const restored = Array.from(tv.viewers).find((x: any) => x.type === 'Histogram') as any;
        if (layoutApplied && restored?.props) break;
        await new Promise(res => setTimeout(res, 50));
      }
      sub.unsubscribe();

      const h2 = Array.from(tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const r: any[] = [];
      r.push(h2 != null);
      r.push(h2 ? h2.props.valueColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.bins : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitStack : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('WEIGHT');
    expect(result[2]).toBe(15);
    expect(result[3]).toBe('RACE');
    expect(result[4]).toBe(true);
  });

  v.finishSpec();
});
