/* ---
realizes: [correlationplot.cp.matrix-values-scope-persist]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

// The server lane of the mixed matrix-values scenario: a layout saved through dapi.layouts and a
// project saved, closed and reopened. The matrix values themselves are proven in the local lane
// (correlation-plot-matrix-values-spec.ts); here the configured viewer is set up through the API.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;

const CP_STATE = {
  correlationType: 'Spearman',
  showPearsonR: false,
  xColumnNames: ['AGE', 'HEIGHT', 'WEIGHT'],
  yColumnNames: ['AGE', 'HEIGHT'],
  filter: '',
  rowSource: 'Filtered',
};

test('Correlation Plot — Layout and Project Persistence', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);

  try {
    await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    await page.evaluate((state) => {
      const cp = grok.shell.tv.addViewer('Correlation plot');
      for (const k of Object.keys(state)) cp.props[k] = (state as any)[k];
    }, CP_STATE);
    await page.locator('[name="viewer-Correlation-plot"]').waitFor({timeout: 10000});

    const readBack = () => page.evaluate(() => {
      const cp = grok.shell.tv?.viewers?.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv?.dataFrame;
      const ready = !!cp && !!df && df.rowCount > 0;
      return {
        present: !!cp,
        ready,
        type: cp?.props.correlationType ?? null,
        showR: cp?.props.showPearsonR ?? null,
        xCols: cp?.props.xColumnNames.slice() ?? null,
        yCols: cp?.props.yColumnNames.slice() ?? null,
        filter: cp?.props.filter ?? null,
        spot: ready ? cp.getCorrelation(df.col('AGE'), df.col('HEIGHT')) : null,
        ref: ready ? DG.Stats.fromColumn(df.col('AGE')).spearmanCorr(df.col('HEIGHT')) : null,
      };
    });

    await softStep('Scenario 5 Step 5 — saved layout restores full config + value spot-check', async () => {
      const viewerSetBefore: string[] = await page.evaluate(() => grok.shell.tv.viewers.map((x: any) => x.type).sort());
      const layoutId: string = await page.evaluate(async () => {
        const layout = grok.shell.tv.saveLayout();
        await grok.dapi.layouts.save(layout);
        return layout.id;
      });
      try {
        await page.evaluate(async () => {
          const w = window as any;
          grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').close();
          grok.shell.tv.addViewer('Scatter plot');
          await w.__poll(() => {
            const types = grok.shell.tv.viewers.map((x: any) => x.type);
            return types.includes('Scatter plot') && !types.includes('Correlation plot');
          }, (ok: boolean) => ok, 1000, 50);
        });
        await page.evaluate(async (id) => {
          const w = window as any;
          const saved = await w.__findSaved(() => grok.dapi.layouts.find(id), 1000);
          grok.shell.tv.loadLayout(saved);
        }, layoutId);
        const r = await v.pollValue(readBack, (x) => x.ready, 3000, 50);
        const viewerSetAfter: string[] = await page.evaluate(() => grok.shell.tv.viewers.map((x: any) => x.type).sort());
        console.log(`[S5] before=${JSON.stringify(viewerSetBefore)} after=${JSON.stringify(viewerSetAfter)} spot=${r.spot} spearmanRef=${r.ref}`);

        expect(r.present).toBe(true);
        expect(viewerSetAfter).toEqual(viewerSetBefore);
        expect(Number.isFinite(r.spot)).toBe(true);
        expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);

        expect(r.type).toBe('Spearman');
        expect(r.showR).toBe(false);
        expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
        expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
        expect(r.filter === '' || r.filter == null).toBe(true);
      } finally {
        // detached, like deleteProjectWithCleanup: nothing reads the layout again and the
        // worker fixture drains __pendingDeletes before it closes the page
        await page.evaluate((id) => {
          const w = window as any;
          w.__pendingDeletes = w.__pendingDeletes ?? [];
          w.__pendingDeletes.push((async () => {
            try {
              const saved = await grok.dapi.layouts.find(id);
              if (saved) await grok.dapi.layouts.delete(saved);
            } catch (_) {}
          })());
        }, layoutId);
      }
    });

    await softStep('Scenario 6 Step 4/5 — reopened project restores config + value spot-check', async () => {
      const projectName = 'zz-cp-matrix-persist-' + Date.now();
      let savedProjectId: string | null = null;
      try {
        const saved = await saveProjectViaApi(page, projectName);
        savedProjectId = saved.projectId;

        await v.closeAllAndWait(page);
        await page.evaluate(async (id) => {
          const proj = await grok.dapi.projects.find(id);
          await proj.open();
        }, savedProjectId);
        const r = await v.pollValue(readBack, (x) => x.ready, 3000, 50);
        console.log(`[S6] present=${r.present} spot=${r.spot} spearmanRef=${r.ref}`);

        expect(r.present).toBe(true);
        expect(Number.isFinite(r.spot)).toBe(true);
        expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);

        expect(r.type).toBe('Spearman');
        expect(r.showR).toBe(false);
        expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
        expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
      } finally {
        if (savedProjectId)
          await deleteProjectWithCleanup(page, {projectId: savedProjectId});
      }
    });
  } finally {
    await v.closeAllAndWait(page);
  }

  v.finishSpec();
});
