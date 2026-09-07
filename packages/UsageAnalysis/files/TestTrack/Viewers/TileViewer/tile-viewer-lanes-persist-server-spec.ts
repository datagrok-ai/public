/* ---
realizes: [tileviewer.cp.lanes-cells-and-layout-persist, tileviewer.cp.selection-classes-and-form-editor, tileviewer.cp.scroll-survives-added-viewer]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

// The server lane of the section: the steps whose subject is a layout or a project surviving a
// round-trip through the server — Scenario 4 of tile-viewer-lanes-persist.md and Scenario 3
// Step 5 of tile-viewer-selection-form-editor.md — plus the scroll-survives-added-viewer step
// of tile-viewer.md, which is lane-bound: after a dock resize the lane is rebuilt at scrollTop 0
// and only the authenticated client puts the position back (measured 2026-09-03 on dev).
// The ladders that prove the viewer itself are the local-lane siblings.
// expect.poll's default cadence (100/250/500/1000ms) overshoots a settled UI by up to a
// second per call; this ladder is fast early and still coarse for the long product waits
const POLL = {intervals: [50, 50, 100, 200, 400, 800]};

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const ROOT = '[name="viewer-Tile-Viewer"]';

test('Tile Viewer — layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'tile-viewer', 'Tile-Viewer', 10000, 'Tile Viewer');
  await page.locator(`${ROOT} .d4-tile-viewer-form`).nth(5).waitFor({timeout: 10000});

  const AMBIENT = /Permissions policy violation: compute-pressure/i;
  const consoleErrors: string[] = [];
  const onConsole = (m: any) => { if (m.type() === 'error') consoleErrors.push(m.text()); };
  const onPageError = (e: any) => { const t = String(e); if (!AMBIENT.test(t)) consoleErrors.push(t); };
  page.on('console', onConsole);
  page.on('pageerror', onPageError);
  const productErrors = (from: number): string[] => consoleErrors.slice(from).filter((t) => !AMBIENT.test(t));

  const openEditForm = async (): Promise<void> => {
    await page.locator(`${ROOT} .d4-tile-viewer-form .d4-sketch`).first().focus();
    await page.keyboard.press('ContextMenu');
    await page.locator('.d4-menu-popup[name="viewer"] .d4-menu-item[name="div-Edit-Form..."]').click();
    await page.locator('.grok-view-sketch').waitFor({timeout: 15000});
    await page.waitForFunction(() =>
      document.querySelectorAll('.grok-view-sketch .d4-host[name^="div-"]').length > 0,
      null, {timeout: 15000});
  };

  const designerValueHosts = (): Promise<string[]> => page.evaluate(() =>
    Array.from(document.querySelectorAll('.grok-view-sketch .d4-host[name^="div-"]'))
      .filter((h) => h.querySelector('input[name^="input-"]'))
      .map((h) => (h.getAttribute('name') || '').replace('div-', '')).sort());

  await softStep('Scroll position (tile-viewer): a scrolled lane keeps its position and its rows when another viewer is added', async () => {
    await page.evaluate(() => {
      (window as any).__laneRead = () => {
        const l = document.querySelector('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-content') as HTMLElement;
        const first = l?.querySelector('.d4-tile-viewer-form');
        const val = (n: string) => (first?.querySelector(`input[name="input-${n}"]`) as HTMLInputElement)?.value ?? null;
        return JSON.stringify({scrollTop: l?.scrollTop ?? -1, age: val('AGE'), sex: val('SEX'), weight: val('WEIGHT')});
      };
    });
    const readLane = () => page.evaluate(() => {
      const w = window as any;
      return w.__settledFor(w.__laneRead, 250, 1000, 25).then((s: string) => JSON.parse(s));
    });
    // one round trip per wheel: the previous scrollTop is left in the page, so the read after
    // the wheel is also the arming read for the next one
    const laneScroll = (armOnly: boolean) => page.evaluate((arm: boolean) => {
      const w = window as any;
      const read = () => {
        const l = document.querySelector('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-content') as HTMLElement;
        return {top: l.scrollTop, max: l.scrollHeight - l.clientHeight};
      };
      const done = (s: {top: number; max: number}) => { w.__laneTop = s.top; return s; };
      if (arm) return Promise.resolve(done(read()));
      const was = w.__laneTop;
      return w.__poll(read, (s: any) => s.top !== was, 120, 20).then(done);
    }, armOnly);
    const lane = page.locator(`${ROOT} .d4-tile-viewer-lane-content`).first();
    const box = await lane.boundingBox();
    await page.mouse.move(box!.x + box!.width / 2, box!.y + box!.height / 2);
    // a paced sequence, like a drag, stopping mid-range: a lane parked at its scroll clamp is
    // moved by the browser when docking another viewer shortens it, which is not the re-basing
    // this step is about
    let scrolled = await laneScroll(true);
    const target = Math.min(1500, scrolled.max / 2);
    for (let i = 0; i < 10 && scrolled.top < target; i++) {
      await page.mouse.wheel(0, 300);
      scrolled = await laneScroll(false);
    }
    const before = await readLane();
    // docking the histogram rebuilds the lane at scrollTop 0 and the viewer puts the position
    // back ~800 ms later, so the wait is for the lane read to return to its pre-dock value
    await page.evaluate(async (was) => {
      const w = window as any;
      await w.__settled('grok.events.onViewerAdded', () => grok.shell.tv.addViewer('Histogram'), 2500);
      await w.__poll(w.__laneRead, (s: string) => s === was, 2000, 50);
    }, JSON.stringify(before));
    const after = await readLane();

    expect(before.scrollTop).toBeGreaterThan(0);
    expect(before.age).not.toBeNull();

    expect(after.scrollTop).toBeGreaterThan(0);
    expect(Math.abs(after.scrollTop - before.scrollTop)).toBeLessThanOrEqual(2);
    expect(after.age).toBe(before.age);
    expect(after.sex).toBe(before.sex);
    expect(after.weight).toBe(before.weight);

    await page.evaluate(() => grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram')?.close());
    await expect.poll(() => page.evaluate(() =>
      grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram').length), {timeout: 15_000, ...POLL}).toBe(0);
    // the resize rebuilds the lane a beat later; the next step focuses a tile, so leave the lane
    // rebuilt, at the top, and holding still
    await page.evaluate(() => {
      const w = window as any;
      const lane = () => document.querySelector('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-content') as HTMLElement | null;
      const shape = () => {
        const l = lane();
        return l ? `${l.querySelectorAll('.d4-tile-viewer-form').length}|${l.scrollHeight}|${l.scrollTop}` : '';
      };
      return w.__settledFor(shape, 200, 3000, 25).then(() => {
        const l = lane();
        if (l) l.scrollTop = 0;
        return w.__settledFor(shape, 200, 3000, 25);
      });
    });
  });

  await softStep('Scenario 3 Step 5 (selection-form-editor): designed field set survives a layout save + re-apply', async () => {
    const errBefore = consoleErrors.length;

    await openEditForm();
    const values = await designerValueHosts();
    expect(values.length).toBeGreaterThan(0);
    const target = values[0];
    const host = page.locator(`.grok-view-sketch .d4-host[name="div-${target}"]`)
      .filter({has: page.locator(`input[name="input-${target}"]`)}).first();
    await host.click();
    await page.keyboard.press('Delete');
    await host.waitFor({state: 'hidden', timeout: 2000});
    await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
    await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
    await page.locator(`${ROOT} input[name="input-${target}"]`).first().waitFor({state: 'detached', timeout: 15000});
    expect(await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer').props.autoGenerate)).toBe(false);

    // awaiting layouts.save IS the completion signal, the same call the Save to Gallery menu makes
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-tileviewer-form-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await page.evaluate(() => {
        const tv = grok.shell.tv;
        tv.viewers.find((x: any) => x.type === 'Tile Viewer')?.close();
        tv.addViewer('Grid');
      });
      await expect.poll(() => page.evaluate(() =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Tile Viewer').length), {timeout: 15_000, ...POLL}).toBe(0);

      const r = await page.evaluate(async ({id, t}) => {
        const w = window as any;
        const saved = await w.__findSaved(() => grok.dapi.layouts.find(id), 3000);
        grok.shell.tv.loadLayout(saved);
        const read = () => {
          const root = document.querySelector('[name="viewer-Tile-Viewer"]');
          const tiles = Array.from(root?.querySelectorAll('.d4-tile-viewer-form') ?? []);
          return {
            hasTiles: tiles.length > 0,
            removedAbsent: tiles.every((tl) => !tl.querySelector(`input[name="input-${t}"]`)),
            gridGone: !grok.shell.tv.viewers.find((x: any) => x.type === 'Grid' && x !== grok.shell.tv.grid),
          };
        };
        return w.__poll(read, (x: any) => x.hasTiles && x.gridGone, 10000, 50);
      }, {id: layoutId, t: target});

      expect(r.hasTiles).toBe(true);
      expect(r.gridGone).toBe(true);
      expect(r.removedAbsent).toBe(true);
      expect(productErrors(errBefore)).toEqual([]);
    } finally {
      // detached: nothing later in the spec reads this layout, and awaiting the round trip
      // cost 1-4s of the test
      await page.evaluate((id) => {
        const drop = async () => {
          try { const l = await grok.dapi.layouts.find(id); if (l) await grok.dapi.layouts.delete(l); } catch (_) {  }
        };
        drop();
      }, layoutId).catch(() => {});
    }
  });

  let probeLayoutId = '';
  const probeProject: {name?: string; id?: string} = {};
  try {
    await softStep('Scenario 4 Step 1: configure the peak on a fresh viewer (RACE lanes, explicit list, showEmptyLanes)', async () => {
      await page.evaluate(() => grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer')?.close());
      await expect.poll(() => page.evaluate(() =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Tile Viewer').length), {timeout: 15_000, ...POLL}).toBe(0);
      await v.addViewerByIcon(page, 'tile-viewer', 'Tile-Viewer', 10000, 'Tile Viewer');

      const r = await page.evaluate(() => {
        const w = window as any;
        const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        viewer.props.lanesColumnName = 'RACE';
        viewer.props.lanes = ['Black', 'Asian', 'Caucasian'];
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        const read = () => ({
          autoGenerate: viewer.props.autoGenerate,
          laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
          headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
          showEmptyLanes: viewer.props.showEmptyLanes,
        });
        return w.__poll(read, (x: any) => x.laneCount === 3, 2000, 25);
      });
      expect(r.autoGenerate).toBe(true);
      expect(r.laneCount).toBe(3);
      expect(r.headers).toEqual(['Black', 'Asian', 'Caucasian']);

      expect(r.showEmptyLanes).toBe(true);
    });

    await softStep('Scenario 4 Step 2: save the current view layout and confirm it reads back', async () => {
      probeLayoutId = await page.evaluate(async () => {
        const w = window as any;
        const layout = grok.shell.tv.saveLayout();
        layout.name = 'zz-tileviewer-lanes-' + Date.now();
        const id = String((await grok.dapi.layouts.save(layout)).id);
        const back = await w.__findSaved(() => grok.dapi.layouts.find(id), 3000);
        return back ? id : '';
      });
      expect(probeLayoutId.length).toBeGreaterThan(0);
    });

    await softStep('Scenario 4 Step 3: modify the view — close the Tile Viewer and add a Grid', async () => {
      await page.evaluate(() => {
        const tv = grok.shell.tv;
        tv.viewers.find((x: any) => x.type === 'Tile Viewer').close();
        tv.addViewer('Grid');
      });

      await page.waitForFunction(() =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Tile Viewer').length === 0,
        null, {timeout: 15_000});
      const r = await page.evaluate(() =>
        ({tileViewers: grok.shell.tv.viewers.filter((x: any) => x.type === 'Tile Viewer').length}));
      expect(r.tileViewers).toBe(0);
    });

    await softStep('Scenario 4 Step 4: re-applying the saved layout restores the lanes with a clean console (GROK-18230)', async () => {
      const errBefore = consoleErrors.length;
      const isLayoutCrash = (t: string) => /method not found|aPa/i.test(t);

      await page.evaluate(async (id) => {
        const w = window as any;
        const tv = grok.shell.tv;
        const saved = await grok.dapi.layouts.find(id);
        await w.__settled('grok.events.onViewLayoutApplied', () => tv.loadLayout(saved), 3000);
      }, probeLayoutId);

      await page.waitForFunction(() => {
        const root = document.querySelector('[name="viewer-Tile-Viewer"]');
        if (!root) return false;
        const h = Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((x) => x.textContent);
        return h.length === 3;
      }, null, {timeout: 20_000}).catch(() => {});
      const r = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const root = document.querySelector('[name="viewer-Tile-Viewer"]');
        const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        return {
          tilePresent: !!root,
          lanesColumnName: viewer?.props.lanesColumnName ?? null,
          headers: root ? Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent) : [],
        };
      });
      expect(r.tilePresent).toBe(true);
      expect(r.lanesColumnName).toBe('RACE');
      expect(r.headers).toEqual(['Black', 'Asian', 'Caucasian']);
      expect(productErrors(errBefore).filter(isLayoutCrash)).toEqual([]);
    });

    await softStep('Scenario 4 Step 6: save the project, reopen, and spot-check a restored tile value', async () => {

      probeProject.name = `TileLanesPersist_${Date.now()}`;
      const saved = await saveProjectViaApi(page, probeProject.name);
      probeProject.id = saved.projectId;

      await v.closeAllAndWait(page);

      await page.evaluate(async (id) => {
        const p = await grok.dapi.projects.find(id);
        await p.open();
      }, saved.projectId);
      await page.locator(ROOT).waitFor({timeout: 30000});

      // a reopened project has no current row, and with lanes the first tile in the DOM is the
      // first of the "Black" lane rather than row 0, so the row is pinned and its tile is the
      // one compared
      await page.evaluate(() => { grok.shell.tv.dataFrame.currentRowIdx = 0; });
      await page.waitForFunction(() => {
        const root = document.querySelector('[name="viewer-Tile-Viewer"]');
        const tile = root?.querySelector('.d4-tile-viewer-form.d4-current');
        const inp = tile?.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement | null;
        return !!inp && !!inp.value;
      }, null, {timeout: 20_000}).catch(() => {});

      const r = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        const df = tv.dataFrame;
        const tile = root.querySelector('.d4-tile-viewer-form.d4-current')!;
        const heightInput = tile?.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement | null;
        const idx = df.currentRowIdx;
        let gridText: string | null = null;
        try { gridText = tv.grid.cell('HEIGHT', idx).cell.valueString; } catch (_) { gridText = null; }
        return {
          lanesColumnName: viewer?.props.lanesColumnName ?? null,
          headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
          tileHeight: heightInput?.value,
          gridText,
        };
      });
      expect(r.lanesColumnName).toBe('RACE');
      expect(r.headers).toEqual(['Black', 'Asian', 'Caucasian']);
      expect(r.tileHeight).toBeTruthy();
      expect(r.tileHeight).toBe(r.gridText);
    });
  } finally {

    if (probeLayoutId) {
      await page.evaluate((id) => {
        const drop = async () => {
          try { const l = await grok.dapi.layouts.find(id); if (l) await grok.dapi.layouts.delete(l); } catch (_) {  }
        };
        drop();
      }, probeLayoutId).catch(() => {});
    }
    if (probeProject.id)
      await deleteProjectWithCleanup(page, {projectId: probeProject.id});
  }

  page.off('console', onConsole);
  page.off('pageerror', onPageError);
  await v.cleanupShell(page);
  v.finishSpec('Tile Viewer persistence failures');
});
