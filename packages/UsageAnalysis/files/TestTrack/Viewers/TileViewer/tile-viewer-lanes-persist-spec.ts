/* ---
realizes: [tileviewer.cp.lanes-cells-and-layout-persist, tileviewer.int.lane-drag-writes-dataframe-cell, tileviewer.int.lanes-rebuild-vs-restyle-scope, tileviewer.int.column-rename-rewrites-sketch-state]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// Scenarios 1-3 of the lanes-persist scenario. Scenario 4 (layout save / re-apply, project
// save / reopen) is tile-viewer-lanes-persist-server-spec.ts, the section's server sibling.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const ROOT = '[name="viewer-Tile-Viewer"]';

test('Tile Viewer — lanes ladder and tile-content mirroring', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'tile-viewer', 'Tile-Viewer', 10000, 'Tile Viewer');
  await page.locator(`${ROOT} .d4-tile-viewer-form`).nth(5).waitFor({timeout: 10000});

  await softStep('Setup: auto-generated entry state (autoGenerate true, formDesigned false)', async () => {
    const s = await page.evaluate(() => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      return {autoGenerate: viewer.props.autoGenerate, formDesigned: viewer.props.sketchState['formDesigned']};
    });
    expect(s.autoGenerate).toBe(true);
    expect(s.formDesigned).toBe(false);
  });

  let filterProbe: {caucRow: number; blackRow: number; usubCauc: string; usubBlack: string; blackCount: number} =
    {caucRow: -1, blackRow: -1, usubCauc: '', usubBlack: '', blackCount: 0};

  await softStep('Scenario 1 Step 2: single-lane baseline carries .d4-tile-viewer-lane-single', async () => {
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const single = root.querySelectorAll('.d4-tile-viewer-lane-single').length;
      const visibleHeaders = Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header'))
        .filter((h) => getComputedStyle(h as HTMLElement).display !== 'none').length;
      return {
        lanes: root.querySelectorAll('.d4-tile-viewer-lane').length,
        single,
        lanesColumnName: grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer').props.lanesColumnName,
        visibleHeaders,
      };
    });
    expect(r.lanesColumnName).toBeNull();
    expect(r.lanes).toBe(1);
    expect(r.single).toBe(1);
    expect(r.visibleHeaders).toBe(0);
  });

  await softStep('Scenario 1 Step 3: lanesColumnName=RACE yields four .d4-tile-viewer-lane-multi lanes in category order', async () => {
    const r = await page.evaluate(async () => {
      const w = window as any;
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const raceCats = grok.shell.tv.dataFrame.col('RACE').categories;
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { viewer.props.lanesColumnName = 'RACE'; }, 1200);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      return {
        lanes: root.querySelectorAll('.d4-tile-viewer-lane').length,
        multi: root.querySelectorAll('.d4-tile-viewer-lane-multi').length,
        single: root.querySelectorAll('.d4-tile-viewer-lane-single').length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
        raceCats,
      };
    });
    expect(r.lanes).toBe(4);
    expect(r.multi).toBe(4);
    expect(r.single).toBe(0);
    expect(r.headers).toEqual(['Asian', 'Black', 'Caucasian', 'Other']);
    expect(r.headers).toEqual(r.raceCats);
  });

  await softStep('Scenario 1 Step 4: explicit lanes list ["Black","Asian"] renders exactly those two lanes in order', async () => {

    const r = await page.evaluate(async () => {
      const w = window as any;
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { viewer.props.lanes = ['Black', 'Asian']; }, 1200);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      return {
        lanes: root.querySelectorAll('.d4-tile-viewer-lane').length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
      };
    });
    expect(r.headers).toEqual(['Black', 'Asian']);
    expect(r.lanes).toBe(2);
  });

  await softStep('Scenario 1 Step 5: filter RACE to Black keeps both listed lanes; Asian lane holds zero tiles', async () => {

    const ids = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE'); const usub = df.col('USUBJID');
      let caucRow = -1; let blackRow = -1;
      for (let i = 0; i < df.rowCount; i++) {
        if (caucRow < 0 && race.get(i) === 'Caucasian') caucRow = i;
        if (blackRow < 0 && race.get(i) === 'Black') blackRow = i;
        if (caucRow >= 0 && blackRow >= 0) break;
      }
      return {caucRow, blackRow, usubCauc: usub.get(caucRow), usubBlack: usub.get(blackRow), rowCount: df.rowCount};
    });
    filterProbe = {...filterProbe, caucRow: ids.caucRow, blackRow: ids.blackRow, usubCauc: ids.usubCauc, usubBlack: ids.usubBlack};

    expect(filterProbe.caucRow).toBeGreaterThanOrEqual(0);
    expect(filterProbe.usubCauc).toBeTruthy();

    const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', ['Black']);
    filterProbe.blackCount = filteredCount;

    const r = await page.evaluate(() => {
      const w = window as any;
      const tv = grok.shell.tv;
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const read = () => {
        const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
        const perLane = lanes.map((l) => ({
          header: l.querySelector('.d4-tile-viewer-lane-header')?.textContent,
          tiles: l.querySelectorAll('.d4-tile-viewer-form').length,
        }));
        return {showEmptyLanes: viewer.props.showEmptyLanes, laneCount: lanes.length, perLane};
      };
      return w.__poll(read, (x: any) => x.perLane.some((l: any) => l.header === 'Asian' && l.tiles === 0) &&
        x.perLane.some((l: any) => l.header === 'Black' && l.tiles > 0), 1500, 50);
    });
    expect(r.showEmptyLanes).toBe(true);
    expect(r.laneCount).toBe(2);
    const black = r.perLane.find((l: any) => l.header === 'Black')!;
    const asian = r.perLane.find((l: any) => l.header === 'Asian')!;
    expect(black.tiles).toBeGreaterThan(0);
    expect(asian.tiles).toBe(0);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(ids.rowCount);
  });

  await softStep('Scenario 1 Step 6: Show Empty Lanes off drops the empty Asian lane; on restores it in order (GROK-20096)', async () => {

    await v.openViewerGear(page, 'Tile Viewer');
    await v.ensurePropertyCategory(page, 'Tile Viewer', 'data', 'show-empty-lanes');
    await v.setPropertyGridCheckbox(page, 'show-empty-lanes', false, 'data');

    await page.waitForFunction(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]');
      return !!root && root.querySelectorAll('.d4-tile-viewer-lane').length === 1;
    }, null, {timeout: 15_000});
    const off = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      return {
        laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
      };
    });
    expect(off.laneCount).toBe(1);
    expect(off.headers).toEqual(['Black']);

    await v.setPropertyGridCheckbox(page, 'show-empty-lanes', true, 'data');

    await page.waitForFunction(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]');
      return !!root && root.querySelectorAll('.d4-tile-viewer-lane').length === 2;
    }, null, {timeout: 15_000});
    const on = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      return {
        laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
      };
    });
    expect(on.laneCount).toBe(2);
    expect(on.headers).toEqual(['Black', 'Asian']);
  });

  await softStep('Scenario 1 Step 7: clear lanesColumnName returns to a single lane holding the FILTERED (Black) row set', async () => {

    const r = await page.evaluate(async (probe) => {
      const w = window as any;
      const tv = grok.shell.tv;
      const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => {
        viewer.props.lanesColumnName = null;
        viewer.props.lanes = null;
      }, 1500);
      const df = tv.dataFrame;
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      const tiles = Array.from(lanes[0].querySelectorAll('.d4-tile-viewer-form'));
      const races = tiles.map((t) => (t.querySelector('input[name="input-RACE"]') as HTMLInputElement)?.value);
      const usubs = tiles.map((t) => (t.querySelector('input[name="input-USUBJID"]') as HTMLInputElement)?.value);
      return {
        laneCount: lanes.length,
        single: root.querySelectorAll('.d4-tile-viewer-lane-single').length,
        tileCount: tiles.length,
        allBlack: races.length > 0 && races.every((x) => x === 'Black'),
        caucAbsent: !usubs.includes(probe.usubCauc),
        trueCount: df.filter.trueCount,
        rowCount: df.rowCount,
        caucBit: df.filter.get(probe.caucRow),
      };
    }, filterProbe);
    try {
      expect(r.laneCount).toBe(1);
      expect(r.single).toBe(1);
      expect(r.tileCount).toBeGreaterThan(0);
      expect(r.allBlack).toBe(true);
      expect(r.caucAbsent).toBe(true);
      expect(r.caucBit).toBe(false);
      expect(r.trueCount).toBe(filterProbe.blackCount);
      expect(r.trueCount).toBeLessThan(r.rowCount);
    } finally {

      await v.resetFilters(page);
    }
  });

  async function laneDropPoint(header: string, tileNth = 0): Promise<{x: number; y: number; usub: string} | null> {
    return page.evaluate(({h, nth}) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      const lane = lanes.find((l) => l.querySelector('.d4-tile-viewer-lane-header')?.textContent === h);
      if (!lane) return null;
      const content = lane.querySelector('.d4-tile-viewer-lane-content')!;
      const cr = content.getBoundingClientRect();
      const tiles = Array.from(lane.querySelectorAll('.d4-tile-viewer-form'));
      const tile = tiles[nth] as HTMLElement | undefined;
      if (!tile) return null;
      const tr = tile.getBoundingClientRect();

      const x = Math.round(cr.x + cr.width / 2);
      const y = Math.round(tr.y + 30);
      const usub = (tile.querySelector('input[name="input-USUBJID"]') as HTMLInputElement)?.value;

      const hit = document.elementFromPoint(x, y)?.closest('.d4-tile-viewer-lane');
      const hitHeader = hit?.querySelector('.d4-tile-viewer-lane-header')?.textContent;
      if (hitHeader !== h) return null;
      return {x, y, usub};
    }, {h: header, nth: tileNth});
  }

  async function rowOfUsubjid(usub: string): Promise<number> {
    return page.evaluate((u) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('USUBJID');
      for (let i = 0; i < df.rowCount; i++) if (col.get(i) === u) return i;
      return -1;
    }, usub);
  }

  const raceOf = (rowIdx: number) => page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);

  await softStep('Scenario 2 Step 1: lanesColumnName=RACE renders four lanes each with tiles', async () => {
    const r = await page.evaluate(async () => {
      const w = window as any;
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      viewer.props.allowDragBetweenLanes = true;
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { viewer.props.lanesColumnName = 'RACE'; }, 1500);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      return {laneCount: lanes.length, allWithTiles: lanes.every((l) => l.querySelectorAll('.d4-tile-viewer-form').length > 0)};
    });
    expect(r.laneCount).toBe(4);
    expect(r.allWithTiles).toBe(true);
  });

  await softStep('Scenario 2 Step 2: drag a tile from Asian to Black writes the RACE cell and relocates the tile', async () => {

    const src = await laneDropPoint('Asian', 0);
    const dst = await laneDropPoint('Black', 0);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    expect(rowIdx).toBeGreaterThanOrEqual(0);

    expect(await raceOf(rowIdx)).toBe('Asian');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move((src!.x + dst!.x) / 2, (src!.y + dst!.y) / 2, {steps: 8});
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();

    await page.waitForFunction((i) =>
      grok.shell.tv.dataFrame.col('RACE').get(i) === 'Black', rowIdx, {timeout: 15_000}).catch(() => {});

    const r = await page.evaluate((i) => {
      const df = grok.shell.tv.dataFrame;
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      const findLane = (h: string) => lanes.find((l) => l.querySelector('.d4-tile-viewer-lane-header')?.textContent === h)!;
      const usub = df.col('USUBJID').get(i);
      const inLane = (lane: Element) => Array.from(lane.querySelectorAll('.d4-tile-viewer-form'))
        .some((t) => (t.querySelector('input[name="input-USUBJID"]') as HTMLInputElement)?.value === usub);
      const asianLane = findLane('Asian');
      return {
        race: df.col('RACE').get(i),
        currentRowIdx: df.currentRowIdx,
        inBlack: inLane(findLane('Black')),
        inAsian: inLane(asianLane),
        asianTiles: asianLane.querySelectorAll('.d4-tile-viewer-form').length,
      };
    }, rowIdx);
    expect(r.race).toBe('Black');
    expect(r.currentRowIdx).toBe(rowIdx);
    expect(r.inBlack).toBe(true);

    expect(r.asianTiles).toBeGreaterThan(0);
    expect(r.inAsian).toBe(false);
  });

  await softStep('Scenario 2 Step 3: dragging a tile onto its own lane leaves the RACE cell unchanged', async () => {
    const src = await laneDropPoint('Black', 0);
    const dst = await laneDropPoint('Black', 1);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    expect(await raceOf(rowIdx)).toBe('Black');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();

    // a hold that proves the cell does NOT move: capped at the sleep it replaces
    const after = await v.pollValue(() => raceOf(rowIdx), (x) => x !== 'Black', 800, 50);
    expect(after).toBe('Black');
  });

  await softStep('Scenario 2 Step 4: with allowDragBetweenLanes=false the drag leaves the RACE cell unchanged', async () => {

    await v.openViewerGear(page, 'Tile Viewer');
    await v.ensurePropertyCategory(page, 'Tile Viewer', 'misc', 'allow-drag-between-lanes');
    await v.setPropertyGridCheckbox(page, 'allow-drag-between-lanes', false, 'misc');

    const dragDisabled = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer').props.allowDragBetweenLanes);
    expect(dragDisabled).toBe(false);

    const src = await laneDropPoint('Black', 0);
    const dst = await laneDropPoint('Caucasian', 0);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    expect(await raceOf(rowIdx)).toBe('Black');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();

    const after = await v.pollValue(() => raceOf(rowIdx), (x) => x !== 'Black', 800, 50);
    expect(after).toBe('Black');
  });

  await softStep('Scenario 2 Step 5: restore allowDragBetweenLanes and clear lanesColumnName', async () => {
    await v.setPropertyGridCheckbox(page, 'allow-drag-between-lanes', true, 'misc');
    await page.evaluate(() => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      return (window as any).__settled('viewer:Tile Viewer.onViewerRendered', () => { viewer.props.lanesColumnName = null; }, 1000);
    });
  });

  await softStep('Scenario 3 Step 1: record the row-0 AGE tile value baseline', async () => {
    const r = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.tv.dataFrame;
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { df.currentRowIdx = 0; }, 400);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('AGE', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        tileAge: (tile.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value,
        gridText,
      };
    });
    expect(r.tileAge).toBeTruthy();
    expect(r.tileAge).toBe(r.gridText);
  });

  await softStep('Scenario 3 Step 2: editing the AGE grid cell updates the tile value without reopening Edit Form (GROK-17775)', async () => {

    const r = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.tv.dataFrame;
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { df.set('AGE', 0, 99); }, 1200);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('AGE', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        tileAge: (tile.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value,
        gridText,
      };
    });
    expect(r.tileAge).toBeTruthy();
    expect(r.tileAge).toBe(r.gridText);
    expect(r.tileAge).toBe('99');
  });

  await softStep('Scenario 3 Step 3: renaming AGE→AGE_YRS updates the tile label, value selector, and host; an unrenamed column is untouched (GROK-20207)', async () => {

    const r = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.tv.dataFrame;
      const heightBefore = (() => {
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        const t = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
        return (t.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement)?.value;
      })();
      await w.__settled('viewer:Tile Viewer.onViewerRendered', () => { df.columns.byName('AGE').name = 'AGE_YRS'; }, 1200);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      const label = Array.from(tile.querySelectorAll('.d4-host[name="div-AGE-YRS"] input.d4-sketch-column-name'))
        .map((i) => (i as HTMLInputElement).value)[0];
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('AGE_YRS', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        hostPresent: !!tile.querySelector('.d4-host[name="div-AGE-YRS"]'),
        newValue: (tile.querySelector('input[name="input-AGE-YRS"]') as HTMLInputElement)?.value,
        gridText,
        oldGone: tile.querySelector('input[name="input-AGE"]') === null,
        label,
        heightBefore,
        heightAfter: (tile.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement)?.value,
      };
    });
    expect(r.hostPresent).toBe(true);
    expect(r.label).toBe('AGE_YRS');
    expect(r.newValue).toBeTruthy();
    expect(r.newValue).toBe(r.gridText);
    expect(r.oldGone).toBe(true);

    expect(r.heightBefore).toBeTruthy();
    expect(r.heightAfter).toBeTruthy();
    expect(r.heightAfter).toBe(r.heightBefore);

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return (window as any).__settled('viewer:Tile Viewer.onViewerRendered', () => {
        df.columns.byName('AGE_YRS').name = 'AGE';
        df.set('AGE', 0, 53);
      }, 600);
    });
  });

  await softStep('Scenario 3 Step 4: a promoted calculated float column mirrors the grid cell text, not the raw value (GROK-20376)', async () => {

    const before = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const t = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      return Array.from(t.querySelectorAll('input[name^="input-"]')).map((i) => i.getAttribute('name'));
    });

    const r = await page.evaluate(async () => {
      const w = window as any;
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;

      const settle = (act: () => void, capMs: number) => w.__settled('viewer:Tile Viewer.onViewerRendered', act, capMs);
      await settle(() => df.columns.addNewCalculated('COMPUTED_H', '${HEIGHT} * 1.0'), 1200);

      await settle(() => df.columns.remove('DEMOG'), 800);
      await settle(() => df.columns.remove('SEVERITY'), 1500);
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      const after = Array.from(tile.querySelectorAll('input[name^="input-"]')).map((i) => i.getAttribute('name'));
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('COMPUTED_H', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        after,
        tileComputed: (tile.querySelector('input[name="input-COMPUTED-H"]') as HTMLInputElement)?.value,
        gridText,
        raw: String(df.col('COMPUTED_H').get(0)),
      };
    });

    try {

      expect(before).not.toContain('input-COMPUTED-H');
      expect(r.after).toContain('input-COMPUTED-H');
      expect(r.after).not.toEqual(before);

      expect(r.tileComputed).toBeTruthy();
      expect(r.tileComputed).toBe(r.gridText);
      expect(r.tileComputed).not.toBe(r.raw);
    } finally {

      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        try { df.columns.remove('COMPUTED_H'); } catch (_) {  }
      });
    }
  });

  await v.cleanupShell(page);
  v.finishSpec('Tile Viewer lanes failures');
});
