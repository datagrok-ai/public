/* ---
realizes: [tileviewer.cp.lanes-cells-and-layout-persist, tileviewer.int.lane-drag-writes-dataframe-cell, tileviewer.int.lanes-rebuild-vs-restyle-scope, tileviewer.int.column-rename-rewrites-sketch-state]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

test('Tile Viewer — lanes ladder, tile-content mirroring, and persistence', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  // ---- Setup ----
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await v.addViewerByIcon(page, 'tile-viewer', 'Tile-Viewer', 10000, 'Tile Viewer');

  await softStep('Setup: auto-generated entry state (autoGenerate true, formDesigned false)', async () => {
    const s = await page.evaluate(() => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      return {autoGenerate: viewer.props.autoGenerate, formDesigned: viewer.props.sketchState['formDesigned']};
    });
    expect(s.autoGenerate).toBe(true);
    expect(s.formDesigned).toBe(false);
  });

  // Row identity captured before the RACE filter is applied (Step 5), read back in
  // Step 7 to prove the FILTERED row set — not just a non-empty window — is shown.
  let filterProbe: {caucRow: number; blackRow: number; usubCauc: string; usubBlack: string; blackCount: number} =
    {caucRow: -1, blackRow: -1, usubCauc: '', usubBlack: '', blackCount: 0};

  // ===================================================================
  // Scenario 1: Lanes ladder — structure accumulates through the steps
  // ===================================================================

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
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const raceCats = grok.shell.tv.dataFrame.col('RACE').categories;
      viewer.props.lanesColumnName = 'RACE';
      await new Promise((res) => setTimeout(res, 1200));
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
    // The explicit lane list has no UI surface (no prop-grid row, no menu item);
    // set via JS API and asserted through the lane DOM order.
    const r = await page.evaluate(async () => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      viewer.props.lanes = ['Black', 'Asian'];
      await new Promise((res) => setTimeout(res, 1200));
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
    // Capture a Caucasian and a Black row identity BEFORE filtering — Step 7 uses
    // them to prove the single lane holds exactly the filtered (Black) row set.
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

    // Sentinel guard: an unfound Caucasian row leaves the -1/'' seed, making Step 7's
    // caucAbsent/caucBit pass vacuously. Assert the capture is real before any comparison.
    expect(filterProbe.caucRow).toBeGreaterThanOrEqual(0);
    expect(filterProbe.usubCauc).toBeTruthy();

    // RACE categorical filter via the section-sanctioned Filter Panel helper.
    const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', ['Black']);
    filterProbe.blackCount = filteredCount;

    const r = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      const perLane = lanes.map((l) => ({
        header: l.querySelector('.d4-tile-viewer-lane-header')?.textContent,
        tiles: l.querySelectorAll('.d4-tile-viewer-form').length,
      }));
      return {showEmptyLanes: viewer.props.showEmptyLanes, laneCount: lanes.length, perLane};
    });
    expect(r.showEmptyLanes).toBe(true);
    expect(r.laneCount).toBe(2);
    const black = r.perLane.find((l: any) => l.header === 'Black')!;
    const asian = r.perLane.find((l: any) => l.header === 'Asian')!;
    expect(black.tiles).toBeGreaterThan(0);
    expect(asian.tiles).toBe(0);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(ids.rowCount);       // the filter is actually active
  });

  await softStep('Scenario 1 Step 6: Show Empty Lanes off drops the empty Asian lane; on restores it in order (GROK-20096)', async () => {
    // Driven through the REAL property-grid checkbox, not a props write.
    await v.openViewerGear(page, 'Tile Viewer');
    await v.ensurePropertyCategory(page, 'Tile Viewer', 'data', 'show-empty-lanes');
    await v.setPropertyGridCheckbox(page, 'show-empty-lanes', false, 'data');
    const off = await page.evaluate(async () => {
      await new Promise((res) => setTimeout(res, 900));
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      return {
        laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
      };
    });
    expect(off.laneCount).toBe(1);
    expect(off.headers).toEqual(['Black']);

    await v.setPropertyGridCheckbox(page, 'show-empty-lanes', true, 'data');
    const on = await page.evaluate(async () => {
      await new Promise((res) => setTimeout(res, 900));
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
    // The RACE=Black filter is still applied: assert the single lane holds only filtered
    // rows (allBlack, Caucasian absent), cross-checked against df.filter; removed after as cleanup.
    const r = await page.evaluate(async (probe) => {
      const tv = grok.shell.tv;
      const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      viewer.props.lanesColumnName = null;
      viewer.props.lanes = null;
      await new Promise((res) => setTimeout(res, 1500));
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
      expect(r.allBlack).toBe(true);                        // filtered row set shown in the single lane
      expect(r.caucAbsent).toBe(true);                      // the filtered-out Caucasian row is absent
      expect(r.caucBit).toBe(false);
      expect(r.trueCount).toBe(filterProbe.blackCount);     // lane content agrees with the dataframe filter
      expect(r.trueCount).toBeLessThan(r.rowCount);
    } finally {
      // Drop the RACE filter regardless of assert outcome — a failing expect above
      // would otherwise abort the callback and leak the filter into Scenario 2.
      await v.resetFilters(page);
    }
  });

  // ===================================================================
  // Scenario 2: Drag between lanes writes the DataFrame cell
  // Real mouse gestures only — synthetic dispatchEvent does not reach the
  // Dart drag handler. Tiles are clipped to their lane's content box, so a
  // point inside the target lane's hittable x-zone unambiguously targets it.
  // ===================================================================

  // Drop point inside a lane: a tile hittable via elementFromPoint near the lane top.
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
      // x centered in the lane's (narrow) content box; y near the tile's top.
      const x = Math.round(cr.x + cr.width / 2);
      const y = Math.round(tr.y + 30);
      const usub = (tile.querySelector('input[name="input-USUBJID"]') as HTMLInputElement)?.value;
      // confirm the point actually resolves to a tile in THIS lane
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

  await softStep('Scenario 2 Step 1: lanesColumnName=RACE renders four lanes each with tiles', async () => {
    const r = await page.evaluate(async () => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      viewer.props.allowDragBetweenLanes = true;
      viewer.props.lanesColumnName = 'RACE';
      await new Promise((res) => setTimeout(res, 1500));
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      return {laneCount: lanes.length, allWithTiles: lanes.every((l) => l.querySelectorAll('.d4-tile-viewer-form').length > 0)};
    });
    expect(r.laneCount).toBe(4);
    expect(r.allWithTiles).toBe(true);
  });

  await softStep('Scenario 2 Step 2: drag a tile from Asian to Black writes the RACE cell and relocates the tile', async () => {
    // Real mouse drag onto a TILE of the Black lane — the lane header is not a
    // drop target (refdoc tile.md drop section).
    const src = await laneDropPoint('Asian', 0);
    const dst = await laneDropPoint('Black', 0);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    expect(rowIdx).toBeGreaterThanOrEqual(0);

    const raceBefore = await page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);
    expect(raceBefore).toBe('Asian');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move((src!.x + dst!.x) / 2, (src!.y + dst!.y) / 2, {steps: 8});
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();
    await page.waitForTimeout(1000);

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
    // Guard the negative half: the Asian lane must still render tiles, else the
    // dragged row's absence from it would be vacuously true on an empty lane.
    expect(r.asianTiles).toBeGreaterThan(0);
    expect(r.inAsian).toBe(false);
  });

  await softStep('Scenario 2 Step 3: dragging a tile onto its own lane leaves the RACE cell unchanged', async () => {
    const src = await laneDropPoint('Black', 0);
    const dst = await laneDropPoint('Black', 1);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    const before = await page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);
    expect(before).toBe('Black');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();
    await page.waitForTimeout(800);

    const after = await page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);
    expect(after).toBe('Black');
  });

  await softStep('Scenario 2 Step 4: with allowDragBetweenLanes=false the drag leaves the RACE cell unchanged', async () => {
    // Driven through the REAL property-grid checkbox under the Misc category.
    await v.openViewerGear(page, 'Tile Viewer');
    await v.ensurePropertyCategory(page, 'Tile Viewer', 'misc', 'allow-drag-between-lanes');
    await v.setPropertyGridCheckbox(page, 'allow-drag-between-lanes', false, 'misc');

    // Read back that the checkbox gesture drove the property to false — written by the
    // gesture, not the test, so it's not a prop echo and the "drag has no effect" assert is non-vacuous.
    const dragDisabled = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer').props.allowDragBetweenLanes);
    expect(dragDisabled).toBe(false);

    const src = await laneDropPoint('Black', 0);
    const dst = await laneDropPoint('Caucasian', 0);
    expect(src).not.toBeNull();
    expect(dst).not.toBeNull();
    const rowIdx = await rowOfUsubjid(src!.usub);
    const before = await page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);
    expect(before).toBe('Black');

    await page.mouse.move(src!.x, src!.y);
    await page.mouse.down();
    await page.mouse.move(dst!.x, dst!.y, {steps: 8});
    await page.mouse.up();
    await page.waitForTimeout(800);

    const after = await page.evaluate((i) => grok.shell.tv.dataFrame.col('RACE').get(i), rowIdx);
    expect(after).toBe('Black');
  });

  await softStep('Scenario 2 Step 5: restore allowDragBetweenLanes and clear lanesColumnName', async () => {
    await v.setPropertyGridCheckbox(page, 'allow-drag-between-lanes', true, 'misc');
    await page.evaluate(async () => {
      const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      viewer.props.lanesColumnName = null;
      await new Promise((res) => setTimeout(res, 1000));
    });
  });

  // ===================================================================
  // Scenario 3: Tile content mirrors the DataFrame
  // Read channel: input[name="input-<COLUMN>"].value vs grid.cell().cell.valueString.
  // valueString resolves through the SAME Dart cell-format path as col.getString(row),
  // so the grid comparison is that formatting channel renamed, not a stronger one. The
  // GROK-20376 raw-vs-formatted check is carried by expect(tileComputed).not.toBe(raw)
  // in Step 4, not by tile-vs-grid equality.
  // ===================================================================

  await softStep('Scenario 3 Step 1: record the row-0 AGE tile value baseline', async () => {
    const r = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      await new Promise((res) => setTimeout(res, 400));
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('AGE', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        tileAge: (tile.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value,
        gridText,
      };
    });
    expect(r.tileAge).toBeTruthy();          // guard both ends: '' === '' would pass vacuously
    expect(r.tileAge).toBe(r.gridText);
  });

  await softStep('Scenario 3 Step 2: editing the AGE grid cell updates the tile value without reopening Edit Form (GROK-17775)', async () => {
    // The edit is written via df.set (grid cells are canvas, no DOM cell editor —
    // see Automation notes); it reaches the same debounced onValuesChanged path.
    const r = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.set('AGE', 0, 99);
      await new Promise((res) => setTimeout(res, 1200));
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      let gridText: string | null = null;
      try { gridText = grok.shell.tv.grid.cell('AGE', 0).cell.valueString; } catch (_) { gridText = null; }
      return {
        tileAge: (tile.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value,
        gridText,
      };
    });
    expect(r.tileAge).toBeTruthy();          // guard both ends: '' === '' would pass vacuously
    expect(r.tileAge).toBe(r.gridText);      // tile mirrors the grid cell's formatted text (same format channel)
    expect(r.tileAge).toBe('99');
  });

  await softStep('Scenario 3 Step 3: renaming AGE→AGE_YRS updates the tile label, value selector, and host; an unrenamed column is untouched (GROK-20207)', async () => {
    // The rename is actuated via df.columns.byName().name — the same
    // onColumnNameChanged path the grid header rename uses
    // (refdoc: "Column rename follows the DOM").
    const r = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const heightBefore = (() => {
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        const t = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
        return (t.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement)?.value;
      })();
      df.columns.byName('AGE').name = 'AGE_YRS';
      await new Promise((res) => setTimeout(res, 1200));
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
    expect(r.newValue).toBeTruthy();         // guard both ends: '' === '' would pass vacuously
    expect(r.newValue).toBe(r.gridText);
    expect(r.oldGone).toBe(true);
    // The unrenamed HEIGHT tile is unchanged — but only meaningful if HEIGHT is on
    // the tile at all (10-field cap): guard non-empty on BOTH ends before comparing.
    expect(r.heightBefore).toBeTruthy();
    expect(r.heightAfter).toBeTruthy();
    expect(r.heightAfter).toBe(r.heightBefore);

    // rename back and restore the edited value for a clean state
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.columns.byName('AGE_YRS').name = 'AGE';
      df.set('AGE', 0, 53);
      await new Promise((res) => setTimeout(res, 600));
    });
  });

  await softStep('Scenario 3 Step 4: a promoted calculated float column mirrors the grid cell text, not the raw value (GROK-20376)', async () => {
    // The card is capped at 10 fields and demog fills it, so a new column shows only once
    // a slot is freed. Live recon [DOM 2026-08-12]: regeneration drops the LAST-added
    // column, so adding COMPUTED_H then removing two columns frees the two slots to promote
    // it. The name slug turns COMPUTED_H into input-COMPUTED-H.
    const before = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const t = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
      return Array.from(t.querySelectorAll('input[name^="input-"]')).map((i) => i.getAttribute('name'));
    });

    const r = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      df.columns.addNewCalculated('COMPUTED_H', '${HEIGHT} * 1.0');
      await new Promise((res) => setTimeout(res, 1200));
      // free two slots so the frame is at the 10-field cap and COMPUTED_H is promoted
      df.columns.remove('DEMOG');
      await new Promise((res) => setTimeout(res, 800));
      df.columns.remove('SEVERITY');
      await new Promise((res) => setTimeout(res, 1500));
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
      // (a) the field set shown on the tile changed
      expect(before).not.toContain('input-COMPUTED-H');
      expect(r.after).toContain('input-COMPUTED-H');
      expect(r.after).not.toEqual(before);
      // (b) the promoted float column mirrors the GRID cell text, not the raw value
      expect(r.tileComputed).toBeTruthy();
      expect(r.tileComputed).toBe(r.gridText);
      expect(r.tileComputed).not.toBe(r.raw);
    } finally {
      // teardown: drop the calculated column regardless of assert outcome so it
      // cannot leak into Scenario 4 on a red run.
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        try { df.columns.remove('COMPUTED_H'); } catch (_) { /* best effort */ }
        await new Promise((res) => setTimeout(res, 300));
      });
    }
  });

  // ===================================================================
  // Scenario 4: Peak configuration persists through layout + project save
  // ===================================================================

  let probeLayoutId = '';
  const probeProject: {name?: string; id?: string} = {};
  try {
    await softStep('Scenario 4 Step 1: configure the peak (RACE lanes, explicit list, showEmptyLanes)', async () => {
      const r = await page.evaluate(async () => {
        const viewer = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        viewer.props.lanesColumnName = 'RACE';
        viewer.props.lanes = ['Black', 'Asian', 'Caucasian'];
        await new Promise((res) => setTimeout(res, 1200));
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        return {
          laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
          headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
          showEmptyLanes: viewer.props.showEmptyLanes,
        };
      });
      expect(r.laneCount).toBe(3);
      expect(r.headers).toEqual(['Black', 'Asian', 'Caucasian']);
      // Confirm Show Empty Lanes is on by READING, not writing: it's true by construction
      // (Scenario 1 Step 6 left it on); a props write would silently repair a lost setting.
      expect(r.showEmptyLanes).toBe(true);
    });

    await softStep('Scenario 4 Step 2: save the current view layout via View | Layout | Save to Gallery', async () => {
      // Record the applicable-layout ids present before the save so the freshly-
      // saved one can be found deterministically (the gallery save is silent).
      const beforeIds = await page.evaluate(async () =>
        (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)).map((l: any) => String(l.id)));

      // Drive the REAL menu leaf — no programmatic grok.dapi.layouts.save fallback.
      expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery'])).toBe(true);

      let layoutId = '';
      await expect.poll(async () => {
        layoutId = await page.evaluate(async (prev) => {
          const me = String(grok.shell.user.id);
          const ls = await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame);
          const fresh = ls.filter((l: any) => !prev.includes(String(l.id)) &&
            (l.author && l.author.id ? String(l.author.id) === me : true));
          return fresh.length ? String(fresh[fresh.length - 1].id) : '';
        }, beforeIds);
        return layoutId.length;
      }, {timeout: 20_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThan(0);
      probeLayoutId = layoutId;
    });

    await softStep('Scenario 4 Step 3: modify the view — close the Tile Viewer and add a Grid', async () => {
      const r = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        viewer.close();
        await new Promise((res) => setTimeout(res, 500));
        tv.addViewer('Grid');
        await new Promise((res) => setTimeout(res, 800));
        return {tileViewers: tv.viewers.filter((x: any) => x.type === 'Tile Viewer').length};
      });
      expect(r.tileViewers).toBe(0);
    });

    await softStep('Scenario 4 Step 4: re-applying the saved layout restores the lanes with a clean console (GROK-18230)', async () => {
      // GROK-18230/GROK-18037 guard: the layout-apply crash emits "method not found"/"aPa
      // on null" — match only those signatures so ambient console noise doesn't false-FAIL.
      // The APPLY leg is JS-API (a gallery layout has no captured selector); the SAVE leg used the real menu.
      const errors: string[] = [];
      const isLayoutCrash = (t: string) => /method not found|aPa/i.test(t);
      page.on('console', (m) => { if (m.type() === 'error' && isLayoutCrash(m.text())) errors.push(m.text()); });
      page.on('pageerror', (e) => { if (isLayoutCrash(String(e))) errors.push(String(e)); });

      const r = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((res) => setTimeout(res, 3000));
        const root = document.querySelector('[name="viewer-Tile-Viewer"]');
        const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        return {
          tilePresent: !!root,
          lanesColumnName: viewer?.props.lanesColumnName ?? null,
          headers: root ? Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent) : [],
        };
      }, probeLayoutId);
      expect(r.tilePresent).toBe(true);
      expect(r.lanesColumnName).toBe('RACE');
      expect(r.headers).toEqual(['Black', 'Asian', 'Caucasian']);
      expect(errors).toEqual([]);
    });

    await softStep('Scenario 4 Step 6: save the project, reopen, and spot-check a restored tile value', async () => {
      // Set the name BEFORE the save so the finally block can find the entity even
      // if saveProjectViaUI throws after the server committed it (leak window).
      probeProject.name = `TileLanesPersist_${Date.now()}`;
      const saved = await saveProjectViaUI(page, probeProject.name);
      probeProject.id = saved.projectId;

      await page.evaluate(() => grok.shell.closeAll());
      await page.waitForTimeout(1000);

      await page.evaluate(async (id) => {
        const p = await grok.dapi.projects.find(id);
        await p.open();
      }, saved.projectId);
      await page.locator('[name="viewer-Tile-Viewer"]').waitFor({timeout: 30000});
      await page.waitForTimeout(2000);

      const r = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const viewer = tv.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
        const df = tv.dataFrame;
        const tile = root.querySelector('.d4-tile-viewer-form.d4-current') || root.querySelector('.d4-tile-viewer-form')!;
        const heightInput = tile.querySelector('input[name="input-HEIGHT"]') as HTMLInputElement | null;
        const idx = df.currentRowIdx >= 0 ? df.currentRowIdx : 0;
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
      expect(r.tileHeight).toBeTruthy();     // guard both ends: '' === '' would pass vacuously
      expect(r.tileHeight).toBe(r.gridText);
    });
  } finally {
    // Cleanup lives in finally (a softStep swallows failures, so a probe-delete
    // step inside the try can leak on an earlier failure).
    if (probeLayoutId) {
      await page.evaluate(async (id) => {
        try { const l = await grok.dapi.layouts.find(id); if (l) await grok.dapi.layouts.delete(l); } catch (_) { /* best effort */ }
      }, probeLayoutId).catch(() => {});
    }
    if (probeProject.id)
      await deleteProjectWithCleanup(page, {projectId: probeProject.id});
    else if (probeProject.name)
      await page.evaluate(async (name) => {
        const g = (window as any).grok;
        try {
          let p = null;
          try { p = await g.dapi.projects.filter(`name = "${name}"`).first(); } catch (_) { /* index lag */ }
          if (!p) {
            const recent = await g.dapi.projects.list({pageSize: 50}).catch(() => []);
            p = (recent || []).find((x: any) => x && (x.friendlyName === name || x.name === name)) || null;
          }
          if (p) await g.dapi.projects.delete(p);
        } catch (_) { /* best effort */ }
      }, probeProject.name).catch(() => {});
  }

  v.finishSpec('Tile Viewer lanes/persist failures');
});
