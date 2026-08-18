/* ---
realizes: [tileviewer.cp.selection-classes-and-form-editor, tileviewer.int.show-selected-rows-is-viewer-local, tileviewer.int.autogenerate-is-a-state-flag]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

test('Tile Viewer — Row-state rendering and Edit Form designer', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

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

  const AMBIENT = /Permissions policy violation: compute-pressure/i;
  const consoleErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  page.on('pageerror', (e) => { const t = String(e); if (!AMBIENT.test(t)) consoleErrors.push(t); });
  const productErrors = (from: number): string[] =>
    consoleErrors.slice(from).filter((t) => !AMBIENT.test(t));

  await softStep('Setup: the freshly added viewer is auto-generated on both channels', async () => {
    const r = await page.evaluate(() => {
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      return {
        autoGenerate: tileV.props.autoGenerate === true,
        formNotDesigned: tileV.props.sketchState?.['formDesigned'] === false,
      };
    });

    expect(r.autoGenerate).toBe(true);
    expect(r.formNotDesigned).toBe(true);
  });

  const clickTile = async (displayIdx: number, modifiers: string[]): Promise<void> => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const box = await tiles[displayIdx].boundingBox();

    await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      (window as any).__tileClickSettled = new Promise<void>((resolve) => {
        const subs = [
          df.onSelectionChanged.subscribe(() => { done(); }),
          df.onCurrentRowChanged.subscribe(() => { done(); }),
        ];
        const cap = setTimeout(() => done(), 300);
        function done() { clearTimeout(cap); for (const s of subs) s.unsubscribe(); resolve(); }
      });
    });
    for (const m of modifiers) await page.keyboard.down(m);
    await page.mouse.click(box!.x + 15, box!.y + 15);
    for (const m of [...modifiers].reverse()) await page.keyboard.up(m);
    await page.evaluate(() => (window as any).__tileClickSettled);
  };

  const rowOfTile = async (displayIdx: number): Promise<number> => {
    return page.evaluate((idx) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tile = root.querySelectorAll('.d4-tile-viewer-form')[idx];
      const rd = (c: string) => (tile.querySelector(`input[name="input-${c}"]`) as HTMLInputElement)?.value;
      const sex = rd('SEX'); const wt = rd('WEIGHT'); const ht = rd('HEIGHT');
      for (let i = 0; i < df.rowCount; i++) {
        if (df.col('SEX').getString(i) === sex && df.col('WEIGHT').getString(i) === wt &&
            df.col('HEIGHT').getString(i) === ht) return i;
      }
      return -1;
    }, displayIdx);
  };

  let row2 = -1; let row4 = -1; let row6 = -1;
  await softStep('Scenario 1 Step 1: baseline — no current/selected tile, empty selection', async () => {
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.currentRowIdx = -1;
      return {
        current: root.querySelectorAll('.d4-tile-viewer-form.d4-current').length,
        selected: root.querySelectorAll('.d4-tile-viewer-form.d4-selected').length,
        selCount: df.selection.trueCount,
      };
    });
    expect(r.current).toBe(0);
    expect(r.selected).toBe(0);
    expect(r.selCount).toBe(0);

    row2 = await rowOfTile(1);
    row4 = await rowOfTile(3);
    row6 = await rowOfTile(5);
    expect(row2).toBeGreaterThanOrEqual(0);
    expect(row4).toBeGreaterThanOrEqual(0);
    expect(row6).toBeGreaterThanOrEqual(0);
  });

  await softStep('Scenario 1 Step 2: plain click sets current only, selection untouched', async () => {
    await clickTile(1, []);
    const r = await page.evaluate((row) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const otherCurrent = tiles.filter((t, i) => i !== 1 && t.classList.contains('d4-current')).length;
      return {
        tileCurrent: tiles[1].classList.contains('d4-current'),
        currentIdx: df.currentRowIdx,
        expectedRow: row,
        selCount: df.selection.trueCount,
        currentTiles: root.querySelectorAll('.d4-tile-viewer-form.d4-current').length,
        otherCurrent,
      };
    }, row2);
    expect(r.tileCurrent).toBe(true);
    expect(r.currentIdx).toBe(r.expectedRow);
    expect(r.selCount).toBe(0);
    expect(r.currentTiles).toBe(1);
    expect(r.otherCurrent).toBe(0);
  });

  await softStep('Scenario 1 Step 3: Ctrl-click toggles selection in both directions', async () => {
    await clickTile(3, ['Control']);
    const on = await page.evaluate((row) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tile = root.querySelectorAll('.d4-tile-viewer-form')[3];
      return {tileSelected: tile.classList.contains('d4-selected'), dfSelected: df.selection.get(row)};
    }, row4);
    expect(on.tileSelected).toBe(true);
    expect(on.dfSelected).toBe(true);
    await clickTile(3, ['Control']);
    const off = await page.evaluate((row) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tile = root.querySelectorAll('.d4-tile-viewer-form')[3];
      return {tileSelected: tile.classList.contains('d4-selected'), dfSelected: df.selection.get(row)};
    }, row4);
    expect(off.tileSelected).toBe(false);
    expect(off.dfSelected).toBe(false);
  });

  await softStep('Scenario 1 Step 4: Shift-click is additive for one row, not a range', async () => {
    await clickTile(3, ['Control']);
    await clickTile(5, ['Shift']);
    const r = await page.evaluate(({r4, r6}) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));

      const betweenSelectedAny = [2, 3, 4].some((i) => tiles[i]?.classList.contains('d4-selected') && i !== 3);
      return {
        row6Selected: tiles[5].classList.contains('d4-selected'),
        row6Df: df.selection.get(r6),
        row4StillDf: df.selection.get(r4),
        betweenSelectedAny,

        selTrueCount: df.selection.trueCount,
      };
    }, {r4: row4, r6: row6});
    expect(r.row6Selected).toBe(true);
    expect(r.row6Df).toBe(true);
    expect(r.row4StillDf).toBe(true);
    expect(r.betweenSelectedAny).toBe(false);
    expect(r.selTrueCount).toBe(2);
  });

  await softStep('Scenario 1 Step 5: Ctrl+Shift-click clears just that row', async () => {
    await clickTile(5, ['Control', 'Shift']);
    const r = await page.evaluate(({r4, r6}) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        row6Selected: tiles[5].classList.contains('d4-selected'),
        row6Df: df.selection.get(r6),
        row4Selected: tiles[3].classList.contains('d4-selected'),
        row4Df: df.selection.get(r4),
      };
    }, {r4: row4, r6: row6});
    expect(r.row6Selected).toBe(false);
    expect(r.row6Df).toBe(false);
    expect(r.row4Selected).toBe(true);
    expect(r.row4Df).toBe(true);
  });

  await softStep('Scenario 2 Step 1: two tiles selected, selection count is 2', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      df.selection.setAll(false);
      tileV.props.showSelectedRows = true;
      tileV.props.rowSource = 'Filtered';
    });
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);

    await clickTile(0, []);
    await clickTile(1, ['Control']);
    await clickTile(3, ['Control']);
    const s = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      return {
        selCount: df.selection.trueCount,
        selectedTiles: root.querySelectorAll('.d4-tile-viewer-form.d4-selected').length,
      };
    });
    expect(s.selCount).toBe(2);
    expect(s.selectedTiles).toBe(2);
  });

  await softStep('Scenario 2 Step 2: Show Selected Rows off neutralises the highlight only', async () => {
    const baseline = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const selTile = root.querySelector('.d4-tile-viewer-form.d4-selected') as HTMLElement;
      const unselTile = Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .find((t) => !t.classList.contains('d4-selected')) as HTMLElement;
      return {
        bgSel: getComputedStyle(selTile).backgroundColor,
        bgUnsel: getComputedStyle(unselTile).backgroundColor,
        selBefore: grok.shell.tv.dataFrame.selection.trueCount,
      };
    });

    await v.openViewerGear(page, 'Tile-Viewer');
    await v.ensurePropertyCategory(page, 'Tile-Viewer', 'selection', 'show-selected-rows');
    await v.setPropertyGridCheckbox(page, 'show-selected-rows', false, 'selection');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const selTile = () => root.querySelector('.d4-tile-viewer-form.d4-selected') as HTMLElement;
      const unselTile = () => Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .find((t) => !t.classList.contains('d4-selected')) as HTMLElement;
      let hostHasClass = false; let bgNeutralised = false;
      for (let i = 0; i < 25; i++) {
        await new Promise((res) => setTimeout(res, 200));
        const host = root.querySelector('.d4-tile-viewer-lanes-host')!;
        hostHasClass = host.classList.contains('d4-tile-viewer-hide-selected');
        bgNeutralised = getComputedStyle(selTile()).backgroundColor === getComputedStyle(unselTile()).backgroundColor;
        if (hostHasClass && bgNeutralised) break;
      }
      return {
        hostHasClass,
        bgNeutralised,
        selCountAfter: df.selection.trueCount,
        selectedTilesStill: root.querySelectorAll('.d4-tile-viewer-form.d4-selected').length,
      };
    });
    expect(baseline.bgSel).not.toBe(baseline.bgUnsel);
    expect(r.hostHasClass).toBe(true);
    expect(r.bgNeutralised).toBe(true);
    expect(r.selCountAfter).toBe(baseline.selBefore);
    expect(r.selectedTilesStill).toBe(2);
  });

  await softStep('Scenario 2 Step 3: Show Selected Rows on restores the highlight', async () => {

    await v.ensurePropertyCategory(page, 'Tile-Viewer', 'selection', 'show-selected-rows');
    await v.setPropertyGridCheckbox(page, 'show-selected-rows', true, 'selection');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const selTile = () => root.querySelector('.d4-tile-viewer-form.d4-selected') as HTMLElement;
      const unselTile = () => Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .find((t) => !t.classList.contains('d4-selected')) as HTMLElement;
      let hostHasClass = true; let bgDiffers = false;
      for (let i = 0; i < 25; i++) {
        await new Promise((res) => setTimeout(res, 200));
        const host = root.querySelector('.d4-tile-viewer-lanes-host')!;
        hostHasClass = host.classList.contains('d4-tile-viewer-hide-selected');
        bgDiffers = getComputedStyle(selTile()).backgroundColor !== getComputedStyle(unselTile()).backgroundColor;
        if (!hostHasClass && bgDiffers) break;
      }
      return {hostHasClass, bgDiffers};
    });
    expect(r.hostHasClass).toBe(false);
    expect(r.bgDiffers).toBe(true);
  });

  await softStep('Scenario 2 Step 5: rowSource=Selected forces suppression and dims the property cell', async () => {

    await v.openViewerProperties(page, 'Tile-Viewer', '[name="prop-category-selection"]');
    await v.ensurePropertyCategory(page, 'Tile-Viewer', 'selection', 'show-selected-rows');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');

      const propRow = () => document.querySelector('.property-grid tr[name="prop-show-selected-rows"]') as HTMLElement;

      const rowState = () => {
        const el = propRow();
        if (!el) return null;
        const cells = Array.from(el.querySelectorAll('td'));
        const box = el.querySelector('input[type="checkbox"]') as HTMLInputElement;
        return {
          opacity: parseFloat(getComputedStyle(el).opacity),
          rowDisabled: el.hasAttribute('disabled'),
          cellCount: cells.length,
          allCellsClickable: cells.every((c) => getComputedStyle(c).pointerEvents === 'auto'),
          boxDisabled: box ? (box.disabled || box.hasAttribute('disabled')) : null,
          boxChecked: box ? box.checked : null,
        };
      };
      const hostSuppressed = () => root.querySelector('.d4-tile-viewer-lanes-host')!.classList
        .contains('d4-tile-viewer-hide-selected');
      const selBefore = df.selection.trueCount;
      const before = rowState(); 
      tileV.props.showSelectedRows = true; 
      tileV.props.rowSource = 'Selected';
      let hostHasClass = false; let dimmed = rowState();
      for (let i = 0; i < 25; i++) {
        await new Promise((res) => setTimeout(res, 200));
        hostHasClass = hostSuppressed();
        dimmed = rowState();
        if (hostHasClass && dimmed && dimmed.opacity === 0.5) break;
      }

      const selectedIdx: number[] = [];
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selectedIdx.push(i);
      const tileAges = Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name="input-AGE"]'))
        .map((i) => (i as HTMLInputElement).value).sort();
      const selectedAges = selectedIdx.map((i) => df.col('AGE').getString(i)).sort();
      const visibleTiles = root.querySelectorAll('.d4-tile-viewer-form').length;

      tileV.props.rowSource = 'Filtered';
      let hostClassGone = false; let restored = rowState();
      for (let i = 0; i < 25; i++) {
        await new Promise((res) => setTimeout(res, 200));
        hostClassGone = !hostSuppressed();
        restored = rowState();
        if (hostClassGone && restored && restored.opacity === 1) break;
      }
      return {
        hostHasClass, before, dimmed, restored, tileAges, selectedAges, visibleTiles,
        hostClassGone,
        selBefore, selAfter: df.selection.trueCount,
        selectedCount: selectedIdx.length,
      };
    });
    expect(r.hostHasClass).toBe(true);

    expect(r.selectedCount).toBeGreaterThan(0);
    expect(r.visibleTiles).toBe(r.selectedCount);
    expect(r.tileAges).toEqual(r.selectedAges);
    expect(r.hostClassGone).toBe(true);
    expect(r.selAfter).toBe(r.selBefore);
    expect(r.before?.opacity).toBeCloseTo(1, 1);
    expect(r.dimmed?.opacity).toBeCloseTo(0.5, 1);
    expect(r.restored?.opacity).toBeCloseTo(1, 1);

    expect(r.dimmed?.cellCount).toBeGreaterThan(0);
    expect(r.before?.rowDisabled).toBe(false);
    expect(r.dimmed?.rowDisabled).toBe(false);
    expect(r.restored?.rowDisabled).toBe(false);
    expect(r.dimmed?.allCellsClickable).toBe(true);
    expect(r.restored?.allCellsClickable).toBe(true);
    expect(r.before?.boxDisabled).toBe(false);
    expect(r.dimmed?.boxDisabled).toBe(false);
    expect(r.dimmed?.boxChecked).toBe(r.before?.boxChecked);
    expect(r.restored?.boxChecked).toBe(r.before?.boxChecked);
  });

  await page.evaluate(() => {
    const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
    tileV.props.rowSource = 'All';
    tileV.props.showSelectedRows = true;
    grok.shell.tv.dataFrame.selection.setAll(false);
  });
  await v.waitForViewerRendered(page, 'Tile Viewer', 900);

  const openEditForm = async (rootSelector = '[name="viewer-Tile-Viewer"]'): Promise<void> => {
    await page.locator(`${rootSelector} .d4-tile-viewer-form .d4-sketch`).first().focus();
    await page.keyboard.press('ContextMenu');
    await page.locator('.d4-menu-popup[name="viewer"] .d4-menu-item[name="div-Edit-Form..."]').click();
    await page.locator('.grok-view-sketch').waitFor({timeout: 15000});

    await page.waitForFunction(() =>
      document.querySelectorAll('.grok-view-sketch .d4-host[name^="div-"]').length > 0,
      null, {timeout: 15000});
  };

  const designerHosts = async (): Promise<{values: string[]; labels: string[]; total: number}> => {
    return page.evaluate(() => {
      const hosts = Array.from(document.querySelectorAll('.grok-view-sketch .d4-host[name^="div-"]'));
      const nameOf = (h: Element) => (h.getAttribute('name') || '').replace('div-', '');
      const values = hosts.filter((h) => h.querySelector('input[name^="input-"]'));
      const labels = hosts.filter((h) => !values.includes(h) && h.querySelector('input.d4-sketch-column-name'));
      return {values: values.map(nameOf).sort(), labels: labels.map(nameOf).sort(), total: hosts.length};
    });
  };

  await softStep('Scenario 3 Step 2: the designer opens on the correct table', async () => {
    await openEditForm();
    const r = await page.evaluate(() => {

      const allViewers = Array.from(grok.shell.views)
        .flatMap((view: any) => view.viewers ? Array.from(view.viewers) : []);
      const tileV = allViewers.find((x: any) => x && x.type === 'Tile Viewer') as any;
      const demog = grok.shell.tables.find((t: any) => t.name === 'Table');
      return {
        sketchViewOpen: !!document.querySelector('.grok-view-sketch'),
        table: tileV?.props?.sketchState?.['table'] ?? null,
        frameName: demog?.name ?? null,
      };
    });

    expect(r.sketchViewOpen).toBe(true);
    expect(r.frameName).toBeTruthy();
    expect(r.table).toBe(r.frameName);
  });

  await softStep('Scenario 3 Step 2b: the Select columns dialog counter matches the card field count', async () => {

    const openedValues = (await designerHosts()).values;
    const fieldCount = openedValues.length;
    expect(fieldCount).toBeGreaterThan(0);

    await page.locator('[name="button-EDIT"]').click();
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 15000});
    const d = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Select-columns..."]');
      if (!dlg) return null;
      const search = dlg.querySelector('.d4-column-grid input.d4-search-input')
        ?? dlg.querySelector('input.d4-search-input');
      let counter = -1;
      for (const l of Array.from(dlg.querySelectorAll('label'))) {
        const m = (l.textContent || '').match(/(\d+)\s+checked/);
        if (m) { counter = parseInt(m[1], 10); break; }
      }
      return {
        searchPresent: !!search,
        labelAllPresent: !!dlg.querySelector('[name="label-All"]'),
        labelNonePresent: !!dlg.querySelector('[name="label-None"]'),
        counter,
      };
    });

    await page.locator('[name="dialog-Select-columns..."] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Select-columns..."]').waitFor({state: 'detached', timeout: 10000});
    expect(d).not.toBeNull();
    expect(d!.searchPresent).toBe(true);
    expect(d!.labelAllPresent).toBe(true);
    expect(d!.labelNonePresent).toBe(true);
    expect(d!.counter).toBe(fieldCount);

    const stillOpen = await designerHosts();
    expect(stillOpen.values).toEqual(openedValues);
  });

  await softStep('Scenario 3 Step 3: remove a field, apply, transition to designed state', async () => {
    const errBefore = consoleErrors.length;

    const opened = await designerHosts();
    expect(opened.values.length).toBeGreaterThan(0);
    const removed = opened.values[0];
    await page.locator(`.grok-view-sketch .d4-host[name="div-${removed}"]`)
      .filter({has: page.locator(`input[name="input-${removed}"]`)}).first().click();
    await page.keyboard.press('Delete');
    await page.waitForTimeout(300); 
    await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
    await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
    await page.waitForTimeout(700); 

    const r = await page.evaluate((gone: string) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const removedAbsent = tiles.every((t) => !t.querySelector(`input[name="input-${gone}"]`));
      const tile = tiles[0];
      const fieldCols = Array.from(tile.querySelectorAll('input[name^="input-"]'))
        .map((i) => (i.getAttribute('name') || '').replace('input-', ''));
      const valueOf = (c: string) => (tile.querySelector(`input[name="input-${c}"]`) as HTMLInputElement);

      const probes = fieldCols.filter((c) => {
        const input = valueOf(c);
        return input && input.type !== 'checkbox' && df.col(c);
      }).slice(0, 2);
      let boundRow = -1;
      for (let i = 0; i < df.rowCount; i++)
        if (probes.every((c) => df.col(c).getString(i) === valueOf(c).value)) { boundRow = i; break; }
      const mismatches: string[] = [];
      for (const c of fieldCols) {
        const input = valueOf(c);
        if (!input || input.type === 'checkbox') continue; 
        if (df.col(c) && input.value !== df.col(c).getString(boundRow)) mismatches.push(c);
      }
      return {
        removedAbsent,
        autoGenerateFalse: tileV.props.autoGenerate === false,
        formDesignedTrue: tileV.props.sketchState?.['formDesigned'] === true,
        boundRow,
        probeCount: probes.length,
        mismatches,
      };
    }, removed);
    expect(r.removedAbsent).toBe(true);
    expect(r.probeCount).toBeGreaterThan(0);
    expect(r.boundRow).toBeGreaterThanOrEqual(0);
    expect(r.mismatches).toEqual([]);
    expect(r.autoGenerateFalse).toBe(true);
    expect(r.formDesignedTrue).toBe(true);
    expect(productErrors(errBefore)).toEqual([]);
  });

  await softStep('Scenario 3 Step 4: delete a label host, then RESET reverts to the opening state', async () => {
    const errBefore = consoleErrors.length;
    await openEditForm();
    const opened = await designerHosts();

    const appliedRemovals = opened.labels.filter((n) => !opened.values.includes(n));
    expect(appliedRemovals.length).toBeGreaterThan(0);

    const target = opened.labels.find((n) => opened.values.includes(n));
    expect(target).toBeTruthy();
    expect(opened.values.length + opened.labels.length).toBe(opened.total);
    await page.locator(`.grok-view-sketch .d4-host[name="div-${target}"]`)
      .filter({has: page.locator('input.d4-sketch-column-name')})
      .filter({hasNot: page.locator(`input[name="input-${target}"]`)}).first().click();
    await page.keyboard.press('Delete');
    await page.waitForTimeout(300); 
    const afterDelete = await designerHosts();

    await page.locator('[name="button-RESET"]').click();
    await page.waitForTimeout(700); 
    const afterReset = await designerHosts();

    await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
    await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
    await page.waitForTimeout(700); 

    const applied = await page.evaluate((t: string) => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const tile = tiles[0];
      const wt = (tile.querySelector('input[name="input-WEIGHT"]') as HTMLInputElement)?.value;
      let boundRow = -1;
      for (let i = 0; i < df.rowCount; i++)
        if (df.col('WEIGHT').getString(i) === wt) { boundRow = i; break; }
      const targetInput = tile.querySelector(`input[name="input-${t}"]`) as HTMLInputElement;
      return {
        labelBack: tiles.every((tl) =>
          !!tl.querySelector(`.d4-host[name="div-${t}"] input.d4-sketch-column-name`)),
        valuePresent: tiles.every((tl) => !!tl.querySelector(`input[name="input-${t}"]`)),
        targetEqualsDisplay: df.col(t) ? targetInput?.value === df.col(t).getString(boundRow) : null,
        boundRow,
        autoGenerateFalse: tileV.props.autoGenerate === false,
        formDesignedTrue: tileV.props.sketchState?.['formDesigned'] === true,
      };
    }, target!);

    expect(afterDelete.labels).toEqual(opened.labels.filter((n) => n !== target));
    expect(afterDelete.values).toEqual(opened.values);
    expect(afterReset.labels).toEqual(opened.labels);
    expect(afterReset.values).toEqual(opened.values);

    for (const n of appliedRemovals) expect(afterReset.values).not.toContain(n);
    expect(applied.labelBack).toBe(true);
    expect(applied.valuePresent).toBe(true);
    expect(applied.boundRow).toBeGreaterThanOrEqual(0);
    expect(applied.targetEqualsDisplay).toBe(true);
    expect(applied.autoGenerateFalse).toBe(true);
    expect(applied.formDesignedTrue).toBe(true);

    expect(productErrors(errBefore)).toEqual([]);
  });

  await softStep('Scenario 3 Step 5: designed field set survives a layout save + re-apply', async () => {
    const errBefore = consoleErrors.length;

    await openEditForm();
    const opened = await designerHosts();
    expect(opened.values.length).toBeGreaterThan(0);
    const target = opened.values[0];
    await page.locator(`.grok-view-sketch .d4-host[name="div-${target}"]`)
      .filter({has: page.locator(`input[name="input-${target}"]`)}).first().click();
    await page.keyboard.press('Delete');
    await page.waitForTimeout(300); 
    await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
    await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});

    await page.locator(`[name="viewer-Tile-Viewer"] input[name="input-${target}"]`)
      .first().waitFor({state: 'detached', timeout: 15000});

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise((res) => setTimeout(res, 1000)); 
      return layout.id;
    });

    await page.evaluate(async () => {
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      if (tileV) tileV.close();
      grok.shell.tv.addViewer('Grid');
      await new Promise((res) => setTimeout(res, 600)); 
    });

    const r = await page.evaluate(async ({id, t}) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);

      for (let i = 0; i < 40; i++) {
        const back = document.querySelector('[name="viewer-Tile-Viewer"]');
        const gridGone = !grok.shell.tv.viewers.find((x: any) => x.type === 'Grid' && x !== grok.shell.tv.grid);
        if (back && back.querySelectorAll('.d4-tile-viewer-form').length > 0 && gridGone) break;
        await new Promise((res) => setTimeout(res, 250));
      }
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const removedAbsent = tiles.every((tl) => !tl.querySelector(`input[name="input-${t}"]`));
      const hasTiles = tiles.length > 0;
      await grok.dapi.layouts.delete(saved);
      return {removedAbsent, hasTiles};
    }, {id: layoutId, t: target});

    expect(r.hasTiles).toBe(true);
    expect(r.removedAbsent).toBe(true);
    expect(productErrors(errBefore)).toEqual([]);
  });

  const AG_COLS = ['COLA', 'COLB', 'COLC', 'COLD', 'COLE', 'COLF'];

  const fixtureFields = async (): Promise<{order: string[]; covered: string[]; uncovered: string[]; tiles: number}> => {
    return page.evaluate(() => {
      const view = Array.from(grok.shell.views).find((x: any) => x.name === 'ag-fixture') as any;
      const root = view.root.querySelector('[name="viewer-Tile-Viewer"]') as Element;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const nameOf = (i: Element) => (i.getAttribute('name') || '').replace('input-', '');
      const order = tiles.length ?
        Array.from(tiles[0].querySelectorAll('input[name^="input-"]')).map(nameOf) : [];
      const anywhere = new Set(Array.from(root
        .querySelectorAll('.d4-tile-viewer-form input[name^="input-"]')).map(nameOf));
      const names = view.dataFrame.columns.names();
      return {
        order,
        covered: names.filter((c: string) => anywhere.has(c)).sort(),
        uncovered: names.filter((c: string) => !anywhere.has(c)).sort(),
        tiles: tiles.length,
      };
    });
  };

  const removeColumnByFieldMenu = async (rootSelector: string, slug: string): Promise<void> => {
    await page.locator(`${rootSelector} .d4-tile-viewer-form input[name="input-${slug}"]`).first()
      .click({button: 'right'});
    await page.locator(`.d4-menu-popup[name="${slug}"] .d4-menu-item[name="div-Remove"]`).click();
    await page.waitForTimeout(300); 
  };
  try {
    await softStep('Scenario 4 Step 2: auto-generated state — delete regenerates without the field', async () => {
      const setup = await page.evaluate(async (cols: string[]) => {
        const t = DG.DataFrame.create(6);
        t.name = 'ag-fixture';
        for (const c of cols) t.columns.addNewInt(c).init((i: number) => i * (cols.indexOf(c) + 1));
        const view = grok.shell.addTableView(t);
        const tileV = view.addViewer('Tile Viewer');

        await new Promise((res) => {
          let done = false;
          const finish = () => { if (!done) { done = true; sub?.unsubscribe(); clearTimeout(cap); res(undefined); } };
          const rootReady = () => {
            const r = view.root.querySelector('[name="viewer-Tile-Viewer"]');
            return !!r && r.querySelectorAll('.d4-tile-viewer-form').length > 0;
          };
          let sub: any = tileV.onViewerRendered.subscribe(() => { if (rootReady()) finish(); });
          const cap = setTimeout(finish, 2000);
          if (rootReady()) finish();
        });
        const root = view.root.querySelector('[name="viewer-Tile-Viewer"]') as Element;
        root.setAttribute('data-ag-fixture', '1');
        const hosts = new Set(Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
          .map((i) => i.getAttribute('name')));
        return {
          renderedBefore: t.columns.names().filter((c: string) => hosts.has(`input-${c}`)),

          autoGenerateTrue: tileV.props.autoGenerate === true,
          formDesignedFalse: tileV.props.sketchState?.['formDesigned'] === false,
        };
      }, AG_COLS);

      const victim = setup.renderedBefore[0];
      await removeColumnByFieldMenu('[data-ag-fixture="1"]', victim);

      const r = await page.evaluate(async (v: string) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === 'ag-fixture') as any;
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = view.root.querySelector('[name="viewer-Tile-Viewer"]') as Element;
        const rendered = () => {
          const hosts = new Set(Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
            .map((i) => i.getAttribute('name')));
          return t.columns.names().filter((c: string) => hosts.has(`input-${c}`));
        };

        let renderedAfter = rendered();
        let survivorsOk = false;
        let labelGone = false;
        for (let i = 0; i < 30; i++) {
          if (!renderedAfter.includes(v)) {
            const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
            const tile = tiles[0];
            survivorsOk = renderedAfter.every((c: string) => {
              const input = tile.querySelector(`input[name="input-${c}"]`) as HTMLInputElement;
              return input && input.value === t.col(c).getString(0);
            });

            labelGone = tiles.every((tl) =>
              !tl.querySelector(`.d4-host[name="div-${v}"] input.d4-sketch-column-name`));
            break;
          }
          await new Promise((res) => setTimeout(res, 400));
          renderedAfter = rendered();
        }
        return {
          columnGoneFromDf: !t.columns.names().includes(v),
          fieldGone: !renderedAfter.includes(v),
          labelGone,
          survivorsOk,
          autoGenerateStillTrue: tileV.props.autoGenerate === true,
        };
      }, victim);

      expect(setup.autoGenerateTrue).toBe(true);
      expect(setup.formDesignedFalse).toBe(true);
      expect(setup.renderedBefore.length).toBe(AG_COLS.length);
      expect(r.columnGoneFromDf).toBe(true);
      expect(r.fieldGone).toBe(true);
      expect(r.labelGone).toBe(true);
      expect(r.survivorsOk).toBe(true);
      expect(r.autoGenerateStillTrue).toBe(true);
    });

    await softStep('Scenario 4 Step 4: designed state — the freed slot is not refilled', async () => {

      const EXTRA_COLS = ['COLG', 'COLH', 'COLI', 'COLJ', 'COLK', 'COLL', 'COLM', 'COLN'];
      await page.evaluate(async (cols: string[]) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === 'ag-fixture') as any;
        grok.shell.v = view;
        const t = view.dataFrame;
        for (const c of cols) t.columns.addNewInt(c).init((i: number) => i + cols.indexOf(c));

        const root = view.root.querySelector('[name="viewer-Tile-Viewer"]') as Element;
        const nameOf = (i: Element) => (i.getAttribute('name') || '').replace('input-', '');
        const covered = () => new Set(Array.from(root
          .querySelectorAll('.d4-tile-viewer-form input[name^="input-"]')).map(nameOf));
        for (let i = 0; i < 30; i++) {
          const cov = covered();
          const uncovered = t.columns.names().filter((c: string) => !cov.has(c));
          if (cov.size >= 10 && uncovered.length > 0) break;
          await new Promise((res) => setTimeout(res, 250));
        }
      }, EXTRA_COLS);
      const grown = await fixtureFields();

      await openEditForm('[data-ag-fixture="1"]');
      const hosts = await designerHosts();
      const designerTarget = hosts.values[0];
      await page.locator(`.grok-view-sketch .d4-host[name="div-${designerTarget}"]`)
        .filter({has: page.locator(`input[name="input-${designerTarget}"]`)}).first().click();
      await page.keyboard.press('Delete');
      await page.waitForTimeout(300); 
      await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
      await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
      await page.waitForTimeout(700); 

      const before = await fixtureFields();

      const victim = before.order[0];
      await removeColumnByFieldMenu('[data-ag-fixture="1"]', victim);
      const r = await page.evaluate(async (v: string) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === 'ag-fixture') as any;
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');

        for (let i = 0; i < 24; i++) {
          if (!t.columns.names().includes(v)) break;
          await new Promise((res) => setTimeout(res, 250));
        }
        return {
          columnGoneFromDf: !t.columns.names().includes(v),
          autoGenerateFalse: tileV.props.autoGenerate === false,
          formDesignedTrue: tileV.props.sketchState?.['formDesigned'] === true,
        };
      }, victim);
      const after = await fixtureFields();

      expect(grown.uncovered.length).toBeGreaterThan(0);
      expect(before.tiles).toBeGreaterThan(0);
      expect(before.order.length).toBeGreaterThan(0);
      expect(before.uncovered.length).toBeGreaterThan(0);
      expect(before.order).toContain(victim);

      expect(before.order).not.toContain(designerTarget);
      expect(after.order).not.toContain(designerTarget);
      expect(r.autoGenerateFalse).toBe(true);
      expect(r.formDesignedTrue).toBe(true);
      expect(r.columnGoneFromDf).toBe(true);

      expect(after.covered).not.toContain(victim);

      for (const c of before.uncovered) expect(after.covered).not.toContain(c);
      expect(after.order).toEqual(before.order.filter((n) => n !== victim));
    });
  } finally {
    await page.evaluate(() => {
      const view = Array.from(grok.shell.views).find((x: any) => x.name === 'ag-fixture');
      if (view) view.close();
      const demog = Array.from(grok.shell.views).find((x: any) => x.name === 'Table');
      if (demog) grok.shell.v = demog;
    });
    await page.waitForTimeout(300); 
  }

  v.finishSpec();
});
