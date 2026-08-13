/* ---
realizes: []
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;


test.use(specTestOptions);

test('Tile Viewer tests', async ({page}) => {
  test.setTimeout(300_000);

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

  // Viewer menu needs the ContextMenu key on a focused '.d4-sketch': a real right-click
  // opens the column FIELD menu instead, so [name="viewer"] never resolves.
  const openViewerMenu = async (rootSelector = '[name="viewer-Tile-Viewer"]'): Promise<void> => {
    await page.locator(`${rootSelector} .d4-tile-viewer-form .d4-sketch`).first().focus();
    await page.keyboard.press('ContextMenu');
    await page.locator('.d4-menu-popup[name="viewer"]').waitFor({timeout: 15000});
  };

  // Real gesture: a programmatic removal skips this menu. Right-click the field's VALUE
  // input → column-named FIELD menu → 'Remove' (no 'Delete' on this build), no dialog.
  const columnSlug = (col: string) => col.replace(/[^A-Za-z0-9]/g, '-');
  const removeColumnViaFieldMenu = async (col: string, rootSelector: string): Promise<void> => {
    const slug = columnSlug(col);
    await page.locator(`${rootSelector} .d4-tile-viewer-form input[name="input-${slug}"]`).first()
      .click({button: 'right'});
    const popup = page.locator(`.d4-menu-popup[name="${slug}"]`);
    await popup.first().waitFor({timeout: 15000});
    await popup.first().locator('.d4-menu-item[name="div-Remove"]').first().click();
    await page.waitForTimeout(600);
  };

  const clickTile = async (displayIdx: number, modifiers: string[] = []): Promise<void> => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const box = await tiles[displayIdx].boundingBox();
    for (const m of modifiers) await page.keyboard.down(m);
    await page.mouse.click(box!.x + 15, box!.y + 15);
    for (const m of [...modifiers].reverse()) await page.keyboard.up(m);
    await page.waitForTimeout(300);
  };

  // DataFrame rows matching everything a visible tile displays, resolved from the tile's
  // own field values so a click assert runs DOM → DataFrame. Returns every match.
  const rowsMatchingTile = async (displayIdx: number): Promise<number[]> => page.evaluate((idx: number) => {
    const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
    const df = grok.shell.tv.dataFrame;
    const tile = root.querySelectorAll('.d4-tile-viewer-form')[idx];
    const names = df.columns.names();
    // name= slugifies non-alphanumerics to '-' (DIS_POP → input-DIS-POP).
    const shown = Array.from(tile.querySelectorAll('input[name^="input-"]'))
      .map((i) => [names.find((n: string) => `input-${n.replace(/[^A-Za-z0-9]/g, '-')}` === i.getAttribute('name')),
        (i as HTMLInputElement).value])
      .filter((p) => p[0] != null);
    const matches: number[] = [];
    for (let r = 0; r < df.rowCount; r++) {
      if (shown.every((p) => df.col(p[0]).getString(r) === p[1]))
        matches.push(r);
    }
    return matches;
  }, displayIdx);

  // Idempotent: the Context Panel is shared state, so acting on the table leaves a stale
  // unclickable row that presence checks accept — re-click the gear each pass to drive it
  // VISIBLE. Data rows ship expanded (own tbody) so pass no category; collapsed ones do.
  const propRow = async (prop: string, category?: string) => {
    const row = page.locator(`.property-grid tr[name="prop-${prop}"]`).first();
    for (let attempt = 0; attempt < 4; attempt++) {
      if (await row.isVisible().catch(() => false))
        return row;
      await v.clickViewerTitlebarIcon(page, 'Tile-Viewer', 'icon-font-icon-settings').catch(() => {});
      await page.waitForTimeout(700);
      if (category != null)
        await v.ensurePropertyCategory(page, 'Tile-Viewer', category, prop).catch(() => {});
      await row.waitFor({state: 'visible', timeout: 5000}).catch(() => {});
    }
    // last pass: fail with the honest visibility error, not a click timeout.
    await row.waitFor({state: 'visible', timeout: 5000});
    return row;
  };

  // The filter GROUP can be empty later, leaving the host on screen with no `.d4-filter`;
  // openFilterPanel waits for one, so use it only to bring the panel up when absent.
  const ensureFilterPanel = async (): Promise<void> => {
    const host = page.locator('[name="viewer-Filters"]');
    if (await host.count() === 0) {
      await v.openFilterPanel(page);
      return;
    }
    await host.first().waitFor({state: 'visible', timeout: 10000});
  };

  // ---- Default form rendering ----
  // Direction DOM → DataFrame: resolve the row from the clicked tile's field values
  // BEFORE the click, then cross '.d4-current' on that tile with df.currentRowIdx.
  let firstTileRow = -1;
  let secondTileRow = -1;
  await softStep('Default form rendering: clicking the first tile makes its row current', async () => {
    const matches = await rowsMatchingTile(0);
    // the tile's field values must identify exactly one row, else the cross-check is ambiguous.
    expect(matches.length).toBe(1);
    firstTileRow = matches[0];
    await clickTile(0);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        tileCurrent: tiles[0].classList.contains('d4-current'),
        currentTiles: tiles.filter((t) => t.classList.contains('d4-current')).length,
        idx: grok.shell.tv.dataFrame.currentRowIdx,
      };
    });
    expect(r.tileCurrent).toBe(true);
    expect(r.currentTiles).toBe(1);
    expect(r.idx).toBe(firstTileRow);
  });

  await softStep('Default form rendering: clicking the second tile moves current to that row', async () => {
    const matches = await rowsMatchingTile(1);
    expect(matches.length).toBe(1);
    secondTileRow = matches[0];
    // the two tiles show different rows, so the move is observable.
    expect(secondTileRow).not.toBe(firstTileRow);
    await clickTile(1);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        tileCurrent: tiles[1].classList.contains('d4-current'),
        firstStillCurrent: tiles[0].classList.contains('d4-current'),
        currentTiles: tiles.filter((t) => t.classList.contains('d4-current')).length,
        idx: grok.shell.tv.dataFrame.currentRowIdx,
      };
    });
    expect(r.tileCurrent).toBe(true);
    expect(r.firstStillCurrent).toBe(false);
    expect(r.currentTiles).toBe(1);
    expect(r.idx).toBe(secondTileRow);
  });

  // ---- Row selection ----
  // processRowClick: plain click sets currentRow only; Shift SETS one row's selection
  // (additive, not a range); Ctrl TOGGLES. Both DOM class and DataFrame channel crossed.
  await softStep('Row selection: plain click sets current, Shift adds one row to the selection', async () => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));

    // read the highlight BEFORE the Shift-click, which moves current onto its own tile.
    await clickTile(0);
    const plain = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const t = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        highlighted: t[0].classList.contains('d4-current'),
        currentTiles: t.filter((x) => x.classList.contains('d4-current')).length,
        idx: grok.shell.tv.dataFrame.currentRowIdx,
      };
    });
    expect(plain.highlighted).toBe(true);
    expect(plain.currentTiles).toBe(1);
    expect(plain.idx).toBe(firstTileRow);

    const box2 = await tiles[2].boundingBox();
    await page.keyboard.down('Shift');
    await page.mouse.click(box2!.x + 10, box2!.y + 10);
    await page.keyboard.up('Shift');
    await page.waitForTimeout(300);

    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const selected = root.querySelectorAll('.d4-tile-viewer-form.d4-selected').length;
      // rows between the two clicked tiles must NOT be selected (Shift is additive, not a range).
      const tiles = root.querySelectorAll('.d4-tile-viewer-form');
      const midHasSelected = tiles[1]?.classList.contains('d4-selected');
      return {
        selCount: df.selection.trueCount, selectedTiles: selected, midHasSelected,
        // identity, not just the count: the Shift-clicked tile is the selected one.
        thirdSelected: tiles[2]?.classList.contains('d4-selected'),
      };
    });
    expect(r.selCount).toBe(1);
    expect(r.selectedTiles).toBe(1);
    expect(r.thirdSelected).toBe(true);
    expect(r.midHasSelected).toBe(false);
  });

  await softStep('Row selection: Ctrl-click adds the fifth tile to the selection', async () => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    const box = await tiles[4].boundingBox();
    await page.keyboard.down('Control');
    await page.mouse.click(box!.x + 10, box!.y + 10);
    await page.keyboard.up('Control');
    await page.waitForTimeout(300);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        selCount: df.selection.trueCount,
        selectedTiles: tiles.filter((t) => t.classList.contains('d4-selected')).length,
        // identity: the Ctrl-clicked tile joined the Shift-clicked one.
        fifthSelected: tiles[4]?.classList.contains('d4-selected'),
        thirdStillSelected: tiles[2]?.classList.contains('d4-selected'),
      };
    });
    expect(r.selCount).toBe(before + 1);
    expect(r.selectedTiles).toBe(r.selCount);
    expect(r.fifthSelected).toBe(true);
    expect(r.thirdStillSelected).toBe(true);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
  });

  // ---- Lanes ----
  // ACTUATION HAZARD (atlas): the lanes VALUE is a canvas column grid with no per-column
  // selector, so set the column via JS API. Assert lane STRUCTURE in the DOM (count,
  // ordered header texts, per-lane class) — a different channel, never a prop echo.
  await softStep('Lanes: set RACE → one lane per category, headers match RACE categories', async () => {
    await propRow('lanes');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props.lanesColumnName = 'RACE';
      await new Promise((res) => setTimeout(res, 600));
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      return {
        laneCount: lanes.length,
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
        allMulti: lanes.every((l) => l.classList.contains('d4-tile-viewer-lane-multi')),
        categories: grok.shell.tv.dataFrame.col('RACE').categories.slice(),
      };
    });
    expect(r.laneCount).toBe(r.categories.length);
    expect(r.headers).toEqual(r.categories);
    expect(r.allMulti).toBe(true);
  });

  await softStep('Lanes: set SEX → grouping updates to two lanes (F, M)', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props.lanesColumnName = 'SEX';
      await new Promise((res) => setTimeout(res, 600));
      return {
        headers: Array.from(root.querySelectorAll('.d4-tile-viewer-lane-header')).map((h) => h.textContent),
        laneCount: root.querySelectorAll('.d4-tile-viewer-lane').length,
        categories: grok.shell.tv.dataFrame.col('SEX').categories.slice(),
      };
    });
    expect(r.laneCount).toBe(r.categories.length);
    expect(r.headers).toEqual(expect.arrayContaining(['F', 'M']));
  });

  await softStep('Lanes: clear → a single flat lane', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props.lanesColumnName = null;
      await new Promise((res) => setTimeout(res, 600));
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      return {
        laneCount: lanes.length,
        single: lanes[0]?.classList.contains('d4-tile-viewer-lane-single'),
        // the single lane must actually hold tiles — an EMPTY lane would satisfy the asserts alone.
        tilesInLane: lanes[0]?.querySelectorAll('.d4-tile-viewer-form').length ?? 0,
      };
    });
    expect(r.laneCount).toBe(1);
    expect(r.single).toBe(true);
    // Population floor only: the lane is a virtualized window, so tile count depends on
    // the viewport and cannot be asserted exactly.
    expect(r.tilesInLane).toBeGreaterThan(0);
  });

  // ---- Row source ----
  // `prop-row-source` (Data) is an ordinary choice editor, so commit the value through it —
  // an uncommitted editor fails the step instead of passing.
  await softStep('Row source: Selected → tiles show only the selected rows', async () => {
    // Select rows FIRST: a selection replaces the Context Panel, so bring the property row up after.
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (const i of [0, 1, 2, 3, 4]) df.selection.set(i, true);
    });
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Selected');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const selectedIdx = [0, 1, 2, 3, 4];
      await new Promise((res) => setTimeout(res, 700));
    // Read what the PRODUCT rendered: the tile population must be exactly the selected
    // rows, not a read-back of the selection the test set.
      const tileAges = Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name="input-AGE"]'))
        .map((i) => (i as HTMLInputElement).value).sort();
      const selectedAges = selectedIdx.map((i) => df.col('AGE').getString(i)).sort();
      return {
        visibleTiles: root.querySelectorAll('.d4-tile-viewer-form').length,
        selectedRows: selectedIdx.length,
        tileAges,
        selectedAges,
      };
    });
    // product-rendered population: one tile per selected row, carrying the selected rows' AGE.
    expect(r.visibleTiles).toBe(r.selectedRows);
    expect(r.tileAges).toEqual(r.selectedAges);
  });

  await softStep('Row source: Filtered + SEX=M → only matching rows on the tiles', async () => {
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Filtered');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise((res) => setTimeout(res, 800));
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        filterCount: df.filter.trueCount,
        total: df.rowCount,
        visibleTiles: tiles.length,
        allTilesM: tiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M'),
      };
    });
    expect(r.filterCount).toBeLessThan(r.total);
    expect(r.visibleTiles).toBeGreaterThan(0);
    // every rendered tile is male — the tile population, not just the count.
    expect(r.allTilesM).toBe(true);
  });

  await softStep('Row source: All → every row is shown again, filter cleared', async () => {
    // Entry state (prior step): rowSource='Filtered' + SEX=M, every tile male. Capture that
    // constrained population before switching, so the change is observable.
    const entry = await page.evaluate(() => {
      const tiles = Array.from(document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form'));
      return {beforeAllMale: tiles.length > 0 &&
        tiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M')};
    });
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'All');
    // ISOLATING WINDOW: commit 'All' while SEX=M is STILL applied — the only reading that
    // can falsify the setting itself; once the filter is removed females return regardless.
    const isolated = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      await new Promise((res) => setTimeout(res, 900));
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        tiles: tiles.length,
        hasFemale: tiles.some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F'),
        filterCount: df.filter.trueCount,
        total: df.rowCount,
      };
    });
    // viewer-local: a row the DataFrame filter excludes still appears on the tiles ...
    expect(isolated.tiles).toBeGreaterThan(0);
    expect(isolated.hasFemale).toBe(true);
    // ... while df.filter is untouched, so 'All' bypasses it rather than the filter having been cleared.
    expect(isolated.filterCount).toBeLessThan(isolated.total);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
    // Remove via the panel only — setAll(true) would make trueCount==rowCount a tautology
    // on what the test set rather than what the product restored.
      for (const f of [...fg.filters]) fg.remove(f);
      df.selection.setAll(false);
      await new Promise((res) => setTimeout(res, 900));
      // Under 'All' at least one rendered tile is now female — a product-rendered change.
      const afterTiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const afterHasFemale = afterTiles.some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F');
      return {afterHasFemale, filterCount: df.filter.trueCount, total: df.rowCount};
    });
    expect(entry.beforeAllMale).toBe(true);
    expect(r.afterHasFemale).toBe(true);
    // filter restored by the product after the panel filter was removed (no setAll).
    expect(r.filterCount).toBe(r.total);
  });

  // ---- Tiles font ----
  // Signal: the lane header's inline style.font + lineHeight = size * 1.4 (product-computed;
  // see render-signal-index). getComputedStyle reads '' on a just-restyled header headless,
  // so the inline value is the durable channel; the tile field font is a second surface.
  // ACTUATION: `prop-tiles-font` (Style) is a FontInput — size box, family select
  // (Monospace / Arial / Open Sans / Roboto only), B/I toggles.
  const fontRow = () => page.locator('.property-grid tr[name="prop-tiles-font"]');
  await softStep('Tiles font: size 18px grows the lane headers and tile text', async () => {
    // group by RACE first: lane headers are hidden in single-lane mode (lanes via JS API, hazard above).
    await page.evaluate(async () => {
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props.lanesColumnName = 'RACE';
      await new Promise((res) => setTimeout(res, 600));
    });
    await propRow('tiles-font', 'style');
    const readHeader = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const h = root.querySelector('.d4-tile-viewer-lane-header') as HTMLElement;
      const tileInput = root.querySelector('.d4-tile-viewer-form input[name="input-AGE"]') as HTMLElement;
      return {
        size: h.style.fontSize, line: h.style.lineHeight,
        // the tile field reads reliably via getComputedStyle (not restyled in place), so
        // it stays the tile-body cross-channel.
        tile: getComputedStyle(tileInput).fontSize,
      };
    });
    await fontRow().locator('input.d4-font-size-input').fill('13');
    await page.waitForTimeout(800);
    const base = await readHeader();
    await fontRow().locator('input.d4-font-size-input').fill('18');
    await page.waitForTimeout(900);
    const grown = await readHeader();
    // base: 13px, lineHeight 13 * 1.4 = 18.2px (product-computed).
    expect(base.size).toBe('13px');
    expect(base.line).toBe('18.2px');
    // grown: 18px, lineHeight 18 * 1.4 = 25.2px — the *1.4 is an independent channel, not an echo.
    expect(grown.size).toBe('18px');
    expect(grown.line).toBe('25.2px');
    expect(grown.tile).toBe('18px');
  });

  await softStep('Tiles font: the family choice reaches the lane header and the tiles', async () => {
    await propRow('tiles-font', 'style');
    await fontRow().locator('select').selectOption('Arial');
    await page.waitForTimeout(900);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const h = root.querySelector('.d4-tile-viewer-lane-header') as HTMLElement;
      const tileInput = root.querySelector('.d4-tile-viewer-form input[name="input-AGE"]') as HTMLElement;
      return {inlineFont: h.style.font, tileFamily: getComputedStyle(tileInput).fontFamily};
    });
    // the family reaches the inline font string ...
    expect(r.inlineFont).toContain('Arial');
    // ... and the tile body, the second surface (SketchForm hostFont).
    expect(r.tileFamily).toContain('Arial');
  });

  await softStep('Tiles font: reset to default 13px Roboto restores the header font', async () => {
    await propRow('tiles-font', 'style');
    await fontRow().locator('select').selectOption('Roboto');
    await page.waitForTimeout(600);
    await fontRow().locator('input.d4-font-size-input').fill('13');
    await page.waitForTimeout(600);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
    // Read inline style.font (written synchronously; avoids the getComputedStyle '' race
    // headless). Poll to let the debounced refresh settle.
      let size = '';
      let line = '';
      let font = '';
      for (let i = 0; i < 20; i++) {
        await new Promise((res) => setTimeout(res, 150));
        const h = root.querySelector('.d4-tile-viewer-lane-header') as HTMLElement | null;
        if (!h) continue;
        size = h.style.fontSize;
        line = h.style.lineHeight;
        font = h.style.font;
        if (size === '13px' && /Roboto/.test(font)) break;
      }
      return {size, line, font};
    });
    expect(r.size).toBe('13px');
    expect(r.line).toBe('18.2px');
    expect(r.font).toContain('Roboto');
  });

  // ---- Auto-generate on column change ----
  // Gesture: a COLUMN delete via onColumnsChanged (not a rebind). Both branches run on a
  // demog CLONE in its own view with every DOM read scoped to that viewer's root — a
  // page-wide selector resolves to the main demog viewer — and close the view in finally.
  // See render-signal-index: field cap, rebind-nulls-sketchState, refill discriminator.
  const AG_FIXTURE = 'ag-fixture';
  try {
    await softStep('Auto-generate (auto state): deleting a fielded column removes it and frees a slot for an excluded column', async () => {
      const entry = await page.evaluate(async (frameName: string) => {
        const t = grok.shell.tv.dataFrame.clone();
        t.name = frameName;
        const view = grok.shell.addTableView(t);
        const tileV = view.addViewer('Tile Viewer');
        // tag the viewer's OWN root so Playwright locators (the removal gesture included)
        // address this viewer, not the main demog one.
        const root = tileV.root as Element;
        root.setAttribute('data-ag-fixture', '1');
        // match columns through the same name= slug (see above).
        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;
        // which columns the product actually put on a tile — read from the rendered value
        // inputs, never from a column list the test assumed.
        const rendered = () => {
          const hosts = new Set(Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
            .map((i) => i.getAttribute('name')));
          return t.columns.names().filter((c: string) => hosts.has(slug(c)));
        };
        let renderedBefore: string[] = [];
        for (let i = 0; i < 20; i++) {
          await new Promise((res) => setTimeout(res, 500));
          renderedBefore = rendered();
          if (renderedBefore.length > 0) break;
        }
        return {
          // entry-state guard: a fresh viewer is auto-generated on BOTH channels.
          autoGenerateTrue: tileV.props.autoGenerate === true,
          formNotDesigned: tileV.props.sketchState?.['formDesigned'] === false,
          totalCols: t.columns.length,
          renderedBefore,
          excludedBefore: t.columns.names().filter((c: string) => !renderedBefore.includes(c)),
          // DELETE a column that HAS a field on the card, and one that's CLICKABLE: a
          // boolean column renders as a disabled checkbox whose field-menu click never
          // enables, so it would time out. Composition still counts every field.
          victim: renderedBefore.filter((c: string) => {
            const el = root.querySelector(`.d4-tile-viewer-form input[name="${slug(c)}"]`) as HTMLInputElement | null;
            return el != null && el.type !== 'checkbox' && !el.disabled;
          })[0] ?? null,
        };
      }, AG_FIXTURE);
      // The delete is the field menu's Remove item on the victim's own tile field.
      if (entry.victim != null)
        await removeColumnViaFieldMenu(entry.victim, '[data-ag-fixture="1"]');
      const after = await page.evaluate(async (args: {frameName: string, victim: string | null, excluded: string[]}) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === args.frameName);
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = tileV.root as Element;
        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;
        const victimField = args.victim == null ? null : slug(args.victim);
    // DOM channel — read field names RAW off the tiles. Filtering through the frame's
    // column list can never fail: the removed column cannot appear there.
        const victimHosts = () => victimField == null ? 0
          : root.querySelectorAll(`.d4-tile-viewer-form input[name="${victimField}"]`).length;
        // the intersection of the two channels, used only for the refill/cap reads.
        const rendered = () => {
          const hosts = new Set(Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
            .map((i) => i.getAttribute('name')));
          return t.columns.names().filter((c: string) => hosts.has(slug(c)));
        };
        let renderedAfter = rendered();
        let refilled: string[] = [];
        for (let i = 0; i < 24 && args.victim != null; i++) {
          renderedAfter = rendered();
          refilled = renderedAfter.filter((c: string) => args.excluded.includes(c));
          // wait on the DOM channel, not on the frame-filtered set.
          if (victimHosts() === 0 && refilled.length > 0) break;
          await new Promise((res) => setTimeout(res, 500));
        }
        // a surviving field resolves to its DataFrame display string; the first tile in
        // DOM order is the frame's first row.
        const keeper = renderedAfter[0] ?? null;
        const tile = root.querySelector('.d4-tile-viewer-form');
        return {
          renderedAfter, refilled, keeper,
          keeperValue: keeper == null ? null
            : (tile?.querySelector(`input[name="${slug(keeper)}"]`) as HTMLInputElement)?.value,
          keeperDisplay: keeper == null ? null : t.col(keeper).getString(0),
          // DOM channel: no tile carries the victim's input any more.
          victimHostsOnTiles: victimHosts(),
          victimGoneEverywhere: Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
            .every((tl) => victimField == null || !tl.querySelector(`input[name="${victimField}"]`)),
          renderedTiles: root.querySelectorAll('.d4-tile-viewer-form').length,
          // DataFrame channel: the gesture really took the column out of the frame.
          victimStillInFrame: args.victim != null && t.columns.names().includes(args.victim),
          totalColsAfter: t.columns.length,
        };
      }, {frameName: AG_FIXTURE, victim: entry.victim, excluded: entry.excludedBefore});
      const r = {...entry, ...after};
      expect(r.autoGenerateTrue).toBe(true);
      expect(r.formNotDesigned).toBe(true);
      // the card rendered and the cap left the overflow out — a zero here means the DOM
      // read hit the wrong viewer, making the step vacuous.
      expect(r.renderedBefore.length).toBe(10);
      expect(r.renderedBefore.length).toBeLessThan(r.totalCols);
      expect(r.excludedBefore.length).toBeGreaterThan(0);
      // the deleted column was on the card before the delete (no vacuous victim).
      expect(r.victim).not.toBeNull();
      // DOM CHANNEL: the victim's input is gone from EVERY rendered tile (read raw) — the
      // tiles really rendered, so the absence means something.
      expect(r.renderedTiles).toBeGreaterThan(0);
      expect(r.victimHostsOnTiles).toBe(0);
      expect(r.victimGoneEverywhere).toBe(true);
      // DATAFRAME CHANNEL, independent of the DOM: the column left the frame.
      expect(r.victimStillInFrame).toBe(false);
      expect(r.totalColsAfter).toBe(r.totalCols - 1);
      // rendered set != entry set — kept for the evidence it prints on failure (both name
      // lists), not as a discriminator: the victim left the frame, so it cannot come back equal.
      expect(r.renderedAfter).not.toEqual(r.renderedBefore);
      // ... and regeneration pulled the previously field-less column onto the card,
      // refilling the freed slot up to the cap.
      expect(r.refilled.length).toBeGreaterThan(0);
      expect(r.renderedAfter.length).toBe(10);
      expect(r.keeperValue).toBe(r.keeperDisplay);
    });

  // DESIGNED state: Edit Form → delete one field host → CLOSE AND APPLY (autoGenerate
  // false, formDesigned true), then delete a column. Discriminator is the REFILL, not the
  // deleted column's field (gone in either state). Edit Form opens a view, not a dialog.
    await softStep('Auto-generate (designed state): the same column delete does not refill the freed slot', async () => {
    // both columns read off the card itself; alphanumeric names only, so the host is
    // 'div-<column>' with no slug to guess. Skip boolean columns (disabled checkbox, click
    // times out — see above).
      const picked = await page.evaluate((frameName: string) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === frameName);
        grok.shell.v = view;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const onCard = Array.from((tileV.root as Element)
          .querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
          .filter((i) => (i as HTMLInputElement).type !== 'checkbox' && !(i as HTMLInputElement).disabled)
          .map((i) => i.getAttribute('name')!.replace(/^input-/, ''))
          .filter((n) => /^[A-Za-z0-9]+$/.test(n));
        return {designCol: onCard[0] ?? null, frameCol: onCard[1] ?? null};
      }, AG_FIXTURE);
      expect(picked.designCol).not.toBeNull();
      expect(picked.frameCol).not.toBeNull();

      await openViewerMenu('[data-ag-fixture="1"]');
      await page.locator('.d4-menu-popup[name="viewer"] .d4-menu-item[name="div-Edit-Form..."]').click();
      await page.locator('.grok-view-sketch').waitFor({timeout: 15000});
      // a real click selects AND focuses the field host (selection rides on mousedown) so
      // the next Delete removes it. Label and value hosts share the name — the value host
      // is the one holding the value input.
      await page.locator(`.grok-view-sketch .d4-host[name="div-${picked.designCol}"]`)
        .filter({has: page.locator(`input[name="input-${picked.designCol}"]`)}).first().click();
      await page.keyboard.press('Delete');
      await page.waitForTimeout(300);
      await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
      await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
      await page.waitForTimeout(600);

      const pre = await page.evaluate((args: {frameName: string, designCol: string, frameCol: string}) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === args.frameName);
        grok.shell.v = view;
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = tileV.root as Element;
        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;
        // one tile is the layout unit — every tile renders the same form state, so its
        // field names in DOM order carry composition and order.
        const before = (() => {
          const tile = root.querySelector('.d4-tile-viewer-form');
          return tile == null ? []
            : Array.from(tile.querySelectorAll('input[name^="input-"]')).map((i) => i.getAttribute('name'));
        })();
        return {
          autoGenerateFalse: tileV.props.autoGenerate === false,
          formDesignedTrue: tileV.props.sketchState?.['formDesigned'] === true,
          designFieldRemoved: !before.includes(slug(args.designCol)),
          before,
          // columns with no field on the card — the refill candidates, here including the
          // field just deleted in the designer.
          excludedBefore: t.columns.names().filter((c: string) => !before.includes(slug(c))).map(slug),
          victimField: slug(args.frameCol),
        };
      }, {frameName: AG_FIXTURE, designCol: picked.designCol!, frameCol: picked.frameCol!});
      // same gesture as the auto branch — the field menu's Remove on a fielded column.
      await removeColumnViaFieldMenu(picked.frameCol!, '[data-ag-fixture="1"]');
      const post = await page.evaluate(async (args: {frameName: string, victimField: string, excluded: string[]}) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === args.frameName);
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = tileV.root as Element;
        const fieldsOnTile = () => {
          const tile = root.querySelector('.d4-tile-viewer-form');
          return tile == null ? []
            : Array.from(tile.querySelectorAll('input[name^="input-"]')).map((i) => i.getAttribute('name'));
        };
        let after = fieldsOnTile();
        for (let i = 0; i < 20; i++) {
          after = fieldsOnTile();
          if (!after.includes(args.victimField)) break;
          await new Promise((res) => setTimeout(res, 300));
        }
        return {
          after,
          refilled: after.filter((n) => args.excluded.includes(n!)),
          victimGoneEverywhere: Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
            .every((tl) => !tl.querySelector(`input[name="${args.victimField}"]`)),
        };
      }, {frameName: AG_FIXTURE, victimField: pre.victimField, excluded: pre.excludedBefore});
      const r = {...pre, ...post, survivorsBefore: pre.before.filter((n) => n !== pre.victimField)};
      expect(r.autoGenerateFalse).toBe(true);
      expect(r.formDesignedTrue).toBe(true);
      expect(r.designFieldRemoved).toBe(true);
      // a field bound to a column that no longer exists cannot survive ...
      expect(r.after).not.toContain(r.victimField);
      expect(r.victimGoneEverywhere).toBe(true);
      // ... but the freed slot stays empty: a candidate WAS available (the designer-deleted
      // field), yet the product did not regenerate the card.
      expect(r.excludedBefore.length).toBeGreaterThan(0);
      expect(r.refilled).toEqual([]);
      // the surviving fields keep both composition and order.
      expect(r.after).toEqual(r.survivorsBefore);
    });
  } finally {
    await page.evaluate((frameName: string) => {
      const view = Array.from(grok.shell.views).find((x: any) => x.name === frameName);
      if (view) view.close();
      const demog = Array.from(grok.shell.views).find((x: any) => x.name === 'Table');
      if (demog) grok.shell.v = demog;
    }, AG_FIXTURE);
    await page.waitForTimeout(300);
  }

  // ---- Multiple table switching (via property Table) ----
  // `prop-table` (Data) is a select of the open frames. PRIMARY signal: a tile field value
  // against the spgi-100 cell it must display; the regenerated field set is the second channel.
  await softStep('Multiple table switching: Table property → spgi-100 shows spgi-100 rows', async () => {
    await page.evaluate(async () => {
      const spgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      spgi.name = 'spgi-100';
      grok.shell.addTableView(spgi);
      await new Promise((res) => setTimeout(res, 700));
      const demogView = Array.from(grok.shell.views).find((x: any) => x.name === 'Table');
      if (demogView) grok.shell.v = demogView;
      await new Promise((res) => setTimeout(res, 300));
    });
    await propRow('table');
    await v.selectPropertyGridChoice(page, 'table', 'spgi-100');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      await new Promise((res) => setTimeout(res, 900));
      const fields = Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
        .map((i) => i.getAttribute('name'));
      // a tile field against the cell it must show: resolve the first rendered field to
      // its spgi-100 column and compare displayed strings.
      const bound = tileV.dataFrame;
      const names = bound?.columns?.names() ?? [];
      const shown = names.find((n: string) => fields.includes(`input-${n.replace(/[^A-Za-z0-9]/g, '-')}`)) ?? null;
      const tile = root.querySelector('.d4-tile-viewer-form');
      return {
        boundTable: bound?.name,
        demogFieldAbsent: !fields.includes('input-SEX'),
        hasFields: fields.length > 0,
        shown,
        tileValue: shown == null ? null
          : (tile?.querySelector(`input[name="input-${shown.replace(/[^A-Za-z0-9]/g, '-')}"]`) as HTMLInputElement)?.value,
        cellValue: shown == null ? null : bound.col(shown).getString(0),
      };
    });
    expect(r.boundTable).toBe('spgi-100');
    expect(r.hasFields).toBe(true);
    expect(r.demogFieldAbsent).toBe(true);
    // the tile renders an spgi-100 cell, not just an spgi-100-shaped field set.
    expect(r.shown).not.toBeNull();
    expect(r.tileValue).toBe(r.cellValue);
  });

  await softStep('Multiple table switching: Table property → demog shows demog rows again', async () => {
    await propRow('table');
    await v.selectPropertyGridChoice(page, 'table', 'Table');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      await new Promise((res) => setTimeout(res, 900));
      const tile = root.querySelector('.d4-tile-viewer-form');
      return {
        boundTable: tileV.dataFrame?.name,
        hasSex: !!tile?.querySelector('input[name="input-SEX"]'),
        tileSex: (tile?.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value,
        cellSex: tileV.dataFrame?.col('SEX')?.getString(0),
      };
    });
    expect(r.boundTable).toBe('Table');
    expect(r.hasSex).toBe(true);
    expect(r.tileSex).toBe(r.cellSex);
  });

  // ---- Context menu items ----
  // Presence-only assert on the viewer's three contributions in the [name="viewer"] popup
  // (ContextMenu key, see openViewerMenu).
  await softStep('Context menu: Edit Form..., Lanes, and Show Empty Lanes are present', async () => {
    await page.evaluate(() => document.body.click());
    await openViewerMenu();
    const r = await page.evaluate(() => {
      const popup = document.querySelector('.d4-menu-popup[name="viewer"]');
      const out = {
        popupOpen: !!popup,
        editForm: !!document.querySelector('.d4-menu-item[name="div-Edit-Form..."]'),
        lanes: !!document.querySelector('.d4-menu-item[name="div-Lanes-"]'),
        showEmptyLanes: !!document.querySelector('.d4-menu-item[name="div-Show-Empty-Lanes"]'),
      };
      document.body.click();
      return out;
    });
    expect(r.popupOpen).toBe(true);
    expect(r.editForm).toBe(true);
    expect(r.lanes).toBe(true);
    expect(r.showEmptyLanes).toBe(true);
  });

  // ---- Viewer title and description ----
  // Signals: the panel titlebar text and a '.d4-viewer-description' node — rendered DOM,
  // not a property read-back. All three rows are ordinary Description-category editors.
  await softStep('Viewer title/description: title appears in the header, description below it', async () => {
    await propRow('title', 'description');
    await v.setPropertyGridValue(page, 'title', 'Patient Cards');
    await v.setPropertyGridValue(page, 'description', 'Demographic data per patient');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const panel = root.closest('.panel-base')!;
      await new Promise((res) => setTimeout(res, 500));
      return {
        titlebar: panel.querySelector('.panel-titlebar-text')?.textContent,
        descText: root.querySelector('.d4-viewer-description')?.textContent,
      };
    });
    expect(r.titlebar).toBe('Patient Cards');
    expect(r.descText).toBe('Demographic data per patient');
  });

  await softStep('Viewer title/description: Description Position Bottom moves it below the tiles', async () => {
    // Signal: the description node's ORDER relative to the lanes host — what the position
    // flips (refdoc). Coordinates would be a fragile stand-in.
    const order = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const desc = root.querySelector('.d4-viewer-description')!;
      const lanes = root.querySelector('.d4-tile-viewer-lanes-host')!;
      // DOCUMENT_POSITION_FOLLOWING (4): lanes come after the description.
      return {descBeforeLanes: (desc.compareDocumentPosition(lanes) & 4) !== 0};
    });
    await propRow('description-position', 'description');
    await v.selectPropertyGridChoice(page, 'description-position', 'Top');
    const top = await order();
    await v.selectPropertyGridChoice(page, 'description-position', 'Bottom');
    const bottom = await order();
    expect(top.descBeforeLanes).toBe(true);
    expect(bottom.descBeforeLanes).toBe(false);
  });

    // The titlebar distinguishes empty title from no title: unset shows the viewer type, an
    // empty string renders empty, only null restores the type name. The property grid writes
    // an empty string when cleared. See tile.md.
  await softStep('Viewer title/description: clearing both leaves the header without a title', async () => {
    await propRow('title', 'description');
    // clearing needs the selection deleted: the setter TYPES the value, and typing an empty
    // string would leave the old text in place.
    const clearRow = async (prop: string): Promise<void> => {
      await page.locator(`.property-grid tr[name="prop-${prop}"]`).first().locator('td').last().click();
      await page.waitForTimeout(400);
      await page.keyboard.press('Control+a');
      await page.keyboard.press('Delete');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(700);
    };
    await clearRow('title');
    await clearRow('description');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const panel = root.closest('.panel-base')!;
      await new Promise((res) => setTimeout(res, 400));
      return {
        // the titlebar node must still BE there — an empty read off a vanished node would
        // satisfy the assert for the wrong reason.
        titlebarPresent: !!panel.querySelector('.panel-titlebar-text'),
        titlebar: (panel.querySelector('.panel-titlebar-text')?.textContent ?? '').trim(),
        descPresent: !!root.querySelector('.d4-viewer-description'),
      };
    });
    expect(r.titlebarPresent).toBe(true);
    expect(r.titlebar).toBe('');
    expect(r.descPresent).toBe(false);
  });

  // ---- Filter interaction ----
  // Counts come from the DataFrame at runtime, never hardcoded. Panel via the registered
  // helpers; individual filters via the filter group (canvas widgets, no per-category selector).
  await softStep('Filter interaction: SEX=M reduces the tiles to male patients', async () => {
    await ensureFilterPanel();
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Filtered');
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of [...fg.filters]) fg.remove(f);
      df.filter.setAll(true);
      await new Promise((res) => setTimeout(res, 400));
      const total = df.filter.trueCount;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise((res) => setTimeout(res, 700));
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        total, afterSex: df.filter.trueCount,
        // the added filter is a real widget in the panel, not just a bitset change.
        filterWidgets: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
        allTilesM: tiles.length > 0 && tiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M'),
      };
    });
    expect(r.afterSex).toBeLessThan(r.total);
    expect(r.filterWidgets).toBeGreaterThan(0);
    expect(r.allTilesM).toBe(true);
  });

  await softStep('Filter interaction: adding AGE > 50 reduces the tiles further', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const ageOnTiles = () => Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .map((t) => parseFloat((t.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value ?? ''))
        .filter((n) => !isNaN(n));
      const before = df.filter.trueCount;
      const tilesBefore = root.querySelectorAll('.d4-tile-viewer-form').length;
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 51, max: 200});
      await new Promise((res) => setTimeout(res, 900));
      const ages = ageOnTiles();
      return {
        before, after: df.filter.trueCount,
        tilesBefore, tilesAfter: root.querySelectorAll('.d4-tile-viewer-form').length,
        // the viewer-local population itself, not just the counter: every tile satisfies the filter.
        ageReadCount: ages.length,
        allTilesOver50: ages.length > 0 && ages.every((a) => a > 50),
      };
    });
    expect(r.after).toBeLessThan(r.before);
    expect(r.after).toBeGreaterThan(0);
    expect(r.tilesAfter).toBeLessThanOrEqual(r.tilesBefore);
    expect(r.ageReadCount).toBeGreaterThan(0);
    expect(r.allTilesOver50).toBe(true);
  });

  await softStep('Filter interaction: removing all filters restores every tile', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      // Entry state (two prior steps): SEX=M + AGE>50, every tile male. Capture that
      // constrained population before removing filters.
      const beforeTiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const beforeAllMale = beforeTiles.length > 0 &&
        beforeTiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M');
      // Remove filters through the panel only — do NOT setAll(true). The product recomputes
      // df.filter to the full frame and re-renders the female rows itself; asserting after
      // setAll(true) would only prove setAll works (E-VACUOUS).
      for (const f of [...fg.filters]) fg.remove(f);
      await new Promise((res) => setTimeout(res, 900));
      const afterTiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const afterHasFemale = afterTiles.some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F');
      return {beforeAllMale, afterHasFemale, filterCount: df.filter.trueCount, total: df.rowCount};
    });
    // constrained (male) → full-frame population, and the product restores df.filter itself.
    expect(r.beforeAllMale).toBe(true);
    expect(r.afterHasFemale).toBe(true);
    expect(r.filterCount).toBe(r.total);
  });

  await softStep('Filter interaction: closing the filter panel removes it from the view', async () => {
    // the panel must be open for this step to mean anything — assert entry state, don't skip.
    expect(await page.locator('[name="viewer-Filters"]').count()).toBeGreaterThan(0);
    // Every dock panel titlebar carries a [name="Close"] button (panel_dock_container.dart).
    await v.clickViewerTitlebarIcon(page, 'Filters', 'Close');
    await page.locator('[name="viewer-Filters"]').waitFor({state: 'detached', timeout: 10000});
    const r = await page.evaluate(() => ({
      panelGone: document.querySelector('[name="viewer-Filters"]') == null,
      tilesStillRendered: document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').length,
    }));
    expect(r.panelGone).toBe(true);
    expect(r.tilesStillRendered).toBeGreaterThan(0);
  });

  // ---- Viewer filter formula (viewer-local vs DataFrame filter) ----
  // look.filter narrows the viewer's own population and must leave df.filter alone.
  // `prop-filter` (Data) commits on change. Tile COUNT is meaningless under virtualization,
  // so the signal is what the tiles CARRY.
  await softStep('Viewer filter: a formula narrows the tiles while the DataFrame filter stays put', async () => {
    const readAges = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const ages = Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .map((t) => parseFloat((t.querySelector('input[name="input-AGE"]') as HTMLInputElement)?.value ?? ''))
        .filter((n) => !isNaN(n));
      return {ages, filterCount: df.filter.trueCount, total: df.rowCount};
    });
    const before = await readAges();
    const row = await propRow('filter');
    const box = row.locator('input.property-grid-ellipsis-editor-input').first();
    await box.fill('${AGE} > 50');
    await box.press('Enter');
    await page.waitForTimeout(1200);
    const after = await readAges();
    // entry state carried rows the formula rejects, so its effect is observable.
    expect(before.ages.length).toBeGreaterThan(0);
    expect(before.ages.some((a) => a <= 50)).toBe(true);
    // viewer-local: every rendered tile now satisfies the formula ...
    expect(after.ages.length).toBeGreaterThan(0);
    expect(after.ages.every((a) => a > 50)).toBe(true);
    // ... and the DataFrame filter — the grid's own row source — never moved.
    expect(after.filterCount).toBe(before.filterCount);
    expect(after.filterCount).toBe(after.total);

    // clearing the formula gives the rejected rows back.
    const row2 = await propRow('filter');
    const box2 = row2.locator('input.property-grid-ellipsis-editor-input').first();
    await box2.fill('');
    await box2.press('Enter');
    await page.waitForTimeout(1200);
    const cleared = await readAges();
    expect(cleared.ages.some((a) => a <= 50)).toBe(true);
    expect(cleared.filterCount).toBe(cleared.total);
  });

  // ---- Scroll position survives adding a viewer ----
  // The scroll must come from a real wheel gesture: assigning scrollTop leaves the virtual
  // list's own scroll model untouched. Asserted: scrollTop + the identity of the first
  // rendered row — not clientHeight (dock re-layout) nor tile count (virtualized).
  await softStep('Scroll position: a scrolled lane keeps its position and its rows when another viewer is added', async () => {
    // one flat lane = a single long scrollable list (lanes column via JS API, hazard above).
    await page.evaluate(async () => {
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props.lanesColumnName = null;
      await new Promise((res) => setTimeout(res, 900));
    });
    const lane = page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-content').first();
    const box = await lane.boundingBox();
    await page.mouse.move(box!.x + box!.width / 2, box!.y + box!.height / 2);
    for (let i = 0; i < 10; i++) {
      await page.mouse.wheel(0, 300);
      await page.waitForTimeout(120);
    }
    await page.waitForTimeout(1000);
    const readLane = () => page.evaluate(() => {
      const l = document.querySelector('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-content') as HTMLElement;
      const first = l?.querySelector('.d4-tile-viewer-form');
      const val = (n: string) => (first?.querySelector(`input[name="input-${n}"]`) as HTMLInputElement)?.value ?? null;
      return {scrollTop: l?.scrollTop ?? -1, age: val('AGE'), sex: val('SEX'), weight: val('WEIGHT')};
    });
    const before = await readLane();
    await page.evaluate(async () => {
      grok.shell.tv.addViewer('Histogram');
      await new Promise((res) => setTimeout(res, 2500));
    });
    const after = await readLane();
    // the wheel really moved the lane, and a row was rendered to identify ...
    expect(before.scrollTop).toBeGreaterThan(0);
    expect(before.age).not.toBeNull();
    // ... adding a viewer neither reset the scroll nor re-based the rendered window. A couple
    // pixels of tolerance cover the dock re-layout; the row identity must hold exactly.
    expect(after.scrollTop).toBeGreaterThan(0);
    expect(Math.abs(after.scrollTop - before.scrollTop)).toBeLessThanOrEqual(2);
    expect(after.age).toBe(before.age);
    expect(after.sex).toBe(before.sex);
    expect(after.weight).toBe(before.weight);
  });

  v.finishSpec();
});
