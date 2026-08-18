/* ---
realizes: [tileviewer.cp.tiles-font-applied-to-rendering, tileviewer.cp.auto-generate-on-columns-change, tileviewer.cp.table-rebind-regenerates-form, tileviewer.cp.context-menu-inventory, tileviewer.cp.viewer-local-filter-vs-dataframe-filter, tileviewer.cp.scroll-survives-added-viewer, tileviewer.int.viewer-local-filter-vs-df-filter]
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

  const openViewerMenu = async (rootSelector = '[name="viewer-Tile-Viewer"]'): Promise<void> => {
    await page.locator(`${rootSelector} .d4-tile-viewer-form .d4-sketch`).first().focus();
    await page.keyboard.press('ContextMenu');
    await page.locator('.d4-menu-popup[name="viewer"]').waitFor({timeout: 15000});
  };

  const columnSlug = (col: string) => col.replace(/[^A-Za-z0-9]/g, '-');
  const removeColumnViaFieldMenu = async (col: string, rootSelector: string): Promise<void> => {
    const slug = columnSlug(col);
    await page.locator(`${rootSelector} .d4-tile-viewer-form input[name="input-${slug}"]`).first()
      .click({button: 'right'});
    const popup = page.locator(`.d4-menu-popup[name="${slug}"]`);
    await popup.first().waitFor({timeout: 15000});
    await popup.first().locator('.d4-menu-item[name="div-Remove"]').first().click();

    await page.waitForTimeout(300);
  };

  const clickTile = async (displayIdx: number, modifiers: string[] = []): Promise<void> => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const box = await tiles[displayIdx].boundingBox();
    for (const m of modifiers) await page.keyboard.down(m);
    await page.mouse.click(box!.x + 15, box!.y + 15);
    for (const m of [...modifiers].reverse()) await page.keyboard.up(m);

    await page.waitForTimeout(300);
  };

  const rowsMatchingTile = async (displayIdx: number): Promise<number[]> => page.evaluate((idx: number) => {
    const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
    const df = grok.shell.tv.dataFrame;
    const tile = root.querySelectorAll('.d4-tile-viewer-form')[idx];
    const names = df.columns.names();

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

    await row.waitFor({state: 'visible', timeout: 5000});
    return row;
  };

  const ensureFilterPanel = async (): Promise<void> => {
    const host = page.locator('[name="viewer-Filters"]');
    if (await host.count() === 0) {
      await v.openFilterPanel(page);
      return;
    }
    await host.first().waitFor({state: 'visible', timeout: 10000});
  };

  let firstTileRow = -1;
  let secondTileRow = -1;
  await softStep('Default form rendering: clicking the first tile makes its row current', async () => {
    const matches = await rowsMatchingTile(0);

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

  await softStep('Row selection: plain click sets current, Shift adds one row to the selection', async () => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));

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

      const tiles = root.querySelectorAll('.d4-tile-viewer-form');
      const midHasSelected = tiles[1]?.classList.contains('d4-selected');
      return {
        selCount: df.selection.trueCount, selectedTiles: selected, midHasSelected,

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

  const setTileProp = async (prop: string, value: any, capMs = 900): Promise<void> => {
    await page.evaluate((args: {prop: string, value: any}) => {
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      tileV.props[args.prop] = args.value;
    }, {prop, value});
    await v.waitForViewerRendered(page, 'Tile Viewer', capMs);
  };

  await softStep('Lanes: set RACE → one lane per category, headers match RACE categories', async () => {
    await propRow('lanes');
    await setTileProp('lanesColumnName', 'RACE');
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
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
    await setTileProp('lanesColumnName', 'SEX');
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
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
    await setTileProp('lanesColumnName', null);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const lanes = Array.from(root.querySelectorAll('.d4-tile-viewer-lane'));
      return {
        laneCount: lanes.length,
        single: lanes[0]?.classList.contains('d4-tile-viewer-lane-single'),

        tilesInLane: lanes[0]?.querySelectorAll('.d4-tile-viewer-form').length ?? 0,
      };
    });
    expect(r.laneCount).toBe(1);
    expect(r.single).toBe(true);

    expect(r.tilesInLane).toBeGreaterThan(0);
  });

  await softStep('Row source: Selected → tiles show only the selected rows', async () => {

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (const i of [0, 1, 2, 3, 4]) df.selection.set(i, true);
    });
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Selected');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const selectedIdx = [0, 1, 2, 3, 4];

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

      const filtered = new Promise((res) => {
        const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 1500);
      });
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await filtered;
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

    expect(r.allTilesM).toBe(true);
  });

  await softStep('Row source: All → every row is shown again, filter cleared', async () => {

    const entry = await page.evaluate(() => {
      const tiles = Array.from(document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form'));
      return {beforeAllMale: tiles.length > 0 &&
        tiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M')};
    });
    await propRow('row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'All');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);

    const isolated = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        tiles: tiles.length,
        hasFemale: tiles.some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F'),
        filterCount: df.filter.trueCount,
        total: df.rowCount,
      };
    });

    expect(isolated.tiles).toBeGreaterThan(0);
    expect(isolated.hasFemale).toBe(true);

    expect(isolated.filterCount).toBeLessThan(isolated.total);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();

      const restored = new Promise((res) => {
        const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 1500);
      });
      for (const f of [...fg.filters]) fg.remove(f);
      df.selection.setAll(false);
      await restored;

      const afterTiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const afterHasFemale = afterTiles.some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F');
      return {afterHasFemale, filterCount: df.filter.trueCount, total: df.rowCount};
    });
    expect(entry.beforeAllMale).toBe(true);
    expect(r.afterHasFemale).toBe(true);

    expect(r.filterCount).toBe(r.total);
  });

  const fontRow = () => page.locator('.property-grid tr[name="prop-tiles-font"]');
  await softStep('Tiles font: size 18px grows the lane headers and tile text', async () => {

    await setTileProp('lanesColumnName', 'RACE');
    await propRow('tiles-font', 'style');
    const readHeader = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const h = root.querySelector('.d4-tile-viewer-lane-header') as HTMLElement;
      const tileInput = root.querySelector('.d4-tile-viewer-form input[name="input-AGE"]') as HTMLElement;
      return {
        size: h.style.fontSize, line: h.style.lineHeight,

        tile: getComputedStyle(tileInput).fontSize,
      };
    });
    await fontRow().locator('input.d4-font-size-input').fill('13');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const base = await readHeader();
    await fontRow().locator('input.d4-font-size-input').fill('18');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const grown = await readHeader();

    expect(base.size).toBe('13px');
    expect(base.line).toBe('18.2px');

    expect(grown.size).toBe('18px');
    expect(grown.line).toBe('25.2px');
    expect(grown.tile).toBe('18px');
  });

  await softStep('Tiles font: the family choice reaches the lane header and the tiles', async () => {
    await propRow('tiles-font', 'style');
    await fontRow().locator('select').selectOption('Arial');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const h = root.querySelector('.d4-tile-viewer-lane-header') as HTMLElement;
      const tileInput = root.querySelector('.d4-tile-viewer-form input[name="input-AGE"]') as HTMLElement;
      return {inlineFont: h.style.font, tileFamily: getComputedStyle(tileInput).fontFamily};
    });

    expect(r.inlineFont).toContain('Arial');

    expect(r.tileFamily).toContain('Arial');
  });

  await softStep('Tiles font: reset to default 13px Roboto restores the header font', async () => {
    await propRow('tiles-font', 'style');
    await fontRow().locator('select').selectOption('Roboto');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    await fontRow().locator('input.d4-font-size-input').fill('13');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;

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

  const AG_FIXTURE = 'ag-fixture';
  try {
    await softStep('Auto-generate (auto state): deleting a fielded column removes it and frees a slot for an excluded column', async () => {
      const entry = await page.evaluate(async (frameName: string) => {
        const t = grok.shell.tv.dataFrame.clone();
        t.name = frameName;
        const view = grok.shell.addTableView(t);
        const tileV = view.addViewer('Tile Viewer');

        const root = tileV.root as Element;
        root.setAttribute('data-ag-fixture', '1');

        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;

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

          autoGenerateTrue: tileV.props.autoGenerate === true,
          formNotDesigned: tileV.props.sketchState?.['formDesigned'] === false,
          totalCols: t.columns.length,
          renderedBefore,
          excludedBefore: t.columns.names().filter((c: string) => !renderedBefore.includes(c)),

          victim: renderedBefore.filter((c: string) => {
            const el = root.querySelector(`.d4-tile-viewer-form input[name="${slug(c)}"]`) as HTMLInputElement | null;
            return el != null && el.type !== 'checkbox' && !el.disabled;
          })[0] ?? null,
        };
      }, AG_FIXTURE);

      if (entry.victim != null)
        await removeColumnViaFieldMenu(entry.victim, '[data-ag-fixture="1"]');
      const after = await page.evaluate(async (args: {frameName: string, victim: string | null, excluded: string[]}) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === args.frameName);
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = tileV.root as Element;
        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;
        const victimField = args.victim == null ? null : slug(args.victim);

        const victimHosts = () => victimField == null ? 0
          : root.querySelectorAll(`.d4-tile-viewer-form input[name="${victimField}"]`).length;

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

          if (victimHosts() === 0 && refilled.length > 0) break;
          await new Promise((res) => setTimeout(res, 500));
        }

        const keeper = renderedAfter[0] ?? null;
        const tile = root.querySelector('.d4-tile-viewer-form');
        return {
          renderedAfter, refilled, keeper,
          keeperValue: keeper == null ? null
            : (tile?.querySelector(`input[name="${slug(keeper)}"]`) as HTMLInputElement)?.value,
          keeperDisplay: keeper == null ? null : t.col(keeper).getString(0),

          victimHostsOnTiles: victimHosts(),
          victimGoneEverywhere: Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
            .every((tl) => victimField == null || !tl.querySelector(`input[name="${victimField}"]`)),
          renderedTiles: root.querySelectorAll('.d4-tile-viewer-form').length,

          victimStillInFrame: args.victim != null && t.columns.names().includes(args.victim),
          totalColsAfter: t.columns.length,
        };
      }, {frameName: AG_FIXTURE, victim: entry.victim, excluded: entry.excludedBefore});
      const r = {...entry, ...after};
      expect(r.autoGenerateTrue).toBe(true);
      expect(r.formNotDesigned).toBe(true);

      expect(r.renderedBefore.length).toBe(10);
      expect(r.renderedBefore.length).toBeLessThan(r.totalCols);
      expect(r.excludedBefore.length).toBeGreaterThan(0);

      expect(r.victim).not.toBeNull();

      expect(r.renderedTiles).toBeGreaterThan(0);
      expect(r.victimHostsOnTiles).toBe(0);
      expect(r.victimGoneEverywhere).toBe(true);

      expect(r.victimStillInFrame).toBe(false);
      expect(r.totalColsAfter).toBe(r.totalCols - 1);

      expect(r.renderedAfter).not.toEqual(r.renderedBefore);

      expect(r.refilled.length).toBeGreaterThan(0);
      expect(r.renderedAfter.length).toBe(10);
      expect(r.keeperValue).toBe(r.keeperDisplay);
    });

    await softStep('Auto-generate (designed state): the same column delete does not refill the freed slot', async () => {

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

      await page.locator(`.grok-view-sketch .d4-host[name="div-${picked.designCol}"]`)
        .filter({has: page.locator(`input[name="input-${picked.designCol}"]`)}).first().click();
      await page.keyboard.press('Delete');

      await page.waitForTimeout(300);
      await page.locator('[name="button-CLOSE-AND-APPLY"]').click();
      await page.locator('.grok-view-sketch').waitFor({state: 'detached', timeout: 15000});
      await v.waitForViewerRendered(page, 'Tile Viewer', 900);

      const pre = await page.evaluate((args: {frameName: string, designCol: string, frameCol: string}) => {
        const view = Array.from(grok.shell.views).find((x: any) => x.name === args.frameName);
        grok.shell.v = view;
        const t = view.dataFrame;
        const tileV = view.viewers.find((x: any) => x.type === 'Tile Viewer');
        const root = tileV.root as Element;
        const slug = (c: string) => `input-${c.replace(/[^A-Za-z0-9]/g, '-')}`;

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

          excludedBefore: t.columns.names().filter((c: string) => !before.includes(slug(c))).map(slug),
          victimField: slug(args.frameCol),
        };
      }, {frameName: AG_FIXTURE, designCol: picked.designCol!, frameCol: picked.frameCol!});

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

      expect(r.after).not.toContain(r.victimField);
      expect(r.victimGoneEverywhere).toBe(true);

      expect(r.excludedBefore.length).toBeGreaterThan(0);
      expect(r.refilled).toEqual([]);

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

  await softStep('Multiple table switching: Table property → spgi-100 shows spgi-100 rows', async () => {
    await page.evaluate(async () => {
      const spgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      spgi.name = 'spgi-100';
      grok.shell.addTableView(spgi);

      await new Promise((res) => {
        const sub = spgi.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 2000);
      });
      const demogView = Array.from(grok.shell.views).find((x: any) => x.name === 'Table');
      if (demogView) grok.shell.v = demogView;

      await new Promise((res) => setTimeout(res, 300));
    });
    await propRow('table');
    await v.selectPropertyGridChoice(page, 'table', 'spgi-100');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
      const fields = Array.from(root.querySelectorAll('.d4-tile-viewer-form input[name^="input-"]'))
        .map((i) => i.getAttribute('name'));

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

    expect(r.shown).not.toBeNull();
    expect(r.tileValue).toBe(r.cellValue);
  });

  await softStep('Multiple table switching: Table property → demog shows demog rows again', async () => {
    await propRow('table');
    await v.selectPropertyGridChoice(page, 'table', 'Table');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const tileV = grok.shell.tv.viewers.find((x: any) => x.type === 'Tile Viewer');
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

  await softStep('Viewer title/description: title appears in the header, description below it', async () => {
    await propRow('title', 'description');
    await v.setPropertyGridValue(page, 'title', 'Patient Cards');
    await v.setPropertyGridValue(page, 'description', 'Demographic data per patient');

    await page.locator('[name="viewer-Tile-Viewer"] .d4-viewer-description')
      .filter({hasText: 'Demographic data per patient'}).first().waitFor({timeout: 5000});
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const panel = root.closest('.panel-base')!;
      return {
        titlebar: panel.querySelector('.panel-titlebar-text')?.textContent,
        descText: root.querySelector('.d4-viewer-description')?.textContent,
      };
    });
    expect(r.titlebar).toBe('Patient Cards');
    expect(r.descText).toBe('Demographic data per patient');
  });

  await softStep('Viewer title/description: Description Position Bottom moves it below the tiles', async () => {

    const order = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const desc = root.querySelector('.d4-viewer-description')!;
      const lanes = root.querySelector('.d4-tile-viewer-lanes-host')!;

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

  await softStep('Viewer title/description: clearing both leaves the header without a title', async () => {
    await propRow('title', 'description');

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

    await page.locator('[name="viewer-Tile-Viewer"] .d4-viewer-description')
      .waitFor({state: 'detached', timeout: 5000}).catch(() => {});
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      const panel = root.closest('.panel-base')!;
      return {

        titlebarPresent: !!panel.querySelector('.panel-titlebar-text'),
        titlebar: (panel.querySelector('.panel-titlebar-text')?.textContent ?? '').trim(),
        descPresent: !!root.querySelector('.d4-viewer-description'),
      };
    });
    expect(r.titlebarPresent).toBe(true);
    expect(r.titlebar).toBe('');
    expect(r.descPresent).toBe(false);
  });

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
      const filtered = new Promise((res) => {
        const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 1500);
      });
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await filtered;
      const tiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      return {
        total, afterSex: df.filter.trueCount,

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
      const filtered = new Promise((res) => {
        const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 1500);
      });
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 51, max: 200});
      await filtered;
      const ages = ageOnTiles();
      return {
        before, after: df.filter.trueCount,
        tilesBefore, tilesAfter: root.querySelectorAll('.d4-tile-viewer-form').length,

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

      const beforeTiles = Array.from(root.querySelectorAll('.d4-tile-viewer-form'));
      const beforeAllMale = beforeTiles.length > 0 &&
        beforeTiles.every((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'M');

      const restored = new Promise((res) => {
        const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 1500);
      });
      for (const f of [...fg.filters]) fg.remove(f);
      await restored;

      const hasFemale = () => Array.from(root.querySelectorAll('.d4-tile-viewer-form'))
        .some((t) => (t.querySelector('input[name="input-SEX"]') as HTMLInputElement)?.value === 'F');
      const deadline = Date.now() + 3000;
      while (!hasFemale() && Date.now() < deadline)
        await new Promise((res) => setTimeout(res, 100));
      const afterHasFemale = hasFemale();
      return {beforeAllMale, afterHasFemale, filterCount: df.filter.trueCount, total: df.rowCount};
    });

    expect(r.beforeAllMale).toBe(true);
    expect(r.afterHasFemale).toBe(true);
    expect(r.filterCount).toBe(r.total);
  });

  await softStep('Filter interaction: closing the filter panel removes it from the view', async () => {

    expect(await page.locator('[name="viewer-Filters"]').count()).toBeGreaterThan(0);

    await v.clickViewerTitlebarIcon(page, 'Filters', 'Close');
    await page.locator('[name="viewer-Filters"]').waitFor({state: 'detached', timeout: 10000});
    const r = await page.evaluate(() => ({
      panelGone: document.querySelector('[name="viewer-Filters"]') == null,
      tilesStillRendered: document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').length,
    }));
    expect(r.panelGone).toBe(true);
    expect(r.tilesStillRendered).toBeGreaterThan(0);
  });

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
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const after = await readAges();

    expect(before.ages.length).toBeGreaterThan(0);
    expect(before.ages.some((a) => a <= 50)).toBe(true);

    expect(after.ages.length).toBeGreaterThan(0);
    expect(after.ages.every((a) => a > 50)).toBe(true);

    expect(after.filterCount).toBe(before.filterCount);
    expect(after.filterCount).toBe(after.total);

    const row2 = await propRow('filter');
    const box2 = row2.locator('input.property-grid-ellipsis-editor-input').first();
    await box2.fill('');
    await box2.press('Enter');
    await v.waitForViewerRendered(page, 'Tile Viewer', 900);
    const cleared = await readAges();
    expect(cleared.ages.some((a) => a <= 50)).toBe(true);
    expect(cleared.filterCount).toBe(cleared.total);
  });

  await softStep('Scroll position: a scrolled lane keeps its position and its rows when another viewer is added', async () => {

    await setTileProp('lanesColumnName', null);
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

      const added = new Promise((res) => {
        const sub = grok.events.onViewerAdded.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 2500);
      });
      grok.shell.tv.addViewer('Histogram');
      await added;

      await new Promise((res) => setTimeout(res, 2000));
    });
    const after = await readLane();

    expect(before.scrollTop).toBeGreaterThan(0);
    expect(before.age).not.toBeNull();

    expect(after.scrollTop).toBeGreaterThan(0);
    expect(Math.abs(after.scrollTop - before.scrollTop)).toBeLessThanOrEqual(2);
    expect(after.age).toBe(before.age);
    expect(after.sex).toBe(before.sex);
    expect(after.weight).toBe(before.weight);
  });

  v.finishSpec();
});
