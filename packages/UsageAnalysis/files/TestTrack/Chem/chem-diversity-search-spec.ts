/* ---
realizes: [chem.cp.diversity-search-viewer]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {openChemMenuItem} from '../helpers/chem';
import {finishSpec, ensurePropertyCategory, openViewerGear, propertyGridValue,
  selectPropertyGridChoice} from '../helpers/viewers';
import {addCardViaColumnSelector, cardCaptions} from '../helpers/filter-panel';

declare const grok: any;

// Selector recon-notes (grok-browser references, read 2026-08-22):
//   [name="viewer-Chem-Diversity-Search"], .chem-diversity-search, .chem-similarity-header
//     — references/chem.md:1513-1515; the header class is SHARED with Chem Similarity Search
//     (chem.md:1534-1540), hence the container scoping of HEADER_SEL.
//   Property rows prop-row-source / prop-fingerprint / prop-limit / prop-distance-metric /
//     prop-size and their prop-view- value cells — chem.md:1604-1616; the Diversity grid carries
//     six knob rows and no cutoff / follow-current-row — chem.md:1517-1525.
//   [name="prop-category-misc"] must be clicked before any Misc row is reachable —
//     chem.md:1527-1532 (this viewer) and chem.md:1629-1652 (the general recipe).
//   Enum rows: click the prop-view- value cell, then drive the select the click creates —
//     chem.md:1693-1711. Numeric rows: [name="prop-<name>"] input.property-grid-slider-textbox
//     (click / fill / Enter); the value cell itself is 246x0 — chem.md:1659-1667.
//   .panel-titlebar [name="icon-font-icon-settings"] (gear) — references/viewers/charts.md:9,:71.
//   [name="icon-filter"] in .d4-ribbon-panel opens the Filter Panel —
//     references/viewers/filters.md:52, :608.
//   [name="viewer-Filters"] .d4-filter and .d4-filter-column-name — filters.md:142-152;
//     [name="viewer-Histogram"] inside a card marks it numeric — filters.md:861;
//     [name="icon-list"] switches it to categorical — filters.md:862, :304;
//     [name="icon-search"], input.d4-search-input[placeholder="Search"] and the SECOND
//     input[type="checkbox"].ui-input-editor ("filter by search results", on by default) —
//     filters.md:831-833, with the search-narrows-rows measurement at filters.md:835.
//     Per-card icons are hover-gated and are clicked on the element itself — filters.md:1594, :1643.
//   .d4-balloon.error — the viewer's failure channel is grok.shell.error
//     (public/packages/Chem/src/analysis/chem-diversity-viewer.ts:53-55).

const MOL_COL = 'canonical_smiles';
const DIVERSITY_LEAF = 'Diversity Search...';
const VIEWER_TYPE = 'Chem Diversity Search';
const VIEWER_NAME = 'Chem-Diversity-Search';
const VIEWER_SEL = '[name="viewer-Chem-Diversity-Search"]';
const CARDS_SEL = '.chem-diversity-search';
const HEADER_SEL = `${VIEWER_SEL} .chem-similarity-header`;
const LIMIT_AFTER_EDIT = 6;
const FILTERED_ROWS = 4;

async function readViewer(page: Page): Promise<any> {
  return page.evaluate((t) => {
    const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
    if (!v) return null;
    return {
      ids: Array.from(v.renderMolIds) as number[],
      idsLen: v.renderMolIds.length,
      limit: v.limit,
      distanceMetric: v.distanceMetric,
      fingerprint: v.fingerprint,
      size: v.size,
      rowSource: v.rowSource,
    };
  }, VIEWER_TYPE);
}

async function cardCount(page: Page): Promise<number> {
  return page.locator(`${CARDS_SEL} > div`).count();
}

// The gear is clicked only when the panel is showing some OTHER object's properties: it toggles,
// so a blind click on an already-open Diversity panel would close it.
async function openDiversityProperties(page: Page, prop: string, category: string): Promise<void> {
  const rows = () => page.locator(
    '.property-grid tr[name="prop-distance-metric"], .property-grid tr[name="prop-row-source"]');
  if (await rows().count() < 2) {
    await openViewerGear(page, VIEWER_TYPE);
    await page.locator('.property-grid tr[name="prop-distance-metric"]').first()
      .waitFor({state: 'attached', timeout: 20_000});
  }
  await ensurePropertyCategory(page, VIEWER_NAME, category, prop);
}

async function setLimitViaPropertyPanel(page: Page, value: number): Promise<void> {
  await openDiversityProperties(page, 'limit', 'misc');
  const box = page.locator('[name="prop-limit"] input.property-grid-slider-textbox').first();
  await box.click();
  await box.fill(String(value));
  await box.press('Enter');
}

function sameIdSet(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  const x = [...a].sort((p, q) => p - q);
  const y = [...b].sort((p, q) => p - q);
  return x.every((v, i) => v === y[i]);
}

async function waitIdSetChanged(page: Page, prevIds: number[]): Promise<void> {
  await page.waitForFunction(({t, prev}) => {
    const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
    if (!v) return false;
    const cur = Array.from(v.renderMolIds) as number[];
    if (cur.length !== prev.length) return true;
    const x = [...cur].sort((p: number, q: number) => p - q);
    const y = [...prev].sort((p: number, q: number) => p - q);
    return x.some((val: number, i: number) => val !== y[i]);
  }, {t: VIEWER_TYPE, prev: prevIds}, {timeout: 60_000});
}

async function cardCanvasSizes(page: Page): Promise<{w: number; h: number}[]> {
  return page.evaluate((sel) => Array.from(document.querySelector(sel)?.children ?? [])
    .map((child) => {
      const canvas = child.querySelector('canvas') as HTMLCanvasElement | null;
      return {w: parseFloat(canvas?.style.width ?? '') || 0, h: parseFloat(canvas?.style.height ?? '') || 0};
    }), CARDS_SEL);
}

function distinctSizes(sizes: {w: number; h: number}[]): string[] {
  return [...new Set(sizes.map((s) => `${s.w}x${s.h}`))];
}

async function errorBalloons(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon.error'))
    .map((b) => (b.textContent ?? '').trim()));
}

interface FilterTarget {column: string; category: string; distinct: number}

async function pickFilterTarget(page: Page, wantRows: number, molCol: string): Promise<FilterTarget | null> {
  return page.evaluate(({want, mol}) => {
    const df = grok.shell.t;
    let best: any = null;
    for (const c of df.columns.toList()) {
      if (c.name === mol || (c.type !== 'int' && c.type !== 'string')) continue;
      const counts: Record<string, number> = {};
      for (let i = 0; i < c.length; i++) {
        const key = String(c.get(i));
        counts[key] = (counts[key] ?? 0) + 1;
      }
      const labels = Object.keys(counts);
      if (labels.length > 30) continue;
      for (const label of labels) {
        if (label === '' || counts[label] !== want) continue;
        if (labels.some((other) => other !== label && other.includes(label))) continue;
        if (!best || labels.length < best.distinct)
          best = {column: c.name, category: label, distinct: labels.length};
      }
    }
    return best;
  }, {want: wantRows, mol: molCol});
}

async function trueCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.t.filter.trueCount as number);
}

async function cardIsHistogram(page: Page, column: string): Promise<boolean> {
  return page.evaluate((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name') as HTMLElement | null)?.textContent ?? '')
        .trim() === col);
    return !!card?.querySelector('[name="viewer-Histogram"]');
  }, column);
}

async function switchCardToCategorical(page: Page, column: string): Promise<void> {
  await page.evaluate((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name') as HTMLElement | null)?.textContent ?? '')
        .trim() === col);
    (card?.querySelector('[name="icon-list"]') as HTMLElement | null)?.click();
  }, column);
  await page.waitForFunction((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name') as HTMLElement | null)?.textContent ?? '')
        .trim() === col);
    return !!card?.querySelector('[name="viewer-Grid"]');
  }, column, {timeout: 20_000});
}

function filterCard(page: Page, column: string) {
  return page.locator('[name="viewer-Filters"] .d4-filter')
    .filter({has: page.locator('.d4-filter-column-name', {hasText: new RegExp(`^${column}$`)})}).first();
}

test.use(specTestOptions);

test('Chem: Diversity Search viewer — metric / limit / fingerprint / size / row-source controls', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);

  let total = 0;
  let idsBeforeMetric: number[] = [];
  let idsBeforeFingerprint: number[] = [];

  await softStep('Setup: open smiles.csv, canonical_smiles Molecule column ready', async () => {
    await page.evaluate(async () => {
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
    await waitForMolecule(page);
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30_000});
    const info = await page.evaluate((molCol) => {
      const df = grok.shell.t;
      return {
        total: df.rowCount as number,
        hasMolCol: df.columns.toList().some((c: any) => c.name === molCol && c.semType === 'Molecule'),
      };
    }, MOL_COL);
    total = info.total;
    console.log(`[probe] setup: rowCount=${info.total} hasMolCol=${info.hasMolCol}`);
    expect(info.hasMolCol, `${MOL_COL} must be a Molecule column`).toBe(true);
    expect(total).toBeGreaterThan(12);
  });

  await softStep('Scenario 1, Step 2: top-menu Chem | Search | Diversity Search adds the viewer', async () => {
    await openChemMenuItem(page, DIVERSITY_LEAF);
    await page.locator(VIEWER_SEL).first().waitFor({state: 'visible', timeout: 30_000});
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      return v && v.renderMolIds.length === v.limit;
    }, VIEWER_TYPE, {timeout: 60_000});
    const state = await readViewer(page);
    console.log(`[probe] viewer added: limit=${state?.limit} idsLen=${state?.idsLen} ` +
      `metric=${state?.distanceMetric} fp=${state?.fingerprint} size=${state?.size} rowSource=${state?.rowSource}`);
    expect(state, 'the Diversity Search viewer must be present in the table view').not.toBeNull();
    expect(state.limit, 'default limit is 12').toBe(12);
  });

  await softStep('Scenario 1 Step 4: the viewer renders one card per diverse structure at the current limit', async () => {
    const state = await readViewer(page);
    const cards = await cardCount(page);
    console.log(`[probe] S1S4: idsLen=${state.idsLen} limit=${state.limit} cards=${cards}`);
    expect(state.idsLen, 'the diverse subset holds `limit` molecule ids').toBe(state.limit);
    expect(cards, 'the viewer renders exactly one molecule card per diverse id').toBe(state.limit);
    const header = await page.locator(HEADER_SEL).first().textContent();
    expect(header ?? '', 'the header reflects the active metric + fingerprint').toContain('Tanimoto, Morgan');
    idsBeforeMetric = state.ids;
    expect(idsBeforeMetric.length, 'the initial diverse id set is captured for the metric contrast').toBe(state.limit);
  });

  await softStep('Scenario 1 Step 5: changing the Distance Metric re-runs diversity search cleanly at full limit', async () => {
    await openDiversityProperties(page, 'distance-metric', 'misc');
    const metricBefore = await propertyGridValue(page, 'distance-metric', 'misc');
    console.log(`[probe] metric: property panel shows distance-metric=${JSON.stringify(metricBefore)} before the edit`);
    expect(metricBefore, 'the property panel must show the Distance Metric row before it is driven, ' +
      'or the edit below lands on nothing').toBe('Tanimoto');
    await selectPropertyGridChoice(page, 'distance-metric', 'Cosine', 'misc');
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      return v && v.distanceMetric === 'Cosine';
    }, VIEWER_TYPE, {timeout: 30_000});
    await waitIdSetChanged(page, idsBeforeMetric);
    const header = await page.locator(HEADER_SEL).first().textContent();
    expect(header ?? '', 'the header reflects the newly applied metric').toContain('Cosine, Morgan');
    const cards = await cardCount(page);
    const state = await readViewer(page);
    const balloons = await errorBalloons(page);
    const sameSet = sameIdSet(state.ids, idsBeforeMetric);
    console.log(`[probe] metric: idsLen=${state.idsLen} limit=${state.limit} cards=${cards} ` +
      `sameSet=${sameSet} errorBalloons=${JSON.stringify(balloons)} ` +
      `before=${JSON.stringify(idsBeforeMetric.slice(0, 12))} after=${JSON.stringify(state.ids.slice(0, 12))}`);
    expect(state.idsLen, 'the re-run produces a full diverse subset at the current limit').toBe(state.limit);
    expect(cards, 'the viewer renders one molecule card per diverse id after the metric change').toBe(state.limit);
    expect(sameSet,
      `a new distance metric must re-query into a different diverse id SET, not merely reorder the same ` +
      `molecules; before=${JSON.stringify(idsBeforeMetric)} after=${JSON.stringify(state.ids)}`).toBe(false);
    expect(balloons,
      `the metric re-run must raise no error balloon — read instantaneously, since balloons auto-hide ` +
      `after 5 s; got ${JSON.stringify(balloons)}`).toEqual([]);
  });

  await softStep('Scenario 1 Step 6: setting the limit to 6 shows exactly 6 cards', async () => {
    await setLimitViaPropertyPanel(page, LIMIT_AFTER_EDIT);
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      return v && v.limit === 6 && v.renderMolIds.length === 6;
    }, VIEWER_TYPE, {timeout: 60_000});
    await page.waitForFunction((sel) => document.querySelectorAll(`${sel} > div`).length === 6, CARDS_SEL, {timeout: 30_000});
    const cards = await cardCount(page);
    const state = await readViewer(page);
    console.log(`[probe] limit=6: limit=${state.limit} idsLen=${state.idsLen} cards=${cards}`);
    expect(state.limit, 'the limit edit set the model limit to 6').toBe(LIMIT_AFTER_EDIT);
    expect(state.idsLen, 'the diverse subset now holds exactly 6 ids').toBe(LIMIT_AFTER_EDIT);
    expect(cards, 'the limit edit re-renders the card set to exactly 6 cards').toBe(LIMIT_AFTER_EDIT);
  });

  await softStep('Scenario 2 Step 3: changing the Fingerprint re-runs diversity search cleanly at full limit', async () => {
    idsBeforeFingerprint = (await readViewer(page)).ids;
    await openDiversityProperties(page, 'fingerprint', 'misc');
    await selectPropertyGridChoice(page, 'fingerprint', 'MACCS', 'misc');
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      return v && v.fingerprint === 'MACCS';
    }, VIEWER_TYPE, {timeout: 30_000});
    await waitIdSetChanged(page, idsBeforeFingerprint);
    const header = await page.locator(HEADER_SEL).first().textContent();
    expect(header ?? '', 'the header reflects the newly applied fingerprint').toContain('MACCS');
    const cards = await cardCount(page);
    const state = await readViewer(page);
    const balloons = await errorBalloons(page);
    const sameSet = sameIdSet(state.ids, idsBeforeFingerprint);
    console.log(`[probe] fingerprint: idsLen=${state.idsLen} limit=${state.limit} cards=${cards} ` +
      `sameSet=${sameSet} errorBalloons=${JSON.stringify(balloons)} ` +
      `before=${JSON.stringify(idsBeforeFingerprint.slice(0, 12))} after=${JSON.stringify(state.ids.slice(0, 12))}`);
    expect(state.idsLen, 'the re-run produces a full diverse subset at the current limit').toBe(state.limit);
    expect(cards, 'the viewer renders one molecule card per diverse id after the fingerprint change').toBe(state.limit);
    expect(sameSet,
      `a new fingerprint must re-query into a different diverse id SET, not merely reorder the same ` +
      `molecules; before=${JSON.stringify(idsBeforeFingerprint)} after=${JSON.stringify(state.ids)}`).toBe(false);
    expect(balloons,
      `the fingerprint re-run must raise no error balloon — read instantaneously, since balloons ` +
      `auto-hide after 5 s; got ${JSON.stringify(balloons)}`).toEqual([]);
  });

  await softStep('Scenario 2 Step 4: changing Size resizes each molecule card', async () => {
    await openDiversityProperties(page, 'size', 'misc');
    await selectPropertyGridChoice(page, 'size', 'normal', 'misc');
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      return v && v.size === 'normal';
    }, VIEWER_TYPE, {timeout: 30_000});
    await page.waitForFunction((sel) => {
      const c = document.querySelector(sel)?.children[0]?.querySelector('canvas') as HTMLCanvasElement | null;
      return c != null && parseFloat(c.style.width) === 200 && parseFloat(c.style.height) === 100;
    }, CARDS_SEL, {timeout: 30_000});
    const normalSizes = await cardCanvasSizes(page);
    console.log(`[probe] size=normal: cards=${normalSizes.length} distinct=${JSON.stringify(distinctSizes(normalSizes))}`);
    expect(normalSizes.length, 'every card must be measured, or the size claims below are vacuous')
      .toBe(LIMIT_AFTER_EDIT);
    expect(distinctSizes(normalSizes),
      `every card must carry the normal tile size, got ${JSON.stringify(distinctSizes(normalSizes))}`)
      .toEqual(['200x100']);
    const normalSize = normalSizes[0];

    await selectPropertyGridChoice(page, 'size', 'large', 'misc');
    await page.waitForFunction((sel) => {
      const c = document.querySelector(sel)?.children[0]?.querySelector('canvas') as HTMLCanvasElement | null;
      return c != null && parseFloat(c.style.width) === 300 && parseFloat(c.style.height) === 150;
    }, CARDS_SEL, {timeout: 30_000});
    const largeSizes = await cardCanvasSizes(page);
    console.log(`[probe] size=large: cards=${largeSizes.length} distinct=${JSON.stringify(distinctSizes(largeSizes))} ` +
      `(normal was ${JSON.stringify(distinctSizes(normalSizes))})`);
    expect(largeSizes.length, 'the card set must survive the resize').toBe(normalSizes.length);
    expect(distinctSizes(largeSizes),
      `every card must carry the large tile size, got ${JSON.stringify(distinctSizes(largeSizes))}`)
      .toEqual(['300x150']);
    const largeSize = largeSizes[0];
    expect(largeSize.w, 'large tiles are wider than normal').toBeGreaterThan(normalSize.w);
    expect(largeSize.h, 'large tiles are taller than normal').toBeGreaterThan(normalSize.h);
  });

  await softStep('Scenario 2 Step 5: Row Source = Filtered restricts the diverse subset to the filtered rows', async () => {
    const target = await pickFilterTarget(page, FILTERED_ROWS, MOL_COL);
    console.log(`[probe] filter target: ${JSON.stringify(target)} (want ${FILTERED_ROWS} rows)`);
    expect(target, `no low-cardinality column in smiles.csv has a category holding exactly ` +
      `${FILTERED_ROWS} rows whose label is not a substring of a sibling label, so the Filter Panel ` +
      'cannot be narrowed below the limit through the card search box').not.toBeNull();

    await page.locator('.d4-ribbon-panel [name="icon-filter"]').first().click();
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 30_000});
    if (!(await cardCaptions(page)).includes(target!.column))
      await addCardViaColumnSelector(page, target!.column);
    if (await cardIsHistogram(page, target!.column))
      await switchCardToCategorical(page, target!.column);

    const card = filterCard(page, target!.column);
    await card.locator('[name="icon-search"]').first().evaluate((el: HTMLElement) => el.click());
    const search = card.locator('input.d4-search-input[placeholder="Search"]').first();
    await search.waitFor({state: 'visible', timeout: 20_000});
    const filterBySearch = card.locator('input[type="checkbox"].ui-input-editor').nth(1);
    const filterBySearchOn = await filterBySearch.isChecked();
    console.log(`[probe] card ${target!.column}: filter-by-search checkbox checked=${filterBySearchOn}`);
    expect(filterBySearchOn, 'the card\'s "filter by search results" box must be on, or typing in the ' +
      'search box narrows nothing and the row source below has no filtered set to clamp to').toBe(true);
    await search.fill(target!.category);

    await expect.poll(() => trueCount(page), {
      timeout: 30_000,
      intervals: [250, 500, 1000],
      message: `the ${target!.column} card search for ${JSON.stringify(target!.category)} never narrowed the ` +
        `table to ${FILTERED_ROWS} rows`,
    }).toBe(FILTERED_ROWS);
    const filteredCount = await trueCount(page);
    console.log(`[probe] filter panel: column=${target!.column} category=${target!.category} ` +
      `trueCount=${filteredCount} of ${total}`);
    expect(filteredCount, 'the Filter Panel narrowed the table below the viewer limit, so the clamp ' +
      'below discriminates a viewer that honours rowSource from one that ignores it').toBe(FILTERED_ROWS);

    await openDiversityProperties(page, 'row-source', 'data');
    await selectPropertyGridChoice(page, 'row-source', 'Filtered', 'data');
    await page.waitForFunction((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      const df = grok.shell.t;
      const expected = Math.min(df.filter.trueCount, v.limit);
      return v && v.rowSource === 'Filtered' && v.renderMolIds.length === expected;
    }, VIEWER_TYPE, {timeout: 60_000});

    const probe = await page.evaluate((t) => {
      const v = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === t) as any;
      const df = grok.shell.t;
      return {
        idsLen: v.renderMolIds.length,
        limit: v.limit,
        filteredCount: df.filter.trueCount as number,
        allInFilter: v.renderMolIds.every((id: number) => df.filter.get(id)),
      };
    }, VIEWER_TYPE);
    const cards = await cardCount(page);
    console.log(`[probe] rowSource=Filtered: filtered=${probe.filteredCount} limit=${probe.limit} ` +
      `idsLen=${probe.idsLen} cards=${cards} allInFilter=${probe.allInFilter}`);
    expect(probe.filteredCount, 'the filtered row set must still be the one the Filter Panel produced')
      .toBe(FILTERED_ROWS);
    expect(probe.limit, 'the limit is still the one set in Scenario 1 Step 6, above the filtered row count')
      .toBe(LIMIT_AFTER_EDIT);
    expect(probe.idsLen, 'the diverse subset holds every filtered row and no more — a viewer that ' +
      'ignored rowSource would fill the limit instead').toBe(FILTERED_ROWS);
    expect(cards, 'the viewer renders one molecule card per filtered diverse id').toBe(FILTERED_ROWS);
    expect(probe.allInFilter, 'every displayed molecule id comes from the filtered row set').toBe(true);
  });

  await softStep('Teardown: no error balloon left by diversity search; close views', async () => {
    const balloons = await errorBalloons(page);
    console.log(`[probe] teardown: error balloons=${JSON.stringify(balloons)}`);
    expect(balloons, `error balloons during diversity search: ${JSON.stringify(balloons)}`).toEqual([]);
    await page.evaluate(() => grok.shell.closeAll());
  });

  finishSpec();
});
