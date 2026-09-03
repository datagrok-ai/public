/* ---
realizes: [filters.cp.hierarchical-and-combined-boolean]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addHierarchicalCard, applyBoolState, applyHierarchyState, closeFilterPanelInPage, GLYPH, glyphName,
  hierCaption, hierNode, openDemogWithSexBool, ROW_COUNT, trueCountOf} from './hierarchical-shared';

declare const grok: any;
declare const window: any;

// The tree and combined-boolean gestures on the local lane. The layout and project round-trips
// (Steps 10, 11, 17, 18) live in hierarchical-and-combined-boolean-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function readTreeState(page: any): Promise<{captions: string[], trueCount: number}> {
  return await page.evaluate(() => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
    const captions: string[] = [];
    for (const node of card!.querySelectorAll('.d4-tree-view-node')) {
      const val = node.querySelector('.d4-hierarchical-filter-caption-value');
      if (val && (node as HTMLElement).offsetParent !== null) captions.push(val.textContent!.trim());
    }
    return {captions, trueCount: grok.shell.tv.dataFrame.filter.trueCount};
  });
}

async function holdTrueCount(page: any, expected: number, why: string): Promise<void> {
  const samples: number[] = [];
  for (let i = 0; i < 7; i++) {
    if (i > 0) await page.waitForTimeout(400);
    samples.push(await trueCountOf(page));
  }
  expect(Array.from(new Set(samples)),
    `${why} — the row filter had to read ${expected} at every sample of a 2.4s window, so a regression `
    + 'that moves the filter on a longer debounce, or moves it and self-corrects, cannot slip between '
    + `two reads; samples: ${samples.join(', ')}`)
    .toEqual([expected]);
}

async function typeTreeSearch(page: any, text: string): Promise<string> {
  return await page.evaluate((frag: string) => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
    const input = card?.querySelector('input.d4-search-input[placeholder="Search..."]') as HTMLInputElement;
    if (!input) throw new Error('hierarchical tree-search input not found');
    input.focus();
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, frag);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    return input.value;
  }, text);
}

test('Filter Panel — Hierarchical and Combined Boolean Filters', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBe(ROW_COUNT);

  const derived = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const sex = df.col('SEX').toList();
    const race = df.col('RACE').toList();
    let allF = 0;
    let caucasianFemale = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (sex[i] !== 'F') continue;
      allF++;
      if (race[i] === 'Caucasian') caucasianFemale++;
    }
    return {allF, caucasianFemale};
  });
  const otherFemale = derived.allF - derived.caucasianFemale;
  expect(derived.caucasianFemale,
    `no row is female and Caucasian (derived ${derived.caucasianFemale}) — checking F / Caucasian would `
    + 'narrow to nothing and the count expectations below would hold without the product doing anything')
    .toBeGreaterThan(0);
  expect(derived.caucasianFemale,
    `the female / Caucasian rows (${derived.caucasianFemale}) cover every female row (${derived.allF}) — `
    + 'the F branch and its Caucasian child would then be the same criterion and Step 7 would lose its contrast')
    .toBeLessThan(derived.allF);
  expect(derived.allF,
    `the derived female row count ${derived.allF} covers the whole ${rowCount}-row table — the F branch would not be a narrowing`)
    .toBeLessThan(rowCount);

  await softStep('Step 5: Add hierarchical filter (SEX / RACE), caption reflects columns', async () => {
    await addHierarchicalCard(page);
    expect(await hierCaption(page)).toBe('SEX / RACE');
    const trueCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(trueCount).toBe(ROW_COUNT);
  });

  await softStep('Step 6: Expand F, check Caucasian — Caucasian females filtered', async () => {
    await hierNode(page, ['F'], 'expand');
    await expect.poll(async () => (await hierNode(page, ['F'], 'read')).childCaptions,
      {message: 'the RACE children of F never rendered after the expander click',
        timeout: 10_000, intervals: [200, 400, 800]}).toContain('Caucasian');
    await hierNode(page, ['F', 'Caucasian'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'checking F / Caucasian did not narrow the table to the Caucasian-female rows derived '
        + `from the raw SEX / RACE columns (${derived.caucasianFemale})`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(derived.caucasianFemale);
    const state = await hierNode(page, ['F', 'Caucasian'], 'read');
    expect(state.glyph,
      'a checked leaf must render the checked glyph U+F14A, not the unchecked or indeterminate one')
      .toBe(GLYPH.checked);
    const trueCount = await trueCountOf(page);
    expect(trueCount, 'a checked criterion must narrow the table').toBeLessThan(rowCount);
    expect(trueCount, 'the criterion must not empty the table').toBeGreaterThan(0);
  });

  await softStep('Step 7: Uncheck one child of F — F parent reads indeterminate', async () => {
    await hierNode(page, ['F'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'checking the whole F branch did not select every female row derived from the raw SEX '
        + `column (${derived.allF})`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(derived.allF);
    expect((await hierNode(page, ['F'], 'read')).glyph,
      'a fully checked branch must read checked (U+F14A) before one of its children is unchecked')
      .toBe(GLYPH.checked);

    await hierNode(page, ['F', 'Caucasian'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'unchecking F / Caucasian did not drop the Caucasian-female rows: expected the derived '
        + `${derived.allF} female rows minus the derived ${derived.caucasianFemale} Caucasian-female ones = ${otherFemale}`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(otherFemale);
    const state = await hierNode(page, ['F'], 'read');
    expect(state.glyph,
      'the F parent must read indeterminate (U+F146) after one of its RACE children is unchecked')
      .toBe(GLYPH.indeterminate);
    expect(state.glyph,
      'F must not read unchecked — the branch still holds checked children')
      .not.toBe(GLYPH.unchecked);
    const trueCount = await trueCountOf(page);
    expect(trueCount, 'the partially checked branch must still select rows').toBeGreaterThan(0);
    expect(trueCount, 'the partially checked branch must select fewer rows than the full branch')
      .toBeLessThan(derived.allF);
  });

  await softStep('Step 8: GROK-19968 — tree search hides nodes without moving trueCount', async () => {
    const before = await readTreeState(page);
    expect(before.trueCount).toBeGreaterThan(0);
    expect(before.trueCount).toBeLessThan(rowCount);
    expect(before.captions).toContain('Caucasian');
    expect(before.captions).toContain('Asian');

    await page.evaluate(() => {
      const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
        .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
      if (!card)
        throw new Error('no hierarchical filter card is painted in the Filter Panel — its caption is the only '
          + 'one carrying "/", and no card in the panel carries one, so the tree search cannot be opened');
      const icon = card.querySelector('.d4-filter-header [name="icon-search"]') as HTMLElement | null;
      if (!icon)
        throw new Error('the hierarchical filter card\'s header carries no [name="icon-search"] icon, so the '
          + 'tree search was never opened through the card\'s own control');
      icon.click();
    });
    const searchBox = page
      .locator('[name="viewer-Filters"] .d4-filter input.d4-search-input[placeholder="Search..."]').first();
    await searchBox.waitFor({state: 'visible', timeout: 10000});
    expect(await typeTreeSearch(page, 'Cau')).toBe('Cau');
    await expect.poll(async () => (await readTreeState(page)).captions,
      {message: 'typing "Cau" never hid the non-matching tree nodes — the search term reached the input '
        + 'but the tree was never re-rendered',
      timeout: 15_000, intervals: [200, 400, 800]}).not.toContain('Asian');

    const after = await readTreeState(page);
    expect(after.captions).toContain('Caucasian');
    expect(after.captions).not.toContain('Asian');
    expect(after.captions.length).toBeLessThan(before.captions.length);
    await holdTrueCount(page, before.trueCount,
      'GROK-19968 — the tree search hides nodes and must never move the row filter');

    expect(await typeTreeSearch(page, '')).toBe('');
    await expect.poll(async () => (await readTreeState(page)).captions,
      {message: 'clearing the tree search did not bring back exactly the pre-search visible node set',
        timeout: 15_000, intervals: [200, 400, 800]}).toEqual(before.captions);
    await holdTrueCount(page, before.trueCount,
      'clearing the tree search must leave the row filter exactly where the search found it');
  });

  await softStep('Step 8b: three levels — unchecking a grandchild marks BOTH ancestors indeterminate', async () => {
    const levels = ['SEX', 'RACE', 'SEVERITY'];
    await applyHierarchyState(page, {colNames: levels, allEnabled: true});
    await expect.poll(() => hierCaption(page),
      {message: 'the card caption must name all three levels',
        timeout: 15_000, intervals: [300, 600, 1200]}).toBe(levels.join(' / '));
    await expect.poll(async () => (await hierNode(page, ['F'], 'probe')).found,
      {message: 'the SEX roots never rebuilt after the three-column hierarchy was applied',
        timeout: 15_000, intervals: [300, 600, 1200]}).toBe(true);

    const grandchild = 'None';
    const expected = await page.evaluate((severity: string) => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX').toList();
      const race = df.col('RACE').toList();
      const sev = df.col('SEVERITY').toList();
      let femaleRows = 0;
      let femaleCaucasianGrandchild = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (sex[i] !== 'F') continue;
        femaleRows++;
        if (race[i] === 'Caucasian' && sev[i] === severity) femaleCaucasianGrandchild++;
      }
      return {femaleRows, femaleCaucasianGrandchild};
    }, grandchild);
    const expectedSeverities: string[] = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX').toList();
      const race = df.col('RACE').toList();
      const sev = df.col('SEVERITY').toList();
      const values = new Set<string>();
      for (let i = 0; i < df.rowCount; i++)
        if (sex[i] === 'F' && race[i] === 'Caucasian') values.add(String(sev[i] ?? '').trim());
      return [...values].sort();
    });
    expect(expected.femaleRows,
      `derived female row count is ${expected.femaleRows} — it must be a real subset of the ${rowCount} table rows`)
      .toBeGreaterThan(0);
    expect(expected.femaleRows,
      `derived female row count ${expected.femaleRows} covers the whole table — the F branch would not be a narrowing`)
      .toBeLessThan(rowCount);
    expect(expected.femaleCaucasianGrandchild,
      `no row is female / Caucasian / ${grandchild} — unchecking that grandchild would drop nothing ` +
      'and the post-uncheck count assertion would hold without the product doing anything')
      .toBeGreaterThan(0);
    expect(expected.femaleCaucasianGrandchild,
      `female / Caucasian / ${grandchild} covers all ${expected.femaleRows} female rows — the uncheck ` +
      'would empty the branch and the count assertion would lose its contrast')
      .toBeLessThan(expected.femaleRows);
    expect(expectedSeverities,
      `only ${JSON.stringify(expectedSeverities)} occurs under female / Caucasian — with fewer than ` +
      `two values there is no sibling left once "${grandchild}" is unchecked and the locality claim is vacuous`)
      .toContain(grandchild);
    expect(expectedSeverities.length,
      `the derived SEVERITY set ${JSON.stringify(expectedSeverities)} must hold more than one value`)
      .toBeGreaterThan(1);
    expect(expectedSeverities.filter((s) => s === ''),
      `SEVERITY holds a blank value under female / Caucasian (${JSON.stringify(expectedSeverities)}) — ` +
      'a blank category cannot be addressed by caption in the tree')
      .toEqual([]);

    await hierNode(page, ['F'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'checking the F branch under the three-level hierarchy did not select every female row ' +
        `(expected the ${expected.femaleRows} rows with SEX = F)`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(expected.femaleRows);
    const allF = await trueCountOf(page);
    expect((await hierNode(page, ['F'], 'read')).glyph,
      'a fully checked branch reads checked (U+F14A), not partial').toBe(GLYPH.checked);

    await hierNode(page, ['F'], 'expand');
    await expect.poll(async () => (await hierNode(page, ['F'], 'read')).childCaptions,
      {message: 'the RACE level under F never rendered',
        timeout: 10_000, intervals: [200, 400, 800]}).toContain('Caucasian');
    await hierNode(page, ['F', 'Caucasian'], 'expand');
    await expect.poll(async () => (await hierNode(page, ['F', 'Caucasian'], 'read')).childCaptions.length,
      {message: 'the SEVERITY level under F / Caucasian never rendered — the third level is missing',
        timeout: 10_000, intervals: [200, 400, 800]}).toBeGreaterThan(0);

    await expect.poll(async () =>
      (await hierNode(page, ['F', 'Caucasian'], 'read')).childCaptions.slice().sort(),
    {message: 'the SEVERITY values listed under F / Caucasian must be exactly the ones the data holds ' +
      `for that branch: ${JSON.stringify(expectedSeverities)}`,
    timeout: 10_000, intervals: [200, 400, 800]}).toEqual(expectedSeverities);

    const caucasian = await hierNode(page, ['F', 'Caucasian'], 'read');
    expect(caucasian.childCaptions,
      'the SEVERITY grandchild this step unchecks is not among the children of F / Caucasian')
      .toContain(grandchild);
    expect(caucasian.childHeaders,
      `the SEVERITY level must carry exactly one caption-less row, the "${levels[2]}" column header; ` +
      `dropped rows were ${JSON.stringify(caucasian.childHeaders)} and the value rows ` +
      `${JSON.stringify(caucasian.childCaptions)}`)
      .toEqual([levels[2]]);
    expect(caucasian.glyph,
      'F / Caucasian must read checked before the grandchild is unchecked — otherwise the ' +
      'indeterminate read below is not a transition and proves nothing')
      .toBe(GLYPH.checked);

    const postUncheck = expected.femaleRows - expected.femaleCaucasianGrandchild;
    await hierNode(page, ['F', 'Caucasian', grandchild], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'unchecking the SEVERITY grandchild did not drop its rows from the filter: expected ' +
        `${expected.femaleRows} female rows minus the ${expected.femaleCaucasianGrandchild} that are ` +
        `Caucasian / ${grandchild} = ${postUncheck}`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(postUncheck);
    const after = await trueCountOf(page);

    const afterUncheck = await hierNode(page, ['F', 'Caucasian'], 'read');
    expect(afterUncheck.childCaptions.slice().sort(),
      'unchecking one SEVERITY value must not change which values the level lists: expected ' +
      `${JSON.stringify(expectedSeverities)}, now ${JSON.stringify(afterUncheck.childCaptions)}`)
      .toEqual(expectedSeverities);
    const siblingGlyphs: Record<string, string> = {};
    for (const sibling of expectedSeverities.filter((c) => c !== grandchild))
      siblingGlyphs[sibling] = glyphName((await hierNode(page, ['F', 'Caucasian', sibling], 'read')).glyph);
    expect(Object.keys(siblingGlyphs).length,
      `no SEVERITY sibling of "${grandchild}" was read, so the locality claim asserts nothing ` +
      `(derived values ${JSON.stringify(expectedSeverities)})`)
      .toBe(expectedSeverities.length - 1);
    for (const [sibling, glyph] of Object.entries(siblingGlyphs)) {
      expect(glyph,
        `SEVERITY sibling "${sibling}" must stay checked (${glyphName(GLYPH.checked)}) — only ` +
        `"${grandchild}" was unchecked; glyphs read: ${JSON.stringify(siblingGlyphs)}`)
        .toBe(glyphName(GLYPH.checked));
    }
    expect((await hierNode(page, ['F', 'Caucasian', grandchild], 'read')).glyph,
      `the unchecked grandchild "${grandchild}" must read unchecked (U+F0C8)`).toBe(GLYPH.unchecked);

    expect((await hierNode(page, ['F', 'Caucasian'], 'read')).glyph,
      'the direct RACE parent must move from checked to indeterminate (U+F146) when one of its ' +
      'SEVERITY children is unchecked')
      .toBe(GLYPH.indeterminate);
    expect((await hierNode(page, ['F'], 'read')).glyph,
      'the SEX grandparent must ALSO read indeterminate — the partial state has to propagate up two levels')
      .toBe(GLYPH.indeterminate);
    expect(after, 'unchecking a grandchild must drop rows').toBeLessThan(allF);
    expect(after, 'unchecking one SEVERITY value must not empty the branch').toBeGreaterThan(0);
  });

  await softStep('Step 9: GROK-16528 — reorder columns to RACE / SEX', async () => {
    await applyHierarchyState(page, {colNames: ['RACE', 'SEX'], allEnabled: true});
    await expect.poll(() => hierCaption(page), {timeout: 10_000, intervals: [100, 200, 400]}).toBe('RACE / SEX');
    const state = await page.evaluate(() => {
      const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
        .find(c => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
      const caption = card!.querySelector('.d4-filter-column-name')!.textContent!.trim();
      const topCaptions: string[] = [];
      for (const n of card!.querySelectorAll('.d4-tree-view-node')) {
        const val = n.querySelector('.d4-hierarchical-filter-caption-value');
        if (val) topCaptions.push(val.textContent!.trim());
      }
      return {caption, topCaptions};
    });
    expect(state.caption).toBe('RACE / SEX');
    expect(state.topCaptions).toContain('Caucasian');
    expect(state.topCaptions).toContain('Asian');
    expect(state.topCaptions).not.toContain('F');
  });

  await softStep('Step 12: Open a fresh demog view and add a second boolean column (SEX_bool)', async () => {
    const result = await openDemogWithSexBool(page, datasetPath);
    expect(result.type).toBe('bool');
    expect(result.hasControl).toBe(true);
  });

  await softStep('Step 13 / Step 14: Open Filter Panel — Combined Boolean auto-added, nothing toggled', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-bool-combined-filter').waitFor({timeout: 10000});
    // The header indicator renders a pass after the card does, and the step reads it — waiting on
    // the card alone lands here with the indicator still blank.
    await page.evaluate(() => (window as any).__poll(() => {
      const e = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator');
      return e ? e.textContent!.trim() : '';
    }, (t: string) => t.length > 0, 2000, 25));
    const state = await page.evaluate(() => {
      const boolCards = document.querySelectorAll('.d4-bool-combined-filter').length;
      const ind = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator');
      return {
        boolCards,
        indicator: ind ? ind.textContent!.trim() : null,
        trueCount: grok.shell.tv.dataFrame.filter.trueCount,
      };
    });
    expect(state.boolCards).toBe(1);
    expect(state.indicator).toBe('0');
    expect(state.trueCount).toBe(ROW_COUNT);
  });

  let trueCountOR2 = 0;
  await softStep('Step 15: Toggle first boolean column (OR) — filter narrows, counter reads 1', async () => {
    const expected = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const boolNames = df.columns.toList().filter((c: any) => c.type === 'bool').map((c: any) => c.name);
      const values = df.col(boolNames[0]).toList();
      return {boolNames, trueRows: values.filter((v: any) => v === true).length};
    });
    expect(expected.boolNames.length).toBe(2);
    expect(expected.trueRows).toBeGreaterThan(0);
    expect(expected.trueRows).toBeLessThan(rowCount);

    await page.evaluate(async () => {
      // The step reads the header indicator as well as the count, and the panel repaints it a pass
      // later — stamping the count alone returns while the indicator still reads its old value.
      const stamp = () => {
        const ind = document.querySelector(
          '[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator');
        return `${grok.shell.tv.dataFrame.filter.trueCount}|${ind ? ind.textContent!.trim() : ''}`;
      };
      const was = stamp();
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'bool-columns') {
          window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, false], 'false': [false, false], mode: 'OR'});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
      await (window as any).__moved(stamp, was, 600);
    });
    const state = await page.evaluate(() => {
      const ind = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator');
      return {trueCount: grok.shell.tv.dataFrame.filter.trueCount, indicator: ind ? ind.textContent!.trim() : null};
    });
    expect(state.trueCount).toBe(expected.trueRows);
    expect(state.indicator).toBe('1');
  });

  await softStep('Step 16: Add second column, then switch OR → AND — AND never widens', async () => {
    const expected = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const boolNames = df.columns.toList().filter((c: any) => c.type === 'bool').map((c: any) => c.name);
      const a = df.col(boolNames[0]).toList();
      const b = df.col(boolNames[1]).toList();
      let or = 0;
      let and = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (a[i] === true || b[i] === true) or++;
        if (a[i] === true && b[i] === true) and++;
      }
      return {or, and};
    });
    expect(expected.and).toBeGreaterThan(0);
    expect(expected.and).toBeLessThan(expected.or);
    expect(expected.or).toBeLessThan(rowCount);

    trueCountOR2 = await applyBoolState(page, {'true': [true, true], 'false': [false, false], mode: 'OR'});
    expect(trueCountOR2).toBe(expected.or);
    const trueCountAND = await applyBoolState(page, {'true': [true, true], 'false': [false, false], mode: 'AND'});
    expect(trueCountAND).toBe(expected.and);
    expect(trueCountAND).toBeLessThan(trueCountOR2);
    expect(trueCountAND).toBeGreaterThan(0);
    expect(trueCountAND).toBeLessThan(rowCount);
  });

  await softStep('Step 19: Remove All then reopen — the combined boolean card is recreated on its own', async () => {
    await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
    await expect.poll(async () => page.locator('[name="viewer-Filters"] .d4-filter').count(),
      {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
    const emptied = await page.evaluate(() => ({
      cards: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
      boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
      trueCount: grok.shell.tv.dataFrame.filter.trueCount,
    }));
    expect(emptied.cards, 'Remove All left cards in the panel').toBe(0);
    expect(emptied.boolCards, 'the combined boolean card survived Remove All').toBe(0);
    expect(emptied.trueCount, 'Remove All must release all filtering').toBe(rowCount);

    await page.evaluate(closeFilterPanelInPage);
    await expect.poll(async () => page.locator('[name="viewer-Filters"]').count(),
      {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);

    await page.locator('.d4-ribbon-panel [name="icon-filter"]').first().click();
    await page.locator('.d4-bool-combined-filter').first().waitFor({timeout: 20_000});
    const reopened = await page.evaluate(() => ({
      boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
      boolColumns: grok.shell.tv.dataFrame.columns.toList().filter((c: any) => c.type === 'bool').length,
      trueCount: grok.shell.tv.dataFrame.filter.trueCount,
    }));
    expect(reopened.boolCards, 'reopening must recreate exactly one combined boolean card').toBe(1);
    expect(reopened.boolColumns, 'the table must still carry more than one boolean column').toBeGreaterThan(1);
    expect(reopened.trueCount, 'a freshly recreated card must not be filtering anything').toBe(rowCount);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
