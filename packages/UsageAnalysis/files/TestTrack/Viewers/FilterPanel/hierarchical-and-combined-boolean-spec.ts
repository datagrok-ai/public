/* ---
realizes: [filters.cp.hierarchical-and-combined-boolean]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;
declare const window: any;

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

const GLYPH = {unchecked: 0xF0C8, checked: 0xF14A, indeterminate: 0xF146};

interface HierNodeState {
  found: boolean;
  glyph: number;
  expanded: boolean;
  childCaptions: string[];
  childHeaders: string[];
}

const glyphName = (g: number): string => 'U+' + g.toString(16).toUpperCase().padStart(4, '0');

async function hierNode(page: any, path: string[],
  action: 'probe' | 'read' | 'expand' | 'toggle'): Promise<HierNodeState> {
  return await page.evaluate(({path, action}: {path: string[], action: string}) => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
    if (!card)
      throw new Error('no hierarchical filter card is painted in the Filter Panel — its caption is the '
        + 'only one carrying "/", and no card in the panel carries one');
    const hostOf = (n: HTMLElement) =>
      n.parentElement?.querySelector(':scope > .d4-tree-view-group-host') as HTMLElement | null;
    const directChildren = (host: HTMLElement) => [...host.children].flatMap((c) =>
      c.matches('.d4-tree-view-node') ? [c as HTMLElement]
        : [...c.querySelectorAll(':scope > .d4-tree-view-node')] as HTMLElement[]);
    const captionOf = (n: HTMLElement) =>
      (n.querySelector('.d4-hierarchical-filter-caption-value')?.textContent ?? '').trim();
    const headerLabelOf = (n: HTMLElement) =>
      (n.querySelector('.d4-tree-view-item-label')?.textContent ?? n.textContent ?? '').trim();

    let node: HTMLElement | null = null;
    for (let level = 0; level < path.length; level++) {
      const value = path[level];
      let candidates: HTMLElement[];
      if (level === 0) {
        const rootHost = card.querySelector('.d4-tree-view-root > .d4-tree-view-group-host') as HTMLElement | null;
        if (!rootHost)
          throw new Error('the hierarchical filter card paints no tree root '
            + '(.d4-tree-view-root > .d4-tree-view-group-host), so its level-0 nodes cannot be addressed');
        candidates = directChildren(rootHost).filter((c) => captionOf(c) === value);
      }
      else {
        const host = hostOf(node!);
        candidates = host ? directChildren(host).filter((c) => captionOf(c) === value) : [];
      }
      if (candidates.length === 0) {
        if (action === 'probe')
          return {found: false, glyph: 0, expanded: false, childCaptions: [], childHeaders: []};
        throw new Error(`hierarchical tree node "${value}" not found (path ${path.join(' / ')}) — `
          + `it was looked for ${level === 0 ? 'inside the hierarchical card' : 'among the direct children of "' + path[level - 1] + '"'}`);
      }
      if (candidates.length > 1)
        throw new Error(`hierarchical tree node "${value}" is AMBIGUOUS (path ${path.join(' / ')}): `
          + `${candidates.length} nodes ${level === 0 ? 'inside the hierarchical card' : 'among the direct children of "' + path[level - 1] + '"'} `
          + 'carry that caption, so addressing it by caption would silently pick one of them');
      node = candidates[0];
    }
    const tri = node!.querySelector(':scope > .d4-tree-view-tri') as HTMLElement | null;
    if (action === 'expand') {
      if (!tri) throw new Error(`hierarchical tree node "${path.join(' / ')}" carries no expander`);
      if (!tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
    }
    if (action === 'toggle') {
      const cb = node!.querySelector('input.d4-hierarchical-filter-checkbox') as HTMLElement | null;
      if (!cb) throw new Error(`hierarchical tree node "${path.join(' / ')}" carries no checkbox`);
      cb.click();
    }
    const sub = node!.querySelector('.d4-hierarchical-filter-checkbox-substitute');
    const host = hostOf(node!);
    const children = host ? directChildren(host) : [];
    return {
      found: true,
      glyph: (sub?.textContent ?? '').codePointAt(0) ?? 0,
      expanded: !!tri?.classList.contains('d4-tree-view-tri-expanded'),
      childCaptions: children.map(captionOf).filter((c) => c !== ''),
      childHeaders: children.filter((c) => captionOf(c) === '').map(headerLabelOf),
    };
  }, {path, action});
}

const trueCountOf = async (page: any): Promise<number> =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

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

  const isAmbientNoise = (text: string) =>
    /ProjectMeta\.publish/.test(text) ||
    /project_meta\.dart/.test(text) ||
    /could not be cloned/i.test(text) ||
    /Failed to load resource/.test(text) ||
    /favicon/.test(text);
  const pageErrors: string[] = [];
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isAmbientNoise(msg.text())) pageErrors.push(msg.text());
  });
  page.on('pageerror', (err) => {
    if (!isAmbientNoise(err.message)) pageErrors.push(`pageerror: ${err.message}`);
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBe(5850);

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
    await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Hierarchical');
    await expect.poll(async () => page.evaluate(() =>
      grok.shell.tv.getFiltersGroup().filters.filter((f: any) => f.filterType === 'hierarchical').length),
    {timeout: 15_000, intervals: [400, 800, 1500]}).toBe(1);
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          window.grok_GridFilterBase_ApplyState(f.dart || f, {
            type: 'hierarchical', active: true, colNames: ['SEX', 'RACE'], allEnabled: true,
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await expect.poll(async () => (await hierNode(page, ['F'], 'probe')).found,
      {message: 'the SEX root node "F" never rendered — the hierarchy was not applied to the card',
        timeout: 10_000, intervals: [300, 600, 1200]}).toBe(true);
    const caption = await page.evaluate(() => {
      for (const c of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
        const cn = c.querySelector('.d4-filter-column-name');
        if (cn && cn.textContent!.includes('/')) return cn.textContent!.trim();
      }
      return '';
    });
    expect(caption).toBe('SEX / RACE');
    const trueCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(trueCount).toBe(5850);
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
    await page.evaluate(async (colNames: string[]) => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          window.grok_GridFilterBase_ApplyState(f.dart || f, {colNames, allEnabled: true});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    }, levels);
    const readCaption = async () => page.evaluate(() => {
      for (const c of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
        const cn = c.querySelector('.d4-filter-column-name');
        if (cn && cn.textContent!.includes('/')) return cn.textContent!.trim();
      }
      return '';
    });
    await expect.poll(readCaption,
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

  let hlTrueCount = 0;
  await softStep('Step 6-restore: re-establish SEX / RACE + Caucasian criterion for the layout save', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          window.grok_GridFilterBase_ApplyState(f.dart || f, {colNames: ['SEX', 'RACE'], allEnabled: true});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await expect.poll(async () => (await hierNode(page, ['F'], 'probe')).found,
      {message: 'the SEX roots never rebuilt after the hierarchy was reset to SEX / RACE',
        timeout: 15_000, intervals: [300, 600, 1200]}).toBe(true);
    await hierNode(page, ['F'], 'expand');
    await expect.poll(async () => (await hierNode(page, ['F'], 'read')).childCaptions,
      {message: 'the RACE children of F never rendered after the expander click',
        timeout: 10_000, intervals: [200, 400, 800]}).toContain('Caucasian');
    await hierNode(page, ['F', 'Caucasian'], 'toggle');
    await expect.poll(async () => trueCountOf(page),
      {message: 'the criterion about to be saved is not the Caucasian-female one derived from the raw '
        + `SEX / RACE columns (${derived.caucasianFemale})`,
      timeout: 10_000, intervals: [200, 400, 800]}).toBe(derived.caucasianFemale);
    hlTrueCount = await trueCountOf(page);
  });

  let hlLayoutId = '';
  await softStep('Step 10: Save the layout with the hierarchical criterion active', async () => {
    hlLayoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    expect(hlLayoutId).toBeTruthy();
    await expect.poll(async () => page.evaluate(async (id: string) => {
      try { return !!(await grok.dapi.layouts.find(id)); }
      catch (_) { return false; }
    }, hlLayoutId),
    {message: 'the saved layout cannot be fetched back from the server — saveLayout() stamps the id '
      + 'client-side before the round-trip, so a dapi.layouts.save that silently failed leaves the id '
      + 'just as non-empty and Step 11 would re-apply nothing',
    timeout: 20_000, intervals: [500, 1000, 2000]}).toBe(true);
  });

  await softStep('Step 9: GROK-16528 — reorder columns to RACE / SEX', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          window.grok_GridFilterBase_ApplyState(f.dart || f, {colNames: ['RACE', 'SEX'], allEnabled: true});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await page.waitForTimeout(1200);
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

  try {
    await softStep('Step 11: Layout round-trip (hierarchical) — close panel, re-apply saved layout', async () => {
      const beforeClose = await page.evaluate(async () => {
        const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
          .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.includes('/'));
        if (!card)
          throw new Error('no hierarchical filter card is painted in the Filter Panel — its caption is the only '
            + 'one carrying "/", and no card in the panel carries one');
        let clicked = false;
        for (const n of card.querySelectorAll('.d4-tree-view-node')) {
          const val = n.querySelector('.d4-hierarchical-filter-caption-value');
          if (val?.textContent?.trim() === 'Caucasian') {
            const cb = n.querySelector('input.d4-hierarchical-filter-checkbox') as HTMLElement | null;
            if (!cb)
              throw new Error('the RACE root node "Caucasian" carries no checkbox, so it could not be ticked');
            cb.click();
            clicked = true;
            break;
          }
        }
        if (!clicked) throw new Error('RACE root node "Caucasian" not found in the hierarchical card');
        await new Promise(r => setTimeout(r, 900));
        return grok.shell.tv.dataFrame.filter.trueCount;
      });
      expect(beforeClose).toBeGreaterThan(0);
      expect(beforeClose).toBeLessThan(rowCount);

      const result = await page.evaluate(async (id: string) => {
        const fv = document.querySelector('[name="viewer-Filters"]');
        let el: any = fv;
        while (el && !el.classList.contains('panel-base')) el = el.parentElement;
        const closeBtn = el?.querySelector('[name="Close"]') as HTMLElement | null;
        if (!closeBtn) throw new Error('Filter Panel titlebar [name="Close"] not found');
        closeBtn.click();
        await new Promise(r => setTimeout(r, 1000));
        const afterClose = grok.shell.tv.dataFrame.filter.trueCount;
        const panelsAfterClose = document.querySelectorAll('[name="viewer-Filters"]').length;

        const saved = await grok.dapi.layouts.find(id);
        grok.shell.tv.loadLayout(saved);
        grok.shell.tv.getFiltersGroup();
        await new Promise(r => setTimeout(r, 3500));
        const afterRestore = grok.shell.tv.dataFrame.filter.trueCount;

        let caption = '';
        for (const c of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
          const cn = c.querySelector('.d4-filter-column-name');
          if (cn && cn.textContent!.includes('/')) { caption = cn.textContent!.trim(); break; }
        }
        return {afterClose, panelsAfterClose, afterRestore, caption};
      }, hlLayoutId);
      expect(result.afterClose).toBe(rowCount);
      expect(result.afterClose).toBeGreaterThan(beforeClose);
      expect(result.panelsAfterClose).toBe(0);
      expect(result.afterRestore).toBe(hlTrueCount);
      expect(result.caption).toBe('SEX / RACE');
    });
  } finally {
    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, hlLayoutId).catch(() => {});
  }

  await softStep('Step 12: Open a fresh demog view and add a second boolean column (SEX_bool)', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1500));
      await df.columns.addNewCalculated('SEX_bool', '${SEX} == "F"');
      const col = df.col('SEX_bool');
      return {type: col.type, hasControl: !!df.col('CONTROL')};
    });
    expect(result.type).toBe('bool');
    expect(result.hasControl).toBe(true);
  });

  await softStep('Step 13 / Step 14: Open Filter Panel — Combined Boolean auto-added, nothing toggled', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.locator('.d4-bool-combined-filter').waitFor({timeout: 10000});
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
    expect(state.trueCount).toBe(5850);
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

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'bool-columns') {
          window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, false], 'false': [false, false], mode: 'OR'});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await page.waitForTimeout(600);
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

    trueCountOR2 = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'bool-columns') {
          window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, true], 'false': [false, false], mode: 'OR'});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(trueCountOR2).toBe(expected.or);
    await page.waitForTimeout(500);
    const trueCountAND = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'bool-columns') {
          window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, true], 'false': [false, false], mode: 'AND'});
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(trueCountAND).toBe(expected.and);
    expect(trueCountAND).toBeLessThan(trueCountOR2);
    expect(trueCountAND).toBeGreaterThan(0);
    expect(trueCountAND).toBeLessThan(rowCount);
  });

  const probeProject = 'zz-grok16488-' + Date.now();
  let boolProjectId = '';
  try {
    await softStep('Step 17: GROK-16488 — save project, reopen, remove combined boolean card, count returns to full, no console error', async () => {
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        for (const f of fg.filters) {
          if (f.filterType === 'bool-columns') {
            window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, true], 'false': [true, true], mode: 'OR'});
            grok.shell.tv.dataFrame.rows.requestFilter();
            break;
          }
        }
      });
      await page.waitForTimeout(500);
      expect(await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
        'the combined boolean reset left the table filtered — the project would be saved narrowed')
        .toBe(5850);

      const saved = await saveProjectViaUI(page, probeProject);
      boolProjectId = saved.projectId;

      await page.evaluate(async (id: string) => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 1200));
        const proj = await grok.dapi.projects.find(id);
        await proj.open();
        await new Promise(r => setTimeout(r, 4000));
        grok.shell.tv.getFiltersGroup();
        await new Promise(r => setTimeout(r, 2500));
      }, boolProjectId);

      await page.locator('.d4-bool-combined-filter').waitFor({timeout: 15000});
      const afterReopen = await page.evaluate(() => ({
        boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filterPanel: document.querySelectorAll('[name="viewer-Filters"]').length,
      }));
      expect(afterReopen.filterPanel).toBeGreaterThan(0);
      expect(afterReopen.boolCards).toBe(1);

      const preRemoval = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const boolNames = df.columns.toList().filter((c: any) => c.type === 'bool').map((c: any) => c.name);
        const expectedRows = df.col(boolNames[0]).toList().filter((v: any) => v === true).length;
        const fg = grok.shell.tv.getFiltersGroup();
        for (const f of fg.filters) {
          if (f.filterType === 'bool-columns') {
            window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, false], 'false': [false, false], mode: 'OR'});
            grok.shell.tv.dataFrame.rows.requestFilter();
            break;
          }
        }
        return {expectedRows, trueCount: grok.shell.tv.dataFrame.filter.trueCount};
      });
      expect(preRemoval.expectedRows).toBeGreaterThan(0);
      expect(preRemoval.expectedRows).toBeLessThan(5850);
      expect(preRemoval.trueCount).toBe(preRemoval.expectedRows);

      const errorsBefore = pageErrors.length;

      await page.evaluate(() => {
        const boolCard = document.querySelector('.d4-bool-combined-filter');
        let host: any = boolCard;
        while (host && !host.classList.contains('d4-filter')) host = host.parentElement;
        const x = host?.querySelector('[name="icon-times"]') as HTMLElement | null;
        if (!x) throw new Error('the combined boolean card carries no [name="icon-times"] remove icon — '
          + 'the removal this step is about never happened');
        x.click();
      });

      await expect.poll(async () => page.evaluate(() => ({
        boolCards: document.querySelectorAll('.d4-bool-combined-filter').length,
        trueCount: grok.shell.tv.dataFrame.filter.trueCount,
      })),
      {message: 'removing the combined boolean card did not both take the card away and release its '
        + `narrowing back to the full ${rowCount} rows`,
      timeout: 15_000, intervals: [200, 400, 800]}).toEqual({boolCards: 0, trueCount: rowCount});

      const removalSamples: string[] = [];
      for (let i = 0; i < 6; i++) {
        await page.waitForTimeout(400);
        removalSamples.push(JSON.stringify(pageErrors.slice(errorsBefore)));
      }
      expect(Array.from(new Set(removalSamples)),
        'GROK-16488 — removing the combined boolean card must not throw at any point of a 2.4s window '
        + `after the click, so an async throw arriving late is caught too; samples: ${removalSamples.join(' ; ')}`)
        .toEqual(['[]']);
    });
  } finally {
    await deleteProjectWithCleanup(page, {projectId: boolProjectId});
  }

  let boolLayoutId = '';
  try {
    await softStep('Step 18: Layout round-trip (combined boolean) — save active state, close, re-apply', async () => {
      const removeClicks = await page.evaluate(() => {
        let clicks = 0;
        for (const boolCard of document.querySelectorAll('.d4-bool-combined-filter')) {
          let host: any = boolCard;
          while (host && !host.classList.contains('d4-filter')) host = host.parentElement;
          const x = host?.querySelector('[name="icon-times"]') as HTMLElement | null;
          if (!x) throw new Error('a .d4-bool-combined-filter card carries no [name="icon-times"] remove icon');
          x.click();
          clicks++;
        }
        return clicks;
      });
      await expect.poll(async () => page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      })),
      {message: `the combined boolean card did not go away after ${removeClicks} remove-icon click(s) — ` +
        'the menu drive below would then be satisfied by a card it did not create',
      timeout: 15_000, intervals: [400, 800, 1500]}).toEqual({cards: 0, filters: 0});

      const beforeAdd = await page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      }));
      expect(beforeAdd,
        'a combined boolean filter survived the removal, so the Add Filter > Combined Boolean drive ' +
        `cannot be shown to have created anything; read ${JSON.stringify(beforeAdd)}`)
        .toEqual({cards: 0, filters: 0});

      await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Combined Boolean');
      await expect.poll(async () => page.evaluate(() => ({
        cards: document.querySelectorAll('.d4-bool-combined-filter').length,
        filters: grok.shell.tv.getFiltersGroup().filters
          .filter((f: any) => f.filterType === 'bool-columns').length,
      })),
      {message: 'Add Filter > Combined Boolean did not take the panel from 0 combined boolean cards to ' +
        'exactly 1 — either the leaf created nothing, or it fired twice and left two identical cards',
      timeout: 15_000, intervals: [400, 800, 1500]}).toEqual({cards: 1, filters: 1});
      const result = await page.evaluate(async () => {
        let stage = 'toggling the combined boolean card before the save';
        try {
        const fg = grok.shell.tv.getFiltersGroup();
        for (const f of fg.filters) {
          if (f.filterType === 'bool-columns') {
            window.grok_GridFilterBase_ApplyState(f.dart ?? f, {'true': [true, false], 'false': [false, false], mode: 'OR'});
            grok.shell.tv.dataFrame.rows.requestFilter();
            break;
          }
        }
        await new Promise(r => setTimeout(r, 700));
        const savedCount = grok.shell.tv.dataFrame.filter.trueCount;

        stage = 'saving the layout';
        const layout = grok.shell.tv.saveLayout();
        await grok.dapi.layouts.save(layout);
        const id = layout.id;
        stage = 'fetching the saved layout back from the server';
        let serverFound = false;
        for (let i = 0; i < 20 && !serverFound; i++) {
          try { serverFound = !!(await grok.dapi.layouts.find(id)); }
          catch (_) { serverFound = false; }
          if (!serverFound) await new Promise(r => setTimeout(r, 500));
        }

        stage = 'closing the Filter Panel';
        const fv = document.querySelector('[name="viewer-Filters"]');
        let el: any = fv;
        while (el && !el.classList.contains('panel-base')) el = el.parentElement;
        const closeBtn = el?.querySelector('[name="Close"]') as HTMLElement | null;
        if (!closeBtn) throw new Error('the Filter Panel titlebar carries no [name="Close"] control — '
          + 'the close gesture the restore below is measured against never happened');
        closeBtn.click();
        await new Promise(r => setTimeout(r, 1000));
        const afterClose = grok.shell.tv.dataFrame.filter.trueCount;
        const cardsAfterClose = document.querySelectorAll('.d4-bool-combined-filter').length;

        stage = 're-applying the saved layout';
        const saved = await grok.dapi.layouts.find(id);
        grok.shell.tv.loadLayout(saved);
        grok.shell.tv.getFiltersGroup();
        await new Promise(r => setTimeout(r, 3500));
        const afterRestore = grok.shell.tv.dataFrame.filter.trueCount;
        const hasBoolCard = document.querySelectorAll('.d4-bool-combined-filter').length;

        return {id, serverFound, savedCount, afterClose, cardsAfterClose, afterRestore, hasBoolCard};
        }
        catch (e: any) {
          throw new Error(`the combined boolean layout round-trip failed while ${stage}: ${e?.message ?? e}`);
        }
      });
      boolLayoutId = result.id;
      expect(result.serverFound,
        'the saved layout could not be fetched back from the server — saveLayout() stamps the id '
        + 'client-side before the round-trip, so a dapi.layouts.save that silently failed would leave '
        + 'the re-apply below round-tripping nothing')
        .toBe(true);
      expect(result.savedCount).toBeGreaterThan(0);
      expect(result.savedCount).toBeLessThan(rowCount);
      expect(result.cardsAfterClose).toBe(0);
      expect(result.afterClose).toBe(rowCount);
      expect(result.hasBoolCard).toBeGreaterThan(0);
      expect(result.afterRestore).toBe(result.savedCount);
    });
  } finally {
    await page.evaluate(async (id: string) => {
      if (!id) return;
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, boolLayoutId).catch(() => {});
  }

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

    await page.evaluate(() => {
      const fv = document.querySelector('[name="viewer-Filters"]');
      let el: any = fv;
      while (el && !el.classList.contains('panel-base')) el = el.parentElement;
      const closeBtn = el?.querySelector('[name="Close"]') as HTMLElement | null;
      if (!closeBtn) throw new Error('the Filter Panel titlebar carries no [name="Close"] control — the '
        + 'close gesture the reopen below is measured against never happened');
      closeBtn.click();
    });
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
