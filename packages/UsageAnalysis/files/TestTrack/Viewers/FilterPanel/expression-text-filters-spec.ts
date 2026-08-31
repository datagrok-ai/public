/* ---
realizes: [filters.cp.expression-and-text-ui]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addCardViaColumnSelector, cardCount, clickResetCriteriaIcon, trueCount} from '../../helpers/filter-panel';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const beerPath = 'System:DemoFiles/beer.csv';
const fullCount_demog = 5850;
const fullCount_beer = 118;

async function expressionRuleCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const st = grok.shell.tv.getFiltersGroup().getStates(null, 'expression');
    return st?.[0]?.gridNames?.length ?? 0;
  });
}

async function filterSummaryText(page: Page): Promise<string> {
  return page.evaluate(() => {
    const el = grok.shell.tv.getFiltersGroup().getFilterSummary();
    return (el?.textContent ?? '').trim();
  });
}

async function removeAllCards(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(async () => cardCount(page), {timeout: 20_000, intervals: [400, 800, 1500]}).toBe(0);
}

async function resetToExpressionCard(page: Page): Promise<void> {
  await removeAllCards(page);
  expect(await trueCount(page), 'removing every card leaves demog unfiltered before the next rule is built')
    .toBe(fullCount_demog);
  await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Expression');
  await page.locator('.d4-expression-filter').first().waitFor({state: 'attached', timeout: 15000});
}

async function commitFormRule(page: Page, column: string, op: string, value: string): Promise<void> {
  await page.evaluate(async ({col, operation, val}) => {
    const card = document.querySelector('.d4-expression-filter') as HTMLElement;
    const selects = card.querySelectorAll('.ui-input-choice select');
    const colSel = selects[0] as HTMLSelectElement;
    const opSel = selects[1] as HTMLSelectElement;
    const setSel = Object.getOwnPropertyDescriptor(window.HTMLSelectElement.prototype, 'value')!.set!;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setSel.call(colSel, col);
    colSel.dispatchEvent(new Event('input', {bubbles: true}));
    colSel.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 700));
    setSel.call(opSel, operation);
    opSel.dispatchEvent(new Event('input', {bubbles: true}));
    opSel.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 400));
    const valInput = card.querySelector('.ui-input-text input') as HTMLInputElement;
    setInp.call(valInput, val);
    valInput.dispatchEvent(new Event('input', {bubbles: true}));
    valInput.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 300));
    (card.querySelector('.fal.fa-plus') as HTMLElement).click();
    await new Promise((r) => setTimeout(r, 1300));
  }, {col: column, operation: op, val: value});
}

async function andOrControl(page: Page, column: string | null, toggle: boolean): Promise<string> {
  return page.evaluate(async ({col, doToggle}) => {
    const name = col === null ? 'expression' : col;
    const card = col === null
      ? (document.querySelector('.d4-expression-filter') as HTMLElement | null)?.closest('.d4-filter')
      : [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
        .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent ?? '').trim() === col);
    if (!card) throw new Error(`the panel carries no ${name} card to read its AND/OR control on`);
    const header = card.querySelector('.d4-filter-header');
    if (!header) throw new Error(`the ${name} card carries no .d4-filter-header`);
    const read = () => [...header.querySelectorAll('div')]
      .find((d) => { const t = (d.textContent ?? '').trim(); return t === 'OR' || t === 'AND'; }) as HTMLElement | undefined;
    const before = read();
    if (!before) throw new Error(`the ${name} card header carries no AND/OR control to read or click`);
    if (!doToggle) return (before.textContent ?? '').trim();
    before.click();
    await new Promise((r) => setTimeout(r, 1500));
    return (read()?.textContent ?? '').trim();
  }, {col: column, doToggle: toggle});
}

async function removeQueryRow(page: Page, rowIndex: number, ruleCount: number): Promise<boolean> {
  return page.evaluate(async ({idx, count}) => {
    const card = (document.querySelector('.d4-expression-filter') as HTMLElement).closest('.d4-filter')!;
    const overlay = card.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement;
    const r = overlay.getBoundingClientRect();
    const rowH = r.height / Math.max(count, 1);
    const x = Math.round(r.left + r.width / 2);
    const y = Math.round(r.top + idx * rowH + rowH / 2);
    overlay.dispatchEvent(new PointerEvent('pointermove', {bubbles: true, clientX: x, clientY: y}));
    overlay.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: x, clientY: y}));
    await new Promise((res) => setTimeout(res, 400));
    overlay.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2}));
    const items = () => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])];
    };
    const deadline = Date.now() + 5000;
    let leaf: Element | null = null;
    while (!leaf && Date.now() < deadline) {
      await new Promise((res) => setTimeout(res, 100));
      leaf = items().find((it) => (it.textContent ?? '').trim() === 'Remove Query')?.closest('.d4-menu-item') ?? null;
    }
    if (!leaf) {
      const visible = items().map((it) => (it.textContent ?? '').trim()).join(' | ');
      document.body.click();
      throw new Error(`the rule-row context menu offered no "Remove Query" leaf; visible: ${visible}`);
    }
    (leaf as HTMLElement).click();
    await new Promise((res) => setTimeout(res, 1300));
    return true;
  }, {idx: rowIndex, count: ruleCount});
}

async function toggleRuleCheckbox(page: Page, rowIndex: number, ruleCount: number): Promise<boolean[]> {
  return page.evaluate(async ({idx, count}) => {
    const card = (document.querySelector('.d4-expression-filter') as HTMLElement).closest('.d4-filter')!;
    const overlay = card.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement;
    const r = overlay.getBoundingClientRect();
    const rowH = r.height / Math.max(count, 1);
    const x = Math.round(r.left + 8);
    const y = Math.round(r.top + idx * rowH + rowH / 2);
    for (const type of ['pointermove', 'mousemove'])
      overlay.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: x, clientY: y}));
    await new Promise((res) => setTimeout(res, 300));
    for (const type of ['mousedown', 'mouseup', 'click'])
      overlay.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0}));
    await new Promise((res) => setTimeout(res, 1200));
    const s = grok.shell.tv.getFiltersGroup().getStates(null, 'expression')[0];
    return (s?.gridValues ?? []) as boolean[];
  }, {idx: rowIndex, count: ruleCount});
}

async function toggleExpressionMode(page: Page): Promise<string> {
  return page.evaluate(async () => {
    const card = (document.querySelector('.d4-expression-filter') as HTMLElement).closest('.d4-filter')!;
    (card.querySelector('[name="icon-italic"]') as HTMLElement).click();
    await new Promise((r) => setTimeout(r, 1300));
    return grok.shell.tv.getFiltersGroup().getStates(null, 'expression')[0].expressionMode;
  });
}

async function typeFreeText(page: Page, expr: string): Promise<{value: string; invalid: boolean}> {
  return page.evaluate(async (e) => {
    const card = (document.querySelector('.d4-expression-filter') as HTMLElement).closest('.d4-filter')!;
    const inputs = [...card.querySelectorAll('input.d4-search-input[placeholder="Search"]')] as HTMLInputElement[];
    const input = inputs.find((i) => i.getBoundingClientRect().width > 0)!;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setInp.call(input, '');
    input.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 300));
    input.focus();
    setInp.call(input, e);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    input.dispatchEvent(new KeyboardEvent('keyup', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    await new Promise((r) => setTimeout(r, 1800));
    const live = [...card.querySelectorAll('input.d4-search-input[placeholder="Search"]')]
      .find((i) => (i as HTMLInputElement).getBoundingClientRect().width > 0) as HTMLInputElement | undefined;
    return {value: live?.value ?? '', invalid: live?.classList.contains('d4-invalid') ?? true};
  }, expr);
}

async function heightRuleCount(page: Page, threshold: number): Promise<{passing: number; missing: number}> {
  return page.evaluate((t) => {
    const df = grok.shell.tv.dataFrame;
    const c = df.col('HEIGHT');
    let passing = 0; let missing = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (c.isNone(i)) { missing++; continue; }
      const x = c.get(i);
      if (x === null || x === undefined || Number.isNaN(x)) { missing++; continue; }
      if (x < t) passing++;
    }
    return {passing, missing};
  }, threshold);
}

async function pasteRegexValueAndCommit(page: Page, column: string, pasted: string): Promise<{value: string; count: number; rules: string[]}> {
  return page.evaluate(async ({col, text}) => {
    const card = document.querySelector('.d4-expression-filter') as HTMLElement;
    const selects = card.querySelectorAll('.ui-input-choice select');
    const colSel = selects[0] as HTMLSelectElement;
    const opSel = selects[1] as HTMLSelectElement;
    const setSel = Object.getOwnPropertyDescriptor(window.HTMLSelectElement.prototype, 'value')!.set!;
    setSel.call(colSel, col);
    colSel.dispatchEvent(new Event('input', {bubbles: true}));
    colSel.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 700));
    setSel.call(opSel, 'regex');
    opSel.dispatchEvent(new Event('input', {bubbles: true}));
    opSel.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 400));
    const valInput = card.querySelector('.ui-input-text input') as HTMLInputElement;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    valInput.focus();
    setInp.call(valInput, '');
    valInput.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 250));
    const dt = new DataTransfer();
    dt.setData('text', text);
    valInput.dispatchEvent(new ClipboardEvent('paste', {bubbles: true, cancelable: true, clipboardData: dt}));
    await new Promise((r) => setTimeout(r, 600));
    const value = valInput.value;
    (card.querySelector('.fal.fa-plus') as HTMLElement).click();
    await new Promise((r) => setTimeout(r, 1300));
    const st = grok.shell.tv.getFiltersGroup().getStates(null, 'expression')[0];
    return {value, count: grok.shell.tv.dataFrame.filter.trueCount, rules: st.gridNames};
  }, {col: column, text: pasted});
}

async function resetToAromaTextFilter(page: Page): Promise<void> {
  await removeAllCards(page);
  await addCardViaColumnSelector(page, 'Aroma');
  await page.waitForFunction(() => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent ?? '').trim() === 'Aroma');
    if (!card) return false;
    if (!card.querySelector('.d4-text-filter')) return false;
    if (card.querySelector('.d4-update-shadow')) return false;
    return !!card.querySelector('input.d4-search-input') && !!card.querySelector('input[type="range"]');
  }, null, {timeout: 120_000, polling: 250});
}

async function addAromaTerm(page: Page, term: string): Promise<number> {
  return page.evaluate(async (t) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent || '').trim() === 'Aroma') as HTMLElement;
    const input = card.querySelector('.d4-search-input') as HTMLInputElement;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    input.focus();
    setInp.call(input, t);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    input.dispatchEvent(new KeyboardEvent('keyup', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    await new Promise((r) => setTimeout(r, 1700));
    return grok.shell.tv.dataFrame.filter.trueCount;
  }, term);
}

async function setAromaFuzziness(page: Page, value: number): Promise<number> {
  return page.evaluate(async (val) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent || '').trim() === 'Aroma') as HTMLElement;
    const range = card.querySelector('input[type="range"]') as HTMLInputElement;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    range.focus();
    setInp.call(range, String(val));
    range.dispatchEvent(new Event('input', {bubbles: true}));
    range.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 1800));
    return grok.shell.tv.dataFrame.filter.trueCount;
  }, value);
}

test('Filter Panel — Expression filter and Text filter driven through their own UI', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: demogPath, withFilterPanel: true});

  await softStep('Setup — record the full demog row count', async () => {
    expect(await trueCount(page), 'baseline is the full demog row set').toBe(fullCount_demog);
  });

  let ageGt50 = 0;
  let trueCount_rule1 = 0;
  let trueCount_and = 0;
  let trueCount_or = 0;

  await softStep('Step 2 and Step 3 — commit AGE > 50 via the form; count drops, one rule named', async () => {
    await resetToExpressionCard(page);
    ageGt50 = await page.evaluate(() => {
      const c = grok.shell.tv.dataFrame.col('AGE'); let n = 0;
      for (let i = 0; i < grok.shell.tv.dataFrame.rowCount; i++) if (c.get(i) > 50) n++;
      return n;
    });
    await commitFormRule(page, 'AGE', '>', '50');
    trueCount_rule1 = await trueCount(page);
    expect(trueCount_rule1, 'AGE > 50 gives the independent AGE>50 row count').toBe(ageGt50);
    expect(trueCount_rule1, 'the first rule drops below the full demog count').toBeLessThan(fullCount_demog);
    expect(await expressionRuleCount(page), 'exactly one rule after the first commit').toBe(1);
    expect(await filterSummaryText(page), 'summary names the AGE rule').toContain('AGE');
  });

  await softStep('Step 4 and Step 5 — commit HEIGHT < 160 as a second rule; two rules named', async () => {
    await commitFormRule(page, 'HEIGHT', '<', '160');
    expect(await expressionRuleCount(page), 'two rules after the second commit').toBe(2);
    expect(await filterSummaryText(page), 'summary now names the HEIGHT rule too').toContain('HEIGHT');
    const orUnion = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const age = df.col('AGE');
      const height = df.col('HEIGHT');
      let n = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const overFifty = !age.isNone(i) && age.get(i) > 50;
        const underOneSixty = !height.isNone(i) && height.get(i) < 160;
        if (overFifty || underOneSixty) n++;
      }
      return n;
    });
    expect(orUnion,
      `the second rule must bring in rows the AGE rule alone does not select — HEIGHT < 160 OR AGE > 50 covers ${orUnion} rows against the ${trueCount_rule1} the AGE rule gave on its own, and were those equal a second rule that committed and did nothing would satisfy the count below`)
      .toBeGreaterThan(trueCount_rule1);
    expect(orUnion, 'the OR union must still exclude some row, or any filter at all would satisfy the count below')
      .toBeLessThan(fullCount_demog);
    const twoRule = await trueCount(page);
    expect(twoRule,
      `the card is in OR mode here, so the two committed rules must select exactly the ${orUnion} rows that are AGE > 50 OR HEIGHT < 160, derived independently off the two columns — a rule that committed without biting leaves the count at ${trueCount_rule1} and fails here, and a card that intersected instead of uniting lands below ${orUnion}`)
      .toBe(orUnion);
  });

  await softStep('Step 6 — AND vs OR ordering: AND result is a subset of OR', async () => {
    expect(await andOrControl(page, null, false), 'the two-rule card starts in OR mode').toBe('OR');
    trueCount_or = await trueCount(page);
    const afterToggle = await andOrControl(page, null, true);
    expect(afterToggle, 'the header toggle switches OR to AND').toBe('AND');
    trueCount_and = await trueCount(page);
    expect(trueCount_and, 'AND is an equal-or-stricter subset of OR for the same two rules')
      .toBeLessThanOrEqual(trueCount_or);
    expect(trueCount_and, 'AND is strictly smaller than OR here — the toggle changed the row set')
      .toBeLessThan(trueCount_or);
  });

  await softStep('Step 7 — Remove Query on the first rule; count rises, one rule remains', async () => {
    const before = await expressionRuleCount(page);
    expect(before, 'two rules present before Remove Query').toBe(2);
    const removed = await removeQueryRow(page, 0, before);
    expect(removed, 'the Remove Query context-menu leaf was reachable and clicked').toBe(true);
    expect(await expressionRuleCount(page), 'exactly one rule remains after Remove Query').toBe(1);
    const after = await trueCount(page);
    expect(after, 'removing the first rule raises the row count above the AND value').toBeGreaterThan(trueCount_and);
    expect(await filterSummaryText(page), 'summary now names only the surviving HEIGHT rule').toContain('HEIGHT');
  });

  let singleRuleCount = 0;
  let height160 = 0;

  await softStep('Step 8 — free-text mode reproduces the single-rule count losslessly', async () => {
    singleRuleCount = await trueCount(page);
    const derived160 = await heightRuleCount(page, 160);
    const derived150 = await heightRuleCount(page, 150);
    height160 = derived160.passing;
    const height150 = derived150.passing;
    expect(height160,
      `the rule the card is left holding is HEIGHT < 160, so the count it produces must be the ${height160} rows derived independently off the HEIGHT column, counting only rows that carry a value (${derived160.missing} rows do not, and a missing numeric reads back as a sentinel that would compare below any threshold)`)
      .toBe(singleRuleCount);
    expect(height150, 'HEIGHT < 150 must select some rows, or the typed predicate below proves nothing')
      .toBeGreaterThan(0);
    expect(height150,
      `HEIGHT < 150 must select strictly fewer rows (${height150}) than HEIGHT < 160 (${height160}), or typing the different predicate below could not be told apart from typing nothing`)
      .toBeLessThan(height160);
    expect(await andOrControl(page, null, false),
      'the card is still in the AND logic Step 6 left it in, so a free-text predicate typed alongside the committed HEIGHT < 160 rule intersects with it')
      .toBe('AND');

    const mode = await toggleExpressionMode(page);
    expect(mode, 'the icon-italic toggle switches the card to free-text mode').toBe('free-text');

    const tighter = await typeFreeText(page, '${HEIGHT} < 150');
    expect(tighter.invalid, `the free-text card rejected ${'${HEIGHT}'} < 150 as unparseable`).toBe(false);
    expect(await trueCount(page),
      `the free-text expression ${'${HEIGHT}'} < 150, ANDed with the committed HEIGHT < 160 rule, must move the row set to the ${height150} rows derived off the HEIGHT column — a predicate DIFFERENT from the one the card was already holding, so an ignored input, a swallowed Enter or a silently discarded expression leaves it at ${height160} and fails here`)
      .toBe(height150);

    await clickResetCriteriaIcon(page);
    await expect.poll(async () => trueCount(page), {timeout: 20_000, intervals: [400, 800, 1500]})
      .toBe(fullCount_demog);
    expect(await page.evaluate(() =>
      grok.shell.tv.getFiltersGroup().getStates(null, 'expression')[0].expressionMode),
    'clearing the criteria leaves the card in free-text mode, so the original predicate can be typed back into it')
      .toBe('free-text');

    const narrow = await typeFreeText(page, '${HEIGHT} < 160');
    expect(narrow.invalid, `the free-text card rejected ${'${HEIGHT}'} < 160 as unparseable`).toBe(false);
    expect(await trueCount(page),
      `typed back into a card proved to be carrying no criterion, the original predicate returns exactly the ${height160} rows the form-built rule produced`)
      .toBe(singleRuleCount);
  });

  await softStep('Step 9 — round-trip back to form mode is lossless', async () => {
    const beforeRevert = await trueCount(page);
    expect(beforeRevert, 'the card enters the revert holding the free-text single-rule row set').toBe(singleRuleCount);
    const mode = await toggleExpressionMode(page);
    expect(mode, 'the same toggle returns the card to form mode').toBe('expression');
    expect(await trueCount(page), 'the round-trip revert leaves the row count unchanged').toBe(beforeRevert);
  });

  await softStep('Step 10 — GROK-20242: comma-list paste for regex is a | alternation; trailing comma is a no-op', async () => {
    const alternationMatches = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c = df.col('USUBJID');
      const rex = new RegExp('5|15|25');
      let n = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const raw = c.get(i);
        if (raw === null || raw === undefined) continue;
        if (rex.test(String(raw))) n++;
      }
      return n;
    });
    expect(alternationMatches, 'the alternation must match some rows, or the count comparison proves nothing')
      .toBeGreaterThan(0);
    expect(alternationMatches, 'the alternation must not match every row, or any filter at all would satisfy the check')
      .toBeLessThan(fullCount_demog);

    await resetToExpressionCard(page);
    const noTrail = await pasteRegexValueAndCommit(page, 'USUBJID', '5,15,25');
    expect(noTrail.value, 'GROK-20242: pasted comma list is rewritten to a | alternation').toBe('5|15|25');
    expect(noTrail.rules?.[0], 'the committed rule carries the alternation regex').toContain('5|15|25');
    expect(noTrail.count,
      `the committed alternation must keep exactly the ${alternationMatches} rows whose USUBJID matches 5, 15 or 25, computed independently off the column — a regex matching some-but-not-all rows would otherwise pass`)
      .toBe(alternationMatches);

    await resetToExpressionCard(page);
    const trail = await pasteRegexValueAndCommit(page, 'USUBJID', '5,15,25,');
    expect(trail.value, 'GROK-20242: the trailing comma is trimmed to the same alternation').toBe('5|15|25');
    expect(trail.count, 'appending a trailing comma does not change the row count').toBe(noTrail.count);
  });

  await softStep('Step 11 — the string and date operators build the row sets they name', async () => {
    const cases = [
      {column: 'SEX', op: 'equals', value: 'F'},
      {column: 'RACE', op: 'contains', value: 'an'},
      {column: 'STARTED', op: 'after', value: '01/01/1991'},
    ];
    for (const c of cases) {
      const expected = await page.evaluate(({col, op, val}) => {
        const df = grok.shell.tv.dataFrame;
        const column = df.col(col);
        let n = 0;
        for (let i = 0; i < df.rowCount; i++) {
          const raw = column.get(i);
          if (raw === null || raw === undefined) continue;
          if (op === 'equals') { if (String(raw) === val) n++; }
          else if (op === 'contains') { if (String(raw).toLowerCase().includes(val.toLowerCase())) n++; }
          else {
            const [mm, dd, yyyy] = val.split('/').map(Number);
            const cutoff = new Date(yyyy, mm - 1, dd).getTime();
            const t = raw instanceof Date ? raw.getTime() : new Date(raw).getTime();
            if (!Number.isNaN(t) && t > cutoff) n++;
          }
        }
        return n;
      }, {col: c.column, op: c.op, val: c.value});
      expect(expected, `${c.column} ${c.op} ${c.value} must select a real subset to be a test`)
        .toBeGreaterThan(0);
      expect(expected).toBeLessThan(fullCount_demog);

      await resetToExpressionCard(page);
      await commitFormRule(page, c.column, c.op, c.value);
      expect(await expressionRuleCount(page), `${c.op} committed exactly one rule`).toBe(1);
      expect(await trueCount(page), `${c.column} ${c.op} ${c.value} must select exactly ${expected} rows`)
        .toBe(expected);
      expect(await filterSummaryText(page), 'the summary names the column the rule was built on')
        .toContain(c.column);
    }
  });

  await softStep('Step 12 — unchecking a rule suspends it without removing it', async () => {
    await resetToExpressionCard(page);
    await commitFormRule(page, 'AGE', '>', '50');
    const onlyAge = await trueCount(page);
    await commitFormRule(page, 'HEIGHT', '<', '160');
    expect(await expressionRuleCount(page), 'two rules before the toggle').toBe(2);
    expect(await andOrControl(page, null, false), 'the card starts in OR mode').toBe('OR');
    expect(await andOrControl(page, null, true), 'switch to AND so both rules bite').toBe('AND');
    const bothRules = await trueCount(page);
    expect(bothRules, 'in AND mode two rules are stricter than the AGE rule alone').toBeLessThan(onlyAge);

    const afterUncheck = await toggleRuleCheckbox(page, 1, 2);
    expect(afterUncheck.length, 'the card exposes no per-rule enabled flags to read back').toBe(2);
    expect(afterUncheck[1], 'the second rule must read as switched off').toBe(false);
    expect(await expressionRuleCount(page), 'unchecking must not remove the rule from the list').toBe(2);
    expect(await trueCount(page), 'with the HEIGHT rule off only the AGE rule is left, so its own count returns')
      .toBe(onlyAge);

    const afterRecheck = await toggleRuleCheckbox(page, 1, 2);
    expect(afterRecheck[1], 'the second rule must read as switched back on').toBe(true);
    expect(await trueCount(page), 're-checking the rule restores the exact two-rule count').toBe(bothRules);
  });

  await v.openTable(page, {path: beerPath, withFilterPanel: true});

  let trueCount_or_beer = 0;
  let trueCount_and_beer = 0;

  await softStep('Scenario 2 Step 13 — Aroma text filter: a term drops the beer row count', async () => {
    await resetToAromaTextFilter(page);
    expect(await trueCount(page), 'baseline is the full beer row set').toBe(fullCount_beer);
    await setAromaFuzziness(page, 0);
    const afterMalt = await addAromaTerm(page, 'malt');
    expect(afterMalt, 'typing a term and pressing Enter drops the beer row count').toBeLessThan(fullCount_beer);
    expect(afterMalt).toBeGreaterThan(0);
    expect(await filterSummaryText(page), 'the summary reflects the active Aroma term').toContain('malt');

    const content = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('Aroma');
      let passingWithTerm = 0, passingWithoutTerm = 0, excludedWithTerm = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const has = String(col.get(i) ?? '').toLowerCase().includes('malt');
        if (df.filter.get(i)) has ? passingWithTerm++ : passingWithoutTerm++;
        else if (has) excludedWithTerm++;
      }
      return {passingWithTerm, passingWithoutTerm, excludedWithTerm};
    });
    expect(content.passingWithTerm, 'no row containing the term survived — the term is not the criterion')
      .toBeGreaterThan(0);
    expect(content.passingWithoutTerm, 'rows whose Aroma does not contain the term passed the filter')
      .toBe(0);
    expect(content.excludedWithTerm, 'rows whose Aroma contains the term were filtered out').toBe(0);
  });

  await softStep('Scenario 2 Step 14 — two terms: AND is strictly stricter than OR', async () => {
    await addAromaTerm(page, 'hop');
    expect(await andOrControl(page, 'Aroma', false), 'the two-term card starts in OR mode').toBe('OR');
    trueCount_or_beer = await trueCount(page);
    const afterToggle = await andOrControl(page, 'Aroma', true);
    expect(afterToggle, 'the header toggle switches OR to AND').toBe('AND');
    trueCount_and_beer = await trueCount(page);
    expect(trueCount_and_beer, 'AND is strictly stricter than OR for the same two terms — the control is wired')
      .toBeLessThan(trueCount_or_beer);
  });

  await softStep('Scenario 2 Step 15 — fuzziness 0 yields zero matches; raising the slider raises the count above 0 and grows it', async () => {
    await resetToAromaTextFilter(page);
    expect(await trueCount(page),
      'the fresh Aroma card must be carrying no criterion — without a proven full beer row set here, a leftover filter already sitting at zero rows would satisfy the near-miss check below on its own')
      .toBe(fullCount_beer);
    const atReset = await setAromaFuzziness(page, 0);
    expect(atReset, 'pinning fuzziness to 0 filters nothing on its own').toBe(fullCount_beer);
    const atZero = await addAromaTerm(page, 'maltx');
    expect(atZero, 'at fuzziness 0 a non-matching near-miss term yields zero matches').toBe(0);
    const atMid = await setAromaFuzziness(page, 0.5);
    const atHigh = await setAromaFuzziness(page, 0.8);
    expect(atMid, 'raising fuzziness recovers matches — the count rises above 0').toBeGreaterThan(0);
    expect(atHigh, 'raising fuzziness further grows the matched set').toBeGreaterThan(atMid);
    expect(atHigh).toBeLessThanOrEqual(fullCount_beer);
  });

  await softStep('Teardown', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  v.finishSpec();
});
