import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {knownOpenBug} from '../../helpers/known-open-bug';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Calendar';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Calendar';
const datasetPath = 'System:DemoFiles/demog.csv';
const DATE_COLUMN = 'STARTED';

const RED: v.RgbRange = {rMin: 180, rMax: 255, gMin: 0, gMax: 90, bMin: 0, bMax: 90};

const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);

const selectionCount = (page: Page) =>
  page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);

const clearSelection = async (page: Page) => {
  await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  await expect.poll(() => selectionCount(page), {timeout: 5000}).toBe(0);
};

async function regions(page: Page, withHeader = true) {
  const box = await page.evaluate(() => {
    const canvas = document.querySelector('[name="viewer-Calendar"] canvas[name="canvas"]') as HTMLCanvasElement;
    const r = canvas.getBoundingClientRect();
    return {x: r.x, y: r.y, w: r.width, h: r.height};
  });
  const dy = (box.w - 5) / 8;
  const monthWidth = dy + 5;
  const headerH = withHeader ? dy * 1.1 : 0;
  const daysLeft = box.x + monthWidth;
  const dayWidth = (box.w - monthWidth) / 7;
  return {
    box, dy, monthWidth, headerH,

    weekday: (i: number) => ({x: daysLeft + dayWidth * (i + 0.5), y: box.y + headerH + dy / 2}),

    month: (week = 2.5) => ({x: box.x + monthWidth / 2, y: box.y + headerH + dy * (week + 1)}),

    day: (col: number, week: number) =>
      ({x: daysLeft + dayWidth * (col + 0.5), y: box.y + headerH + dy * (week + 1.5)}),
  };
}

interface Hit { caption: string; rows: number; action: string }

async function hover(page: Page, at: {x: number; y: number}): Promise<Hit | null> {
  await page.mouse.move(at.x - 25, at.y - 25);
  await page.mouse.move(at.x, at.y, {steps: 4});
  const text = await page.evaluate(async () => {
    for (let i = 0; i < 20; i++) {
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      if (t && getComputedStyle(t).display !== 'none' && t.innerText.trim().length > 0)
        return t.innerText;
      await new Promise((r) => setTimeout(r, 100));
    }
    return '';
  });
  const m = text.match(/^(.+?)\s*\n\s*Click to (\w+)\s*\n\s*(\d+) rows/);
  return m ? {caption: m[1].trim(), action: m[2], rows: Number(m[3])} : null;
}

async function findDayCell(page: Page): Promise<{hit: Hit; at: {x: number; y: number}}> {
  const r = await regions(page);
  for (const week of [2, 3, 4, 1, 5, 6])
    for (const col of [2, 3, 4, 1, 5]) {
      const at = r.day(col, week);
      const hit = await hover(page, at);
      if (hit && hit.rows > 0) return {hit, at};
    }
  throw new Error('no day cell with rows found on the calendar');
}

function dateCounts(page: Page, column: string) {
  return page.evaluate((col) => {
    const t = (window as any).grok.shell.t;
    const c = t.col(col);
    const byDay: Record<string, number> = {};
    const byMonth: Record<string, number> = {};
    const byWeekday: number[] = [0, 0, 0, 0, 0, 0, 0];
    for (let i = 0; i < t.rowCount; i++) {
      const value = c.get(i);
      if (value == null) continue;
      const day = value.format ? value.format('YYYY-MM-DD') : String(value).slice(0, 10);
      byDay[day] = (byDay[day] ?? 0) + 1;
      byMonth[day.slice(0, 7)] = (byMonth[day.slice(0, 7)] ?? 0) + 1;
      byWeekday[new Date(day).getUTCDay()]++;
    }
    return {byDay, byMonth, byWeekday};
  }, column);
}

const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June',
  'July', 'August', 'September', 'October', 'November', 'December'];
const WEEKDAYS = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'];

test('Calendar', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add Calendar from the Viewers toolbox', async () => {
    await page.locator('[name="icon-calendar"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE)).total,
      {timeout: 30_000}).toBeGreaterThan(1000);
  });

  await softStep('Day, month and weekday tooltips report the matching rows', async () => {
    const counts = await dateCounts(page, DATE_COLUMN);

    const {hit: day} = await findDayCell(page);

    expect(day.caption).toMatch(/^\d{4}-\d{2}-\d{2}/);
    expect(day.rows).toBe(counts.byDay[day.caption.slice(0, 10)]);

    const r = await regions(page);
    const month = await hover(page, r.month());
    expect(month).not.toBeNull();
    const [monthName, year] = month!.caption.split(' ');
    expect(MONTHS).toContain(monthName);
    const monthKey = `${year}-${String(MONTHS.indexOf(monthName) + 1).padStart(2, '0')}`;
    expect(month!.rows).toBe(counts.byMonth[monthKey]);

    for (let i = 1; i < 7; i++) {
      const weekday = await hover(page, r.weekday(i));
      expect(weekday).not.toBeNull();
      expect(weekday!.caption).toBe(WEEKDAYS[i]);
      expect(weekday!.rows).toBe(counts.byWeekday[i]);
    }

    const sunday = await hover(page, r.weekday(0));
    expect(sunday).not.toBeNull();
    expect(sunday!.caption).toBe(WEEKDAYS[0]);
    await knownOpenBug('GROK-20634', () => {
      expect(sunday!.rows).toBe(counts.byWeekday[0]);
    });
  });

  await softStep('The auto-picked date column is reflected in the property grid', async () => {

    await knownOpenBug('GROK-20636', async () => {
      expect(await page.evaluate(() => (window as any).grok.shell.tv.viewers
        .find((x: any) => x.type === 'Calendar').props.dateColumnName)).toBe(DATE_COLUMN);
    });
  });

  await softStep('Clicking a day, a month and a weekday selects their rows', async () => {
    await clearSelection(page);
    const {hit: day, at} = await findDayCell(page);
    await page.mouse.click(at.x, at.y);
    await expect.poll(() => selectionCount(page), {timeout: 8000}).toBe(day.rows);

    const r = await regions(page);
    await clearSelection(page);
    const month = await hover(page, r.month());
    await page.mouse.click(r.month().x, r.month().y);
    await expect.poll(() => selectionCount(page), {timeout: 8000}).toBe(month!.rows);

    await clearSelection(page);
    const monday = await hover(page, r.weekday(1));
    await page.mouse.click(r.weekday(1).x, r.weekday(1).y);
    await expect.poll(() => selectionCount(page), {timeout: 8000}).toBe(monday!.rows);
  });

  await softStep('Shift+click extends the selection, Ctrl+click toggles it', async () => {
    await clearSelection(page);
    const r = await regions(page);
    const monday = await hover(page, r.weekday(1));
    await page.mouse.click(r.weekday(1).x, r.weekday(1).y);
    await expect.poll(() => selectionCount(page), {timeout: 8000}).toBe(monday!.rows);

    const tuesday = await hover(page, r.weekday(2));
    await page.keyboard.down('Shift');
    await page.mouse.click(r.weekday(2).x, r.weekday(2).y);
    await page.keyboard.up('Shift');
    await expect.poll(() => selectionCount(page), {timeout: 8000})
      .toBe(monday!.rows + tuesday!.rows);

    await page.keyboard.down('Control');
    await page.mouse.click(r.weekday(2).x, r.weekday(2).y);
    await page.keyboard.up('Control');
    await expect.poll(() => selectionCount(page), {timeout: 8000}).toBe(monday!.rows);
    await clearSelection(page);
  });

  await softStep('On Click set to Filter makes the same click filter instead', async () => {
    await v.openViewerProperties(page, VIEWER_NAME, '[name="prop-on-click"]');
    await category(page, 'data', 'on-click');
    await v.selectPropertyGridChoice(page, 'on-click', 'Filter', 'data');
    expect(await v.propertyGridValue(page, 'on-click', 'data')).toBe('Filter');

    const r = await regions(page);
    const month = await hover(page, r.month());

    expect(month!.action).toBe('filter');

    await page.mouse.click(r.month().x, r.month().y);
    await expect.poll(() => page.evaluate(() => (window as any).grok.shell.t.filter.trueCount),
      {timeout: 8000}).toBe(month!.rows);
    expect(await selectionCount(page)).toBe(0);

    await v.resetFilters(page);
    await category(page, 'data', 'on-click');
    await v.selectPropertyGridChoice(page, 'on-click', 'Select', 'data');
    expect(await v.propertyGridValue(page, 'on-click', 'data')).toBe('Select');
  });

  await softStep('Red Weekends colours the weekend day numbers', async () => {
    await category(page, 'misc', 'red-weekends');
    expect((await v.countCanvasPixels(page, VIEWER_TYPE, {rgbRange: RED})).matched).toBeGreaterThan(0);

    await v.setPropertyGridCheckbox(page, 'red-weekends', false, 'misc');
    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE, {rgbRange: RED})).matched,
      {timeout: 10_000}).toBe(0);

    await v.setPropertyGridCheckbox(page, 'red-weekends', true, 'misc');
    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE, {rgbRange: RED})).matched,
      {timeout: 10_000}).toBeGreaterThan(0);
  });

  await softStep('Show Header hides the year caption and gives the days its space', async () => {
    await category(page, 'misc', 'show-header');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    const withHeader = (await v.countCanvasPixels(page, VIEWER_TYPE)).total;
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.setPropertyGridCheckbox(page, 'show-header', false, 'misc');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    expect((await v.countCanvasPixels(page, VIEWER_TYPE)).total).toBeGreaterThan(withHeader);

    const shifted = await regions(page, false);
    const weekday = await hover(page, shifted.weekday(1));
    expect(weekday!.caption).toBe('Monday');

    await v.setPropertyGridCheckbox(page, 'show-header', true, 'misc');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});
  });

  await softStep('Filtering the table reshapes the calendar', async () => {
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeLessThan(total);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 200});

    const counts = await dateCounts(page, DATE_COLUMN);
    const {hit} = await findDayCell(page);
    expect(hit.rows).toBe(counts.byDay[hit.caption.slice(0, 10)]);

    await category(page, 'misc', 'show-filtered-only');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.setPropertyGridCheckbox(page, 'show-filtered-only', false, 'misc');
    await knownOpenBug('GROK-20635', async () => {
      await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 200, timeoutMs: 8000});
    });
    await v.setPropertyGridCheckbox(page, 'show-filtered-only', true, 'misc');

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.resetFilters(page);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 200});
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
