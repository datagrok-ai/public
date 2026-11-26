/* What the how-to clips of the Viewers section need and the vocabulary lacks: a viewer docked beside the
   grid is given a larger share of the view (set-up, written as a Given so it is off camera), and the
   aggregation a viewer's on-chart value selector offers next to the column (count, avg, sum...) is picked
   the way a person does, from its list. */
import {type Page} from '@playwright/test';
import {Given, When, element} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, gestures, viewers} from '@datagrok-libraries/bdd/runtime';

// the platform names none: a bare <input placeholder="Search"> (ui.searchInput) above the property grid
element('context panel search', {selector: '.grok-prop-panel input.d4-search-input', aliases: ['settings search'],
  description: 'the Search box above a viewer\'s properties in the context panel: it keeps the properties whose row holds the text, and opens their sections'});

/** dock-spawn puts a splitbar right before the panel docked beside or under another one (`.splitbar-vertical`
 * to the right of the grid, `.splitbar-horizontal` under it, where a line chart goes), and its container is
 * the space the two share; the bar is dragged to where the viewer gets the share asked for. */
export const viewerShare = Given('{widget} takes {int}% of the space it shares with the grid',
  async (page: Page, target: ElementRef, percent: number) => {
    const viewer = await viewers.viewerLocator(page, target);
    const panel = viewer.locator('xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]');
    // the panel itself, or the column or row of panels it sits in, comes right after the bar
    const isBar = 'contains(concat(" ", normalize-space(@class), " "), " splitbar-vertical ") or ' +
      'contains(concat(" ", normalize-space(@class), " "), " splitbar-horizontal ")';
    const bar = panel.locator(`xpath=ancestor-or-self::*[preceding-sibling::*[1][${isBar}]][1]/preceding-sibling::*[1]`);
    if (await bar.count() === 0)
      throw new Error(`${target.phrase} has no dock splitter between it and the grid`);
    const across = (await bar.getAttribute('class') ?? '').includes('splitbar-vertical');
    const barBox = await bar.boundingBox();
    const shared = await bar.locator('xpath=..').boundingBox();
    if (!barBox || !shared)
      throw new Error(`the dock splitter of ${target.phrase} is not on screen`);
    const from = {x: barBox.x + barBox.width / 2, y: barBox.y + barBox.height / 2};
    const to = across ? {x: shared.x + shared.width * (1 - percent / 100), y: from.y} :
      {x: from.x, y: shared.y + shared.height * (1 - percent / 100)};
    await page.mouse.move(from.x, from.y);
    await page.mouse.down();
    await page.mouse.move((from.x + to.x) / 2, (from.y + to.y) / 2, {steps: 4});
    await page.mouse.move(to.x, to.y, {steps: 4});
    await page.mouse.up();
    // off the bar, which stays lit while the pointer is on it: over the grid, which shows no tooltip
    await page.mouse.move(across ? to.x - 200 : to.x, across ? to.y : to.y - 100, {steps: 2});
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'drags the dock splitter between the grid and the viewer until the viewer takes that share of the space they share'});

/** The column selector drawn on a viewer (`div-column-combobox-<property>`) carries the aggregation as a native
 * `<select>` in front of the column name ("count ⌄ AGE ⌄" under a bar chart). */
export const pickAggregation = When('user picks the {string} aggregation in the {string} column selector of {widget}',
  async (page: Page, aggregation: string, selector: string, target: ElementRef) => {
    const viewer = await viewers.viewerLocator(page, target);
    const select = viewer.locator(`[name="div-column-combobox-${selector}"] select`).first();
    await gestures.selectNative(page, select, aggregation);
    await expect(select, `the aggregation of the "${selector}" selector of ${target.phrase}`).toHaveValue(aggregation);
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'picks from the aggregation list in front of the column of an on-viewer column selector'});

/** A column picker some other step opened (the + of a pivot table's row): the name typed into its search and
 * the row it leaves clicked, as `pickInColumnGrid` films it. */
export const pickInOpenColumnPicker = When('user picks {string} in the open column picker', async (page: Page, column: string) => {
  await gestures.pickInColumnGrid(page, column, 'the open column picker');
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'types the column into the search of the column picker already open and takes it'});

/** A line chart with several Y columns in one chart (Multi Axis) and a split shows "Y Aggr Type:" above the plot,
 * one native <select> for every line (`d4-line-chart-aggr-type-selector`, bound to Aggr Type). */
export const pickLineChartAggregation = When('user picks {string} in the Y Aggr Type selector of {widget}',
  async (page: Page, aggregation: string, target: ElementRef) => {
    const viewer = await viewers.viewerLocator(page, target);
    const select = viewer.locator('select.d4-line-chart-aggr-type-selector').first();
    await gestures.selectNative(page, select, aggregation);
    await expect(select, `the Y Aggr Type selector of ${target.phrase}`).toHaveValue(aggregation);
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'picks from the "Y Aggr Type" list a multi-axis line chart with a split shows above its plot'});
