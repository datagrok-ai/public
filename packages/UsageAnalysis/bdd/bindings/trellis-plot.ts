/* The steps only the trellis plot needs. Everything else its features use is the library's
   `viewers` tier and the platform's data steps, read from the hit areas and readings the trellis
   reports (`cell F | Caucasian`, `cell body …`, `x plus`, `x range slider max handle`, `cells`,
   `cells drawn`, `distinct cell signatures`, `cell signature <cell>` and the rest — see
   `core/client/d4/lib/src/viewers/trellis_plot/CLAUDE.md`).

   What is left here: the inner viewer's own look, which is a property OF a property
   (`innerViewerLook`) and so has no caption the library's property steps can name; the two-number
   geometry read the features make a dozen times; and the four paging icons, whose `aria-disabled`
   is a state of an element, not of a hit area. */
import {Page} from '@playwright/test';
import {Then, When, element} from '@datagrok-libraries/bdd';
import {ElementRef, exactText, viewers} from '@datagrok-libraries/bdd/runtime';

// --- the paging icons ---------------------------------------------------------------------------

/* The (+)/(-) icons that page one category row of an axis in and out. They are hit areas as well
   (`x plus`, so a click lands on them by name), but "disabled" is an attribute of the element the
   trellis writes it on. */
element('x plus icon', {selector: '[name="x-axis-icons"] .d4-trellis-plot-add-cat-icon'});
element('x minus icon', {selector: '[name="x-axis-icons"] .d4-trellis-plot-remove-cat-icon'});
element('y plus icon', {selector: '[name="y-axis-icons"] .d4-trellis-plot-add-cat-icon'});
element('y minus icon', {selector: '[name="y-axis-icons"] .d4-trellis-plot-remove-cat-icon'});

// --- the inner viewer's look --------------------------------------------------------------------

const coerce = (text: string): unknown => {
  if (/^(true|false)$/i.test(text))
    return /^true$/i.test(text);
  if (text !== '' && !Number.isNaN(Number(text)))
    return Number(text);
  return text;
};

/** A setting of the viewer inside every cell: `viewer.setOptions({innerViewerLook: {...}})`, which
 * is the only channel that reaches the cells — a mutation of the look object itself is not read
 * back by them. The name is the look's own field (`allowZoom`, `xColumnName`, `colorColumnName`),
 * as `viewer.getOptions().look.innerViewerLook` spells it. */
export const setInnerProperty = When('user sets {string} inner property of {widget} to {string}',
  async (page: Page, name: string, target: ElementRef, value: string) => {
    const loc = await viewers.viewerLocator(page, target);
    await viewers.snapshot(page, target);
    await loc.evaluate((el, [n, v]) => {
      (window as any).__bdd.viewerOf(el).setOptions({innerViewerLook: {[n]: v}});
    }, [name, coerce(value)] as [string, unknown]);
    await loc.evaluate((el) => (window as any).__bdd.settle(el, 3000));
  }, {tier: 'api', description: 'the inner viewer\'s own look, by the field name the look serializes (allowZoom, xColumnName, colorColumnName)'});

// --- geometry -------------------------------------------------------------------------------------

/** The viewport in cells, from the `columns` and `rows` readings: the two numbers the category,
 * tiling and scrolling scenarios claim together. */
export const cellsWideTall = Then('the cells of {widget} should be {int} wide and {int} tall',
  async (page: Page, target: ElementRef, wide: number, tall: number) => {
    await viewers.expectReading(page, target, 'columns', 'equal', wide);
    await viewers.expectReading(page, target, 'rows', 'equal', tall);
  }, {description: 'the viewport of cells on screen — "columns" by "rows"'});

// --- the inner viewer type, through the control panel ----------------------------------------------

/** The type selector of the control panel (a `ComboPopup` named "viewer selector"): a click opens
 * the list of trellisable viewers, and the pick goes through `setViewerType`, which is what
 * announces `d4-trellis-plot-viewer-type-changed` — writing the property does not. */
export const pickInnerViewer = When('user picks {string} in the viewer selector of {widget}',
  async (page: Page, type: string, target: ElementRef) => {
    const loc = await viewers.viewerLocator(page, target);
    const selector = loc.locator('[name="viewer selector"]').first();
    await selector.waitFor({state: 'visible', timeout: 5000});
    await viewers.snapshot(page, target);
    await selector.click();
    const popup = page.locator('.d4-combo-popup-expanded').last();
    await popup.waitFor({state: 'visible', timeout: 5000});
    const item = popup.locator('.d4-list-item').filter({hasText: exactText(type)}).first();
    if (await item.count() === 0)
      throw new Error(`no "${type}" in the viewer selector; it offers: ${(await popup.locator('.d4-list-item').allTextContents()).join(', ')}`);
    await item.click();
    await loc.evaluate((el) => (window as any).__bdd.settle(el, 3000));
  }, {tier: 'ui', description: 'the control panel must be shown; the pick is the gesture that fires the type-changed event'});
