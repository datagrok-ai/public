/* What the how-to clips of the Annotate and explain charts section need and the vocabulary lacks: gestures at a
   point given in the units of a scatter plot's axes. A region is drawn around chosen points by a drag whose corners
   are data values, and the region is then pointed at and clicked at a point inside it that no marker covers (a
   marker under the pointer takes precedence over a region). The plot's `view` area and its `x axis min/max`,
   `y axis min/max` readings give the mapping on linear axes. */
import {type Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import {type ElementRef, gestures, viewers} from '@datagrok-libraries/bdd/runtime';

type Point = {x: number; y: number};

/** Lights the `lit` area for the gesture that follows and returns where a data point is on the page. */
async function axesOf(page: Page, target: ElementRef, lit: string): Promise<(x: number, y: number) => Point> {
  const shown = await viewers.hitArea(page, target, lit, true);
  const view = lit === 'view' ? shown : await viewers.hitArea(page, target, 'view');
  const reading = async (name: string): Promise<number> => Number(await viewers.readValue(page, target, name));
  const [x0, x1] = [await reading('x axis min'), await reading('x axis max')];
  const [y0, y1] = [await reading('y axis min'), await reading('y axis max')];
  return (x, y) => ({x: view.x + (x - x0) / (x1 - x0) * view.width, y: view.y + view.height - (y - y0) / (y1 - y0) * view.height});
}

async function dragOn(page: Page, target: ElementRef, x1: number, y1: number, x2: number, y2: number): Promise<void> {
  const at = await axesOf(page, target, 'view');
  const [a, b] = [at(x1, y1), at(x2, y2)];
  await page.mouse.move(a.x, a.y);
  await page.mouse.down();
  await page.mouse.move((a.x + b.x) / 2, (a.y + b.y) / 2, {steps: 4});
  await page.mouse.move(b.x, b.y, {steps: 4});
  await page.mouse.up();
}

export const dragOnAxes = When('user drags across {widget} from {float}, {float} to {float}, {float} on its axes',
  (page: Page, target: ElementRef, x1: number, y1: number, x2: number, y2: number) => dragOn(page, target, x1, y1, x2, y2),
  {tier: 'ui', description: 'a plain drag between two points given as X and Y values of the plot (linear axes): draws a region in the plot\'s drawing mode'});

export const selectOnAxes = When('user drags a selection box on {widget} from {float}, {float} to {float}, {float} on its axes',
  async (page: Page, target: ElementRef, x1: number, y1: number, x2: number, y2: number) => {
    await gestures.withKeys(page, ['Shift'], () => dragOn(page, target, x1, y1, x2, y2));
    await viewers.settle(page, target);
  }, {tier: 'ui', description: 'a Shift-drag between two points given as X and Y values of the plot (linear axes): selects the rows inside'});

export const pointAtOnAxes = When('user points at {float}, {float} inside the {string} area of {widget}',
  async (page: Page, x: number, y: number, area: string, target: ElementRef) => {
    const p = (await axesOf(page, target, area))(x, y);
    await page.mouse.move(p.x, p.y, {steps: 4});
    await viewers.settle(page, target);
  }, {tier: 'ui', description: 'the pointer moved to a point given as X and Y values of the plot, inside the area named (lit for it)'});

export const clickOnAxes = When('user clicks at {float}, {float} inside the {string} area of {widget}',
  async (page: Page, x: number, y: number, area: string, target: ElementRef) => {
    const p = (await axesOf(page, target, area))(x, y);
    await page.mouse.click(p.x, p.y);
    await viewers.settle(page, target);
  }, {tier: 'ui', description: 'a click at a point given as X and Y values of the plot, inside the area named (lit for it)'});
