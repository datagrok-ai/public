/* Steps of the structure demo clips (features/demos/05-08): the molecule a context panel pane draws, its
   full-screen dialog and the magnifying lens over it, and the 3D view turned with the pointer. */
import type {Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {expect, pollMs, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The current row is drawn green, which tints the molecules a demo colors by monomer. */
export const noCurrentRow = Given('no row of the table is current', async (page: Page) => {
  silent(page);
  await page.evaluate(() => { grok.shell.tv.dataFrame.currentRowIdx = -1; });
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx), {message: 'the current row'}).toBe(-1);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the current row of the current view\'s table cleared'});

/** A context panel pane by its title (the Dart accordion names its header `div-section--<Title>`, spaces dashed),
 * a nested one too (Oligo Structures > Sense). */
function pane(page: Page, title: string): Locator {
  return page.locator(`.d4-accordion-pane:has(> [name="div-section--${title.replace(/ /g, '-')}"])`).filter({visible: true}).last();
}

/** A pane of an accordion inside a context panel pane (SequenceTranslator's Oligo Structures > Sense): its own
 * header is clicked. The library's `user expands "Sense" section` passes there without expanding it. */
export const expandInnerPane = When('user expands the {string} pane inside the {string} section',
  async (page: Page, inner: string, outer: string) => {
    const header = pane(page, outer).locator(`[name="div-section--${inner.replace(/ /g, '-')}"]`).filter({visible: true}).first();
    await expect(header, `the ${inner} pane inside the ${outer} section`).toBeVisible();
    if (await header.getAttribute('aria-expanded') !== 'true')
      await header.click();
    await expect(header, `the ${inner} pane inside the ${outer} section`).toHaveAttribute('aria-expanded', 'true');
  }, {tier: 'ui', description: 'a click on the header of the inner pane unless it is expanded; done when it says it is'});

/** Bio's Molecular Structure pane draws the molecule on a canvas that opens it full screen on a click. */
export const clickPaneMolecule = When('user clicks on the molecule in the {string} section', async (page: Page, title: string) => {
  const canvas = pane(page, title).locator('canvas').filter({visible: true}).first();
  await expect(canvas, `the molecule of the ${title} section`).toBeVisible();
  await canvas.click();
  await expect(page.locator('[name="dialog-Molecule"] canvas').filter({visible: true}), 'the full-screen molecule')
    .toHaveCount(1, {timeout: pollMs(30000)});
}, {tier: 'ui', description: 'a click on the drawing of the section (Click to expand); done when the full-screen Molecule dialog has its drawing'});

/** The full-screen dialog magnifies the part of the molecule under the pointer in a lens centred on it. */
export const moveLens = When('user hovers over the molecule in the {string} dialog at {int}% across and {int}% down',
  async (page: Page, title: string, across: number, down: number) => {
    const dialog = page.locator(`[name="dialog-${title}"]`).filter({visible: true});
    const drawing = dialog.locator('canvas').filter({visible: true}).first();
    const box = await drawing.boundingBox();
    if (!box)
      throw new Error(`the ${title} dialog shows no molecule`);
    const x = box.x + box.width * across / 100;
    const y = box.y + box.height * down / 100;
    await page.mouse.move(x, y, {steps: 8});
    await expect.poll(() => dialog.evaluate((el, [px, py]) => {
      const lens = Array.from(el.querySelectorAll('canvas')).find((c) => getComputedStyle(c).pointerEvents === 'none');
      if (!lens || lens.style.display === 'none')
        return 'no lens shown';
      const r = lens.getBoundingClientRect();
      const off = Math.hypot(r.x + r.width / 2 - px, r.y + r.height / 2 - py);
      return off <= 2 ? 'under the pointer' : `${Math.round(off)} px from the pointer`;
    }, [x, y]), {message: `the lens of the ${title} dialog`}).toBe('under the pointer');
  }, {tier: 'ui', description: 'the pointer on that point of the drawing; done when the magnifying lens is centred on it'});

/** A drawn canvas: a 2D one has ink that is not the white of its paper; a WebGL one (NGL) cannot be read back,
 * so its picture is taken twice and must hold more than a blank square and have stopped changing. */
export const paneShowsDrawing = Then('the {string} section should show a drawn molecule', async (page: Page, title: string) => {
  const canvas = pane(page, title).locator('canvas').filter({visible: true}).first();
  await expect(canvas, `the drawing of the ${title} section`).toBeVisible({timeout: pollMs(60000)});
  await expect(pane(page, title), `the ${title} section`).not.toContainText('Couldn\'t get');
  let last: Buffer | null = null;
  await expect.poll(async () => {
    const ink = await canvas.evaluate((c: HTMLCanvasElement) => {
      const g = c.getContext('2d');
      if (!g)
        return -1;
      const d = g.getImageData(0, 0, c.width, c.height).data;
      let n = 0;
      for (let i = 0; i < d.length; i += 4)
        if (d[i + 3] > 128 && (d[i] < 200 || d[i + 1] < 200 || d[i + 2] < 200))
          n++;
      return n;
    });
    if (ink >= 0)
      return ink > 200 ? 'drawn' : `${ink} inked pixels`;
    const shot = await canvas.screenshot();
    const steady = last != null && shot.equals(last);
    last = shot;
    return shot.length > 12000 && steady ? 'drawn' : `a ${shot.length} byte picture${steady ? '' : ', still changing'}`;
  }, {message: `the molecule of the ${title} section`, timeout: pollMs(60000)}).toBe('drawn');
}, {description: 'the first canvas of the section is drawn: ink on a 2D canvas, a steady non-blank picture of a WebGL one; no "Couldn\'t get" message'});

/** SequenceTranslator's oligo duplex renderer keeps itself in the grid column's temp (`oligoNucleotideRenderer`)
 * and answers where it drew a monomer (`getMonomerPoint`: a nucleotide's or a conjugate's chip) or the arc of the
 * 3' linkage after it (`getLinkagePoint`), px from the cell's top left; monomers count from 1 over the sense strand,
 * then the antisense one. Done when the tooltip shows with its structures drawn. */
async function hoverDuplex(page: Page, target: ElementRef, area: string, monomer: number, linkage: boolean): Promise<void> {
  const cellOf = /^cell (\d+) of (.+)$/.exec(area);
  if (!cellOf)
    throw new Error(`"${area}" is no grid cell: a duplex is hovered in "cell <row> of <column>"`);
  const cell = await viewers.hitArea(page, target, area, true);
  const point = await viewers.onViewer(page, target, (el, [row, column, pos, link]) => {
    const grid = (window as any).__bdd.viewerOf(el);
    const gridCell = grid.cell(column, grid.tableRowToGrid(row - 1));
    const back = gridCell.gridColumn.temp.oligoNucleotideRenderer;
    if (!back)
      return 'it is not drawn as an oligo duplex';
    return link ? back.getLinkagePoint(gridCell, pos - 1) : back.getMonomerPoint(gridCell, pos - 1);
  }, [Number(cellOf[1]), cellOf[2], monomer, linkage]) as {x: number, y: number} | string | null;
  if (point === null || typeof point === 'string')
    throw new Error(`${area} of ${target.phrase} draws no ${linkage ? 'linkage after ' : ''}monomer ${monomer}${point ? `: ${point}` : ''}`);
  await page.mouse.move(cell.x + point.x - 3, cell.y + point.y - 3);
  await page.mouse.move(cell.x + point.x, cell.y + point.y);
  const tooltip = page.locator('.d4-tooltip').filter({visible: true});
  await expect(tooltip, `the tooltip of ${linkage ? 'the linkage after ' : ''}monomer ${monomer}`).toHaveCount(1);
  await expect(tooltip, 'the structures of the tooltip').not.toContainText('Loading structures', {timeout: pollMs(30000)});
}

export const hoverDuplexMonomer = When('user hovers over monomer {int} of the duplex in the {string} area of {widget}',
  (page: Page, monomer: number, area: string, target: ElementRef) => hoverDuplex(page, target, area, monomer, false),
  {tier: 'ui', description: 'a nucleotide\'s or a conjugate\'s chip of an oligo duplex cell, counted from 1 over the sense strand then the antisense one; done when its tooltip has drawn its structures'});

export const hoverDuplexLinkage = When('user hovers over the linkage after monomer {int} of the duplex in the {string} area of {widget}',
  (page: Page, monomer: number, area: string, target: ElementRef) => hoverDuplex(page, target, area, monomer, true),
  {tier: 'ui', description: 'the arc of the 3\' linkage after that monomer of an oligo duplex cell; done when its tooltip has drawn the linkage\'s structure'});

/** NGL turns the structure as the pointer drags across it. */
export const turn3d = When('user turns the 3D structure in the {string} section', async (page: Page, title: string) => {
  const canvas = pane(page, title).locator('.d4-ngl-viewer canvas').filter({visible: true}).first();
  const box = await canvas.boundingBox();
  if (!box)
    throw new Error(`the ${title} section shows no 3D view`);
  const y = box.y + box.height * 0.5;
  await page.mouse.move(box.x + box.width * 0.2, y);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.8, y + box.height * 0.08, {steps: 16});
  await page.mouse.up();
}, {tier: 'ui', description: 'a drag across the middle of the 3D view, left to right'});
