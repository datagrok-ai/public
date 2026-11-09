/* Steps of the antibody demo clips (features/demos/09-11b): off camera, the window of positions a sequence
   column's header shows, a width the header does not reset, balloons and a stale tooltip out of the way; on
   camera, the header dragged along its ruler and its region band under the pointer. */
import type {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {expect, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The header and the cells of a long sequence column show a window of its positions, from the column's
 * `.positionShift` tag (Bio's sequence-scrolling widget), named by its position names or numbered from 1;
 * a numbered column starts in the N-terminal flank, so a clip moves it to the scheme's first position. */
export const headerStartsAt = When('the header of the {string} column starts at position {string}',
  async (page: Page, column: string, position: string) => {
    silent(page);
    await page.evaluate(([c, p]) => {
      const col = grok.shell.tv?.dataFrame.col(c);
      if (!col)
        throw new Error(`the table has no "${c}" column`);
      const tag: string = col.getTag('.positionLabels') || col.getTag('.positionNames') || '';
      const names = tag ? tag.split(',').map((n: string) => n.trim()) : [];
      const index = names.length ? names.indexOf(p) : Number(p) - 1;
      if (!(index >= 0))
        throw new Error(`"${c}" has no position "${p}"; it starts: ${names.slice(0, 30).join(', ')}`);
      col.setTag('.positionShift', String(index));
      grok.shell.tv.grid.invalidate();
    }, [column, position] as [string, string]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a filmed demo leaves out: the column\'s position window starts at the named position'});

/** Bio's header sets a new long sequence column 400 px wide 300 ms after it first paints the column's header
 * (`sequence-scrolling-widget.ts`), over a width set meanwhile: the width is set once the header shows
 * and set again whenever it changes in the next 400 ms. */
export const gridColumnStaysWide = When('the {string} column of the grid stays {int} px wide',
  async (page: Page, column: string, width: number) => {
    silent(page);
    await page.evaluate(async ([c, w]) => {
      const col = grok.shell.tv?.grid?.columns.byName(c);
      if (!col)
        throw new Error(`the grid of the current view has no "${c}" column`);
      const until = Date.now() + 400;
      while (Date.now() < until) {
        if (col.width !== w)
          col.width = w;
        await new Promise((r) => requestAnimationFrame(r));
      }
    }, [column, width] as [string, number]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a filmed demo leaves out: the column\'s width, held past the header\'s own 400 px default'});

/** A sequence header is dragged along its positions like a strip (a position per 20 px the pointer moves);
 * on the ruler (the position names, the band above the header's slider) no track shows a tooltip on the
 * way. The ruler's labels sit 15 px above the 8 px slider at the bottom of the grid's column header. A
 * drag, not the wheel: the header moves a position per wheel event, and a filmed wheel turn is split into
 * notches. */
export const dragRuler = When('user drags the ruler under the {string} area of {widget} by {int} pixels to the {word}',
  async (page: Page, area: string, target: ElementRef, px: number, direction: string) => {
    // the WebLogo area only places the press: it is not lit, the press lights the header it lands on
    await viewers.onViewer(page, target, (el) => (window as any).__bdd.settle(el, 300));
    const logo = await viewers.hitArea(page, target, area);
    const header = await viewers.onViewer(page, target, (el) => {
      const grid = (window as any).__bdd.viewerOf(el);
      return {top: grid.canvas.getBoundingClientRect().top, height: grid.props.colHeaderHeight};
    });
    await viewers.dragFrom(page, {x: logo.x + logo.width / 2, y: header.top + header.height - 8 - 15},
      viewers.dragDelta(px, direction));
  }, {tier: 'ui', description: 'pressed on the position names under that WebLogo stack, moved by the distance, released'});

/** The region band (Bio's annotation track) is the lowest track of the header, right under the WebLogo
 * whose `position <name> of <column> header` areas the grid reports: the pointer goes into the band
 * under that position, and the region's tooltip (name, scheme, range, position) shows. */
export const hoverRegionTrack = When('user hovers over the region track under the {string} area of {widget}',
  async (page: Page, area: string, target: ElementRef) => {
    const logo = await viewers.hitArea(page, target, area, true);
    const x = logo.x + logo.width / 2;
    const y = logo.y + logo.height + 4 + 7;
    await page.mouse.move(x - 3, y - 3);
    await page.mouse.move(x, y);
    await expect(page.locator('.d4-tooltip').filter({visible: true}), `the tooltip of the region under ${area}`).toHaveCount(1);
  }, {tier: 'ui', description: 'the annotation track 4 px under the WebLogo stack of that position; done when its tooltip shows'});
