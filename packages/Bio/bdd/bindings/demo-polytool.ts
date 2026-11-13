/* Steps of the PolyTool demo clip (features/demos/19): a click on a monomer of the HELM drawn in a dialog. */
import type {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

/** A HELM input draws its molecule as an SVG that takes no pointer events (`pointer-events: none`): the click goes to
 * the input's host, which finds the monomer nearest to the point (PolyTool's Enumerate HELM adds that position to its
 * placeholders). The point is the centre of the monomer's shape, `.hw-atom[data-atom-id]` of the hwe editor (the
 * monomer's continuous position from 1), read once the drawing holds still: the editor redraws while its dialog
 * sizes it. */
export const clickHelmMonomer = When('user clicks on monomer {int} of the HELM in the {string} dialog',
  async (page: Page, position: number, title: string) => {
    let point: {x: number, y: number} | string = '';
    await expect.poll(async () => {
      point = await page.evaluate(async ([t, p]) => {
        const dialog = [...document.querySelectorAll('.d4-dialog')].reverse()
          .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === t);
        if (!dialog)
          return `no "${t}" dialog`;
        const box = () => {
          const atom = dialog.querySelector(`.ui-input-helm .hw-atom[data-atom-id="${p}"]`);
          const r = atom?.getBoundingClientRect();
          return r && r.width > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
        };
        const before = box();
        await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
        const after = box();
        if (!before || !after)
          return `no monomer ${p} drawn in the "${t}" dialog`;
        return Math.abs(before.x - after.x) < 0.5 && Math.abs(before.y - after.y) < 0.5 ? after : 'the drawing moves';
      }, [title, position] as [string, number]);
      return typeof point === 'string' ? point : 'still';
    }, {message: `monomer ${position} of the HELM in the "${title}" dialog`}).toBe('still');
    const {x, y} = point as unknown as {x: number, y: number};
    await page.mouse.move(x - 3, y - 3);
    await page.mouse.move(x, y);
    await page.mouse.click(x, y);
  }, {tier: 'ui', description: 'the pointer on the centre of that monomer\'s shape once the drawing holds still, then a click there'});
