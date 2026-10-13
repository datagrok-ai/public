/* Drawing in the Ketcher sketcher as a hand does: a template picked from its bottom toolbar, placed on
   the canvas, and the pointer then drifting on over the canvas. While the template tool is on and the
   pointer is over the canvas, Ketcher shows the template's floating preview under it, and that preview
   is part of the structure Ketcher holds. */
import {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';

export const placeBenzene = When('user places the benzene template on the Ketcher canvas', async (page: Page) => {
  const ketcher = page.locator('.d4-dialog .Ketcher-root').last();
  // the sketcher dialog shows before Ketcher has loaded into it
  await ketcher.locator('[data-testid="template-0"]').click({timeout: 30000});
  const canvas = await ketcher.locator('[data-testid="ketcher-canvas"][data-canvasmode="molecules-mode"]').boundingBox();
  if (canvas == null)
    throw new Error('the Ketcher canvas takes no space');
  const x = canvas.x + canvas.width / 2;
  const y = canvas.y + canvas.height / 2;
  await page.mouse.click(x, y);
  await page.mouse.move(x + 8, y + 8, {steps: 3});
}, {tier: 'ui', description: 'the first template of the bottom toolbar (benzene) placed at the centre of the canvas of the open sketcher dialog; the pointer then drifts a few pixels on, still over the canvas'});
