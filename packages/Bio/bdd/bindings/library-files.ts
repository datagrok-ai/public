/* The readiness of the Manage Monomers view's sketcher. */
import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

/** The Manage Monomers view hosts a monomer editor with a sketcher that mounts seconds after the
 * view (Ketcher: about ten on dev, over a minute on a stand serving a second worker). Closed
 * before that, the sketcher's late mount throws ResizeObserver TypeErrors into whatever runs
 * next, so a feature that opens the view lets the editor finish first — a Ketcher toolbar or,
 * for another backend, its canvas. */
export const monomerSketcherReady = Then('the monomer sketcher of the Manage Monomers view should be ready', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const root = document.querySelector('.monomer-manager-sketcher');
    if (!root)
      return 'no sketcher';
    const ketcher = root.querySelector('.Ketcher-root');
    return ketcher ? (ketcher.querySelectorAll('button').length > 5 ? 'ready' : 'mounting') : (root.querySelector('canvas') ? 'ready' : 'mounting');
  }), {timeout: pollMs(150000), message: 'the monomer editor sketcher'}).toBe('ready');
}, {description: 'the sketcher inside .monomer-manager-sketcher has its toolbar (Ketcher) or its canvas (other backends)'});
