/* Steps of the HELM demo clip (features/demos/03b): a pane of an accordion nested in another one, as the editor's
   Molecular Structure tab lays out the Explore pane (Chem's panels for the whole molecule: Chemistry > Properties). */
import type {Locator, Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

/** An accordion pane by its title (the Dart accordion names its header `div-section--<Title>`, spaces dashed), the
 * innermost one shown when several are. */
function pane(page: Page, title: string): Locator {
  return page.locator(`.d4-accordion-pane:has(> [name="div-section--${title.replace(/ /g, '-')}"])`).filter({visible: true}).last();
}

/** The library's `user expands "Chemistry" section` finds a nested pane twice (the pane and its header) and passes
 * without opening it: the header of the inner pane is clicked, inside the outer one, unless it says it is open. */
export const expandInnerPane = When('user expands the {string} pane inside the {string} section',
  async (page: Page, inner: string, outer: string) => {
    const header = pane(page, outer).locator(`[name="div-section--${inner.replace(/ /g, '-')}"]`).filter({visible: true}).first();
    await expect(header, `the ${inner} pane inside the ${outer} section`).toBeVisible();
    if (await header.getAttribute('aria-expanded') !== 'true')
      await header.click();
    await expect(header, `the ${inner} pane inside the ${outer} section`).toHaveAttribute('aria-expanded', 'true');
  }, {tier: 'ui', description: 'a click on the header of the inner pane unless it is expanded; done when it says it is'});
