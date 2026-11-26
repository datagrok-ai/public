/* The Basics how-to clips (features/howto/01-getting-started): the Home page's search box and its results, which the
   Tutorials project's vocabulary lacks (PowerPack's own bdd names them the same way). Every search path of the
   Home search appends its own panel as its answer arrives; the results host is aria-busy until the last one has
   settled, so what the results show is read after `the Home search should have finished`. */
import {Page} from '@playwright/test';
import {Then, element} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

element('home search', {selector: '.power-search-search-everywhere-input',
  description: 'the "Search everywhere" box at the top of the Home page'});
element('home search results', {selector: '.power-pack-search-host',
  description: 'where a Home search shows what it found, a list per kind under its header; aria-busy until every search path has settled'});
export const homeSearchFinished = Then('the Home search should have finished', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const host = grok.shell.v?.root?.querySelector('.power-pack-search-host');
    return host ? host.getAttribute('aria-busy') ?? 'no aria-busy' : 'no search results on the current view';
  }), {message: 'aria-busy of the search results of the Home page'}).toBe('false');
}, {description: 'every search path has answered (the results host clears aria-busy)'});

/* A typed sentence that matches a query's search pattern runs the query and shows its table in the results; the
   widget shows a loader until the query answers, after the results host has already cleared aria-busy. */
export const homeSearchTableLoaded = Then('the table the Home search answered with should have loaded', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const host = document.querySelector('.power-pack-search-host');
    if (!host)
      return 'no search results';
    const loading = [...host.querySelectorAll('.grok-loader')].filter((e) => (e as HTMLElement).offsetParent !== null);
    if (loading.length > 0)
      return 'still loading';
    return host.querySelector('.d4-grid') ? 'loaded' : 'no table in the results';
  }), {timeout: 60000, message: 'the table under the Home search'}).toBe('loaded');
}, {description: 'the query a typed sentence ran has answered and its table is drawn in the results (up to a minute)'});

export const homeSearchCategories = Then('the Home search results should list the kinds {string}', async (page: Page, list: string) => {
  const wanted = list.split(',').map((n) => n.trim()).filter(Boolean);
  await expect.poll(async () => {
    const shown: string[] = await page.evaluate(() => [...document.querySelectorAll('.power-pack-search-host .power-pack-search-list-container')]
      .filter((c) => (c as HTMLElement).offsetParent !== null)
      .map((c) => (c.querySelector('.power-pack-search-list-header')?.textContent ?? '').trim()));
    const missing = wanted.filter((c) => !shown.includes(c));
    return missing.length === 0 ? 'all there' : `missing ${missing.join(', ')}; shown: ${shown.join(', ') || 'none'}`;
  }, {message: `the result kinds ${list}`}).toBe('all there');
}, {description: 'each named kind (Functions, Help, Files…) has its header and list among the results; others may be there too'});
