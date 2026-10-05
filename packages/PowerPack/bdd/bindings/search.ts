/* The Power Search of the Home page: what a finished search shows. Every search path appends its
   own panel — a category list under an h3 header, a widget with a title — as its answer arrives; the
   results host is aria-busy until the last path has settled, so a claim about what is (or is not)
   there reads after `the search should have finished`. The suggestion menu under the box is a
   Dart popup menu; the highlight the arrow keys move is the `d4-menu-item-hover` class. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

const namesOf = (list: string): string[] => list.split(',').map((n) => n.trim()).filter(Boolean);

function shownCategories(page: Page): Promise<string[]> {
  return page.evaluate(() => [...document.querySelectorAll('.power-pack-search-host .power-pack-search-list-container')]
    .filter((c) => (c as HTMLElement).offsetParent !== null)
    .sort((a, b) => Number((a as HTMLElement).style.order || 0) - Number((b as HTMLElement).style.order || 0))
    .map((c) => (c.querySelector('.power-pack-search-list-header')?.textContent ?? '').trim()));
}

export const searchListsCategories = Then('the search results should list the categories {string}', async (page: Page, list: string) => {
  const wanted = namesOf(list);
  await expect.poll(async () => {
    const shown = await shownCategories(page);
    return wanted.filter((c) => !shown.includes(c)).length === 0 ? 'all there' : `shown: ${shown.join(', ') || 'none'}`;
  }, {message: `the result categories ${list}`}).toBe('all there');
}, {description: 'each named category has a header and a list in the results (others may be there too)'});

export const searchShowsNothing = Then('the search results should show nothing', async (page: Page) => {
  const shown = await page.evaluate(() => {
    const host = (window as any).grok.shell.v.root.querySelector('.power-pack-search-host') as HTMLElement;
    const text = (host.innerText ?? '').trim();
    // every path appends to the host; the only child an empty answer leaves is the empty holder of the category lists
    const others = [...host.children].filter((c) => !(c.classList.contains('power-search-lists-host') && c.childElementCount === 0))
      .map((c) => `<${c.tagName.toLowerCase()} class="${c.className}">${(c.textContent ?? '').trim().slice(0, 60)}`);
    return [...(text ? [`text: ${text.slice(0, 120)}`] : []), ...others];
  });
  expect(shown, 'what the finished search shows').toEqual([]);
}, {description: 'read once the search is over (after "the search should have finished"): no visible text in the results and no element but the empty holder of the category lists — a category, a widget, an evaluated formula or an entity card would each be there'});

export const noSuggestions = Then('no search suggestion should be shown', async (page: Page) => {
  expect(await menuItems(page), 'the suggestions under the search box').toEqual([]);
}, {description: 'read after the search has finished, when the suggestions of that text would have been shown'});

export const categoryListsItem = Then('the {string} category of the search results should list {string}', async (page: Page, category: string, item: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const container = [...document.querySelectorAll('.power-pack-search-host .power-pack-search-list-container')]
      .find((x) => (x.querySelector('.power-pack-search-list-header')?.textContent ?? '').trim() === c);
    return container ? [...container.querySelectorAll('.power-pack-search-list > *')].map((i) => (i.textContent ?? '').trim()) : [];
  }, category), {message: `the items of the ${category} category`}).toContain(item);
});

function menuItems(page: Page): Promise<string[]> {
  return page.evaluate(() => [...document.querySelectorAll('.ui-input-type-ahead .d4-menu-item')]
    .filter((m) => (m as HTMLElement).offsetParent !== null)
    .sort((a, b) => (a as HTMLElement).offsetTop - (b as HTMLElement).offsetTop)
    .map((m) => `${m.classList.contains('d4-menu-item-hover') ? '>' : ''}${(m.querySelector('.d4-menu-item-label')?.textContent ?? m.textContent ?? '').trim()}`));
}

export const suggestionsAre = Then('the search suggestions should be {string}', async (page: Page, list: string) => {
  await expect.poll(() => menuItems(page).then((items) => items.map((i) => i.replace(/^>/, '')).join(' | ')),
    {message: 'the suggestions under the search box, top to bottom'}).toBe(list.split('|').map((n) => n.trim()).join(' | '));
}, {description: 'every suggestion under the box, in the order shown (separated by |)'});

export const suggestionHighlighted = Then('the highlighted search suggestion should be {string}', async (page: Page, text: string) => {
  await expect.poll(() => menuItems(page).then((items) => items.filter((i) => i.startsWith('>')).map((i) => i.slice(1)).join(' | ') || 'none'),
    {message: 'the highlighted suggestion'}).toBe(text);
}, {description: 'exactly one suggestion carries the highlight, and it is this one'});

export const searchFinished = Then('the search should have finished', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const host = (window as any).grok.shell.v?.root?.querySelector('.power-pack-search-host');
    return host ? host.getAttribute('aria-busy') ?? 'no aria-busy' : 'no search results on the current view';
  }), {message: 'aria-busy of the search results of the Home page'}).toBe('false');
}, {description: 'every search path has answered (the results host clears aria-busy); an empty answer has nothing to show, so this is read off the host, not off what is visible'});
