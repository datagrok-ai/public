/* The steps only the Forms viewer needs. Everything else its features say is read from the areas
   and readings it reports (`card <n>`, `current card`, `mouse-over card`, `pinned card <n>`,
   `field <COL> of <label>`, `label <COL>`, `remove <COL>`, `sort indicator <COL>`; `cards`,
   `records shown`, `fields`, `header labels`, `record of card <n>`, `card kind of card <n>`,
   `pinned values`, `pinned by`, `pinned pane shown`, `<COL> of <label>` and its width / height /
   background / align / font) — see `public/libraries/utils/src/viewers/CLAUDE.md`.

   What is here is the one thing the readings cannot say in a sentence: WHICH table rows the record
   cards stand for, in the order they stand in. That is the claim the four specs this replaces made
   by re-deriving `getSortedOrder(...)` inside the page and comparing it with a scrape of the card
   DOM — the test computing the answer it was checking. `record of card <n>` gives the row, and
   `card kind of card <n>` separates the record cards from the two leading ones and from the pinned
   pane, so an expected row order is all a scenario has to write down.

   `every record card of {widget} should show {string} in {string}` exists because the library's
   promoted `every tile of {widget} should show {string} in {string}` reads keys shaped
   `<COL> of row <n>` (the tile viewer's), while the Forms viewer labels its cards `card <n>` /
   `current card` / `pinned card <n>`. The two card viewers do not share the key shape; making them
   share it — in the library step or in one of the two statuses — is the fix, and this is the
   stand-in until then. */
import {expect, Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

interface Card {label: string; kind: string; record: number}

async function cardsOf(page: Page, target: ElementRef): Promise<Card[]> {
  await viewers.installViewerRuntime(page);
  const loc = await viewers.viewerLocator(page, target);
  return loc.evaluate((el) => {
    const values: Record<string, unknown> = (window as any).__bdd.viewerOf(el).getWidgetStatus()?.values ?? {};
    const out: Card[] = [];
    for (const key of Object.keys(values)) {
      const m = /^card kind of (.+)$/.exec(key);
      if (m === null)
        continue;
      const record = values[`record of ${m[1]}`];
      out.push({label: m[1], kind: String(values[key]), record: record === '' || record == null ? -1 : Number(record)});
    }
    const index = (label: string): number => {
      const n = /(\d+)$/.exec(label);
      return n === null ? 0 : Number(n[1]);
    };
    return out.sort((a, b) => index(a.label) - index(b.label));
  });
}

const rowsOf = (cards: Card[], kind: string): number[] =>
  cards.filter((c) => c.kind === kind && c.record > 0).map((c) => c.record);

const wanted = (list: string): number[] =>
  list.trim() === '' ? [] : list.split(/\s*,\s*/).map(Number);

function expectRows(kind: string, phrase: string): (page: Page, target: ElementRef, list: string) => Promise<void> {
  return async (page, target, list) => {
    let now: number[] = [];
    await expect.poll(async () => {
      now = rowsOf(await cardsOf(page, target), kind);
      return now.join(', ');
    }, {timeout: 5000, message: `the ${phrase} of ${target.phrase}`}).toBe(wanted(list).join(', '));
  };
}

export const recordCardRows = Then('the record cards of {widget} should show rows {string}',
  expectRows('record', 'record cards'),
{description: 'the table rows the record cards stand for, in card order, 1-based — a sort claim is an expected order'});

export const pinnedCardRows = Then('the pinned cards of {widget} should show rows {string}',
  expectRows('pinned', 'pinned cards'),
{description: 'the same for the pinned pane, which keeps its own cards'});

export const recordCardsAreSelection = Then('the record cards of {widget} should be exactly the selected rows that pass the filter',
  async (page: Page, target: ElementRef) => {
    let detail = '';
    const expected: number[] = await page.evaluate(() => {
      const df = (grok as any).shell.t;
      const rows: number[] = [];
      for (let i = 0; i < df.rowCount; i++)
        if (df.selection.get(i) && df.filter.get(i))
          rows.push(i + 1);
      return rows;
    });
    if (expected.length === 0)
      throw new Error('no selected row passes the filter, so the claim would hold of an empty card set');
    await expect.poll(async () => rowsOf(await cardsOf(page, target), 'record').sort((a, b) => a - b).join(', '),
      {timeout: 5000, message: `the record cards of ${target.phrase} against the selection and the filter`})
      .toBe(expected.join(', '));
  }, {description: 'the rule the record cards follow, read off the table rather than written down — the virtual view must have laid every one of them out, so keep the selection small enough to fit'});

export const everyRecordCardShows = Then('every record card of {widget} should show {string} in {string}',
  async (page: Page, target: ElementRef, value: string, column: string) => {
    await expect.poll(async () => {
      const cards = (await cardsOf(page, target)).filter((c) => c.kind === 'record' && c.record > 0);
      if (cards.length === 0)
        return 'no record card is laid out';
      await viewers.installViewerRuntime(page);
      const loc = await viewers.viewerLocator(page, target);
      const values: Record<string, unknown> = await loc.evaluate((el) =>
        (window as any).__bdd.viewerOf(el).getWidgetStatus()?.values ?? {});
      return cards.filter((c) => String(values[`${column} of ${c.label}`]) !== value)
        .map((c) => `row ${c.record} shows "${values[`${column} of ${c.label}`]}"`).join(', ');
    }, {timeout: 5000,
      message: `every record card of ${target.phrase} shows "${value}" in "${column}"`}).toBe('');
  }, {description: 'the card\'s own text for a column on every record card — the "all cards are Critical" claim'});
