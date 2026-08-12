/**
 * Shared helpers for the FormsViewer section (PowerGrid multi-form `Forms` viewer, viewer.type
 * `FormsViewer` — NOT the core D4 `Form` viewer; reference references/viewers/formsviewer.md). Field
 * elements carry `[column="<Name>"]`, never `name=`; ordinary cards live under `#vlist`, pinned cards
 * under `.d4-multi-form-pinned-forms`.
 *
 * SECTION VOCABULARY (authoritative — the term means the SAME thing here and in the prose):
 *   - "ordinary cards" = the WHOLE set under `#vlist`, INCLUDING the two leading positions (current-row
 *     card + the always-built mouse-over card). The `ORDINARY` export selects this whole set.
 *   - "selected-row cards" / "cards after the leading offset" = `ORDINARY` sliced past the leading cards
 *     (`fieldValuesByPosition(page, col).slice(offset)`). This is what the p0 claim "the pinned row is
 *     absent from the ordinary cards" ranges over; it holds only because setup keeps the pin out of both leading slots.
 *   - The current-row card is addressed by the `CURRENT` selector — its green indicator, NOT position.
 */

import {Page, expect} from '@playwright/test';

// Section selectors.

/** Host container of the Forms viewer (the viewer's own `.d4-multi-form` root carries no name). */
export const HOST = '[name="viewer-Forms"]';

/**
 * The "ordinary cards" set: the WHOLE `#vlist` set including the two leading positions. `:not(.temp)`
 * is load-bearing — renderHeader() appends a throwaway `.temp` row-0 form to document.body to measure
 * header geometry; host + `#vlist` scoping already excludes it, but the explicit `:not(.temp)` guards
 * a build that attaches the temp form inside the host.
 */
export const ORDINARY = `${HOST} #vlist .d4-multi-form-form:not(.temp)`;

/** The current-row card, identified by its green indicator rather than by position. */
export const CURRENT = `${ORDINARY}:has(.d4-multi-form-form-indicator-current-row)`;

/** The pinned-cards pane (a sibling of `#vlist`, outside the scroller; `display:none` when empty). */
export const PINNED_PANE = `${HOST} .d4-multi-form-pinned-forms`;

/** Pinned cards. */
export const PINNED = `${PINNED_PANE} .d4-multi-form-form`;

// Card / label reads.

/**
 * Names of the drawn header labels IN DRAW ORDER, read independently of the `fieldsColumnNames`
 * property — catches a viewer that restored the property but failed to DRAW a field. Reads the inner
 * `div[name^="div-"]` text so the sort-indicator arrow (a sibling) is excluded.
 */
export async function drawnLabelNames(page: Page): Promise<string[]> {
  return page.evaluate((host) => Array.from(document.querySelectorAll(
    `${host} .d4-multi-form-header .d4-multi-form-column-name div[name^="div-"]`))
    .map((e) => (e.textContent ?? '').trim()).filter((n) => n.length > 0), HOST);
}

/**
 * Balloons of ANY severity — an error-only check would miss the `.d4-balloon.warning` a cap notice uses.
 */
export async function balloonCount(page: Page): Promise<number> {
  return page.locator('.d4-balloon').count();
}

/**
 * Value of a card's `[column]` field by position. `cardSel` defaults to ORDINARY; pass PINNED for a
 * pinned card. Returns `null` when the card or field element is absent (an empty leading card has none).
 *
 * WARNING — reads the INPUT `.value` ONLY. Renderer-backed fields (molecule / curve) are CANVAS with
 * no `.value` → `null`. The field-lifecycle spec keeps a LOCAL reader for that canvas case; do NOT
 * replace it with this helper.
 */
export async function cardFieldValue(
  page: Page, cardIndex: number, column: string, cardSel: string = ORDINARY,
): Promise<string | null> {
  return page.evaluate(({sel, i, col}) => {
    const card = document.querySelectorAll(sel)[i] as HTMLElement | undefined;
    const el = card?.querySelector(`[column="${col}"]`) as HTMLInputElement | null;
    return el ? el.value : null;
  }, {sel: cardSel, i: cardIndex, col: column});
}

/**
 * Values of a `[column]` field for EVERY card BY POSITION. Empty leading cards contribute `null`, so
 * the array stays position-aligned (length = card count). Callers slice by the leading offset FIRST,
 * then drop the nulls; filtering before the slice collapses the array and misaligns the offset.
 * `cardSel` defaults to ORDINARY; pass PINNED for the pinned pane.
 */
export async function fieldValuesByPosition(
  page: Page, column: string, cardSel: string = ORDINARY,
): Promise<(string | null)[]> {
  return page.evaluate(({sel, col}) => Array.from(document.querySelectorAll(sel))
    .map((c) => ((c as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement)?.value ?? null),
  {sel: cardSel, col: column});
}

/**
 * Index of the FIRST card whose `[column]` field value equals `value` — the by-VALUE lookup, so a
 * card is found across a sort/filter re-render that shifted its index. Returns -1 on no match.
 * `cardSel` defaults to ORDINARY; pass PINNED for pinned cards.
 *
 * Resolves duplicates to the FIRST match. WARNING — forms-spec.ts Step 5b keeps an inline loop that
 * does NOT break (resolves to the LAST match); they disagree only on duplicate values. Converting
 * Step 5b to this helper flips the tie-break last→first — do it consciously.
 */
export async function cardIndexByValue(
  page: Page, column: string, value: string, cardSel: string = ORDINARY,
): Promise<number> {
  return page.evaluate(({sel, col, val}) => {
    const cards = Array.from(document.querySelectorAll(sel));
    for (let i = 0; i < cards.length; i++) {
      const el = (cards[i] as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement | null;
      if (el && el.value === val) return i;
    }
    return -1;
  }, {sel: cardSel, col: column, val: value});
}

// Readiness barriers.

/**
 * Count `console.error`s emitted while `fn` runs — the no-error floor for the section's
 * column-lifecycle, canvas-paint and empty-field steps. 400 ms settle window: the viewer swallows a
 * semType-null build error to console.error inside a 50 ms-debounced re-render, so a late error can
 * land a few hundred ms after the action; 400 ms catches it, 300 ms could close first.
 *
 * Returns the COUNT. Pass an optional `texts` array to also receive each caught error's verbatim text
 * (in emission order), so a caller failing on a non-zero count can name the offending text.
 */
export async function withConsoleErrorCount(
  page: Page, fn: () => Promise<void>, settleMs = 400, texts?: string[],
): Promise<number> {
  let count = 0;
  const handler = (msg: {type(): string; text(): string}) => {
    if (msg.type() === 'error') { count++; if (texts) texts.push(msg.text()); }
  };
  page.on('console', handler);
  try {
    await fn();
    await page.waitForTimeout(settleMs);
  } finally {
    page.off('console', handler);
  }
  return count;
}

/**
 * Wait until the card ORDER settles: two consecutive value-sequence reads agree. The card COUNT
 * stabilises before the ORDER does, so a count-only gate can read an anchor off a pre-reorder card
 * while a later right-click lands on the post-reorder card at the same index. `column` is the order
 * key (USUBJID by default); `cardSel` defaults to ORDINARY.
 */
export async function waitForOrderStable(
  page: Page, opts: {column?: string; cardSel?: string; timeoutMs?: number} = {},
): Promise<void> {
  const column = opts.column ?? 'USUBJID';
  const cardSel = opts.cardSel ?? ORDINARY;
  let prev: string | null = null;
  await expect.poll(async () => {
    const cur = await page.evaluate(({sel, col}) => JSON.stringify(
      Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)), {sel: cardSel, col: column});
    const stable = prev !== null && cur === prev;
    prev = cur;
    return stable;
  }, {timeout: opts.timeoutMs ?? 20_000, intervals: [250, 250, 250, 300, 500]}).toBe(true);
}
