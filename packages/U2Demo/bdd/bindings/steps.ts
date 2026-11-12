import type {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, locate} from '@datagrok-libraries/bdd/runtime';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openSubDemo} from './demo.js';

export const openWorkbench = Given('user opens the MSA workbench', async (page: Page) => {
  await openSubDemo(page, 'MSA workbench', '[data-u2-name="msaWorkbench"]');
}, {tier: 'ui', enters: 'MSA workbench', description: 'needs this package published on the stand'});

/** The pages of a feature that already put a kept choice aside: its first value goes back at feature end. */
const keptAside = new WeakMap<Page, Set<string>>();

export const catalogChoicesCleared = Given('the catalog table choices kept under {string} are cleared now and put back at feature end',
  async (page: Page, key: string) => {
    const id = `u2.catalog.${key}`;
    if (!await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false))
      await loggedIn(page);
    const before: string | null = await page.evaluate((k) => {
      const value = localStorage.getItem(k);
      localStorage.removeItem(k);
      return value;
    }, id);
    const aside = keptAside.get(page) ?? new Set<string>();
    keptAside.set(page, aside);
    if (aside.has(id))
      return;
    aside.add(id);
    atFeatureEnd(page, async () => {
      aside.delete(id);
      await page.evaluate(([k, v]) => v === null ? localStorage.removeItem(k) : localStorage.setItem(k, v), [id, before] as [string, string | null]);
    });
  }, {tier: 'api', description: 'the sort, grouping, columns and folds a CatalogTable keeps in the browser (`persist: <key>`) are cleared before the page opens, and the value found first is put back at feature end'});

/* A u2 table row under the pointer is marked as the platform grid marks its own (css/table.css): the grid's
   mouse-over tint over the table, lighter than a section header, and the grid's stripe at the row's start. The
   colours compared are the painted ones: the background colour with its gradient layers over it. */
export const hoverAsGrid = Then('the row under the pointer in {element} should be marked as the platform grid marks it',
  async (page: Page, target: ElementRef) => {
    const table = (await locate(page, target)).filter({visible: true}).first();
    await expect.poll(() => table.evaluate((root) => {
      const ctx = document.createElement('canvas').getContext('2d', {willReadFrequently: true})!;
      const paint = (layers: string[]): number[] => {
        ctx.clearRect(0, 0, 1, 1);
        for (const layer of layers) {
          ctx.fillStyle = layer;
          ctx.fillRect(0, 0, 1, 1);
        }
        return [...ctx.getImageData(0, 0, 1, 1).data.slice(0, 3)];
      };
      const painted = (el: Element): number[] => {
        const cs = getComputedStyle(el);
        // a solid layer is a gradient between two stops of one colour: painted once
        const layers = [...cs.backgroundImage.matchAll(/linear-gradient\((?:[^,()]*,\s*)?((?:rgba?|color)\([^()]*\))/g)]
          .map((m) => m[1]).reverse();
        return paint([cs.backgroundColor, ...layers]);
      };
      const token = (name: string): string => {
        const probe = document.createElement('div');
        probe.style.color = `var(${name})`;
        root.append(probe);
        const color = getComputedStyle(probe).color;
        probe.remove();
        return color;
      };
      const row = root.querySelector('tr.u2-table-row:hover');
      const header = root.querySelector('thead th');
      if (row === null || header === null)
        return 'no hovered row or no header';
      const tint = painted(row);
      const grid = paint([getComputedStyle(header).backgroundColor, token('--dg-grid-mouse-over-row')]);
      const sum = (c: number[]) => c[0] + c[1] + c[2];
      const section = root.querySelector('.u2-table-section th');
      const stripe = getComputedStyle(row.firstElementChild!).boxShadow;
      const problems = [
        ...tint.join() === grid.join() ? [] : [`the row is painted ${tint}, the grid's tint is ${grid}`],
        ...section === null || sum(tint) > sum(painted(section)) ? [] : [`the row (${tint}) is not lighter than a section header (${painted(section!)})`],
        ...stripe.startsWith(token('--dg-grid-mouse-over-stripe')) && stripe.includes('inset') ? [] : [`the row's first cell has no stripe (${stripe})`],
      ];
      return problems.join('; ');
    }), {message: `the hovered row of ${target.phrase} against the platform grid's and its section headers`}).toBe('');
  }, {tier: 'ui', description: 'the hovered row\'s painted background is the grid\'s mouse-over tint over the table, lighter than a section header, and its first cell carries the grid\'s mouse-over stripe'});
