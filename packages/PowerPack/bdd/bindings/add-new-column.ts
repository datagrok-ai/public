/* The Add New Column dialog of PowerPack (src/dialogs/add-new-column.ts) and the Formula pane it
   also hosts: the CodeMirror formula editor, its hint and error lines, the column list (a Dart grid
   whose rows are the table's columns, `__name` holding the name), the functions list (the platform's
   FunctionsWidget: one row per function, a plus icon revealed on hover, a sort icon), the preview
   grid, and the input history menu. */
import {Page} from '@playwright/test';
import {element, kind, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import {dragAreaOntoElement} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';

declare const DG: any;

element('formula editor', {selector: '.add-new-column-dialog-root .add-new-column-dialog-cm-div .cm-content',
  description: 'the CodeMirror field of the Add New Column dialog'});
element('formula pane editor', {selector: '.add-new-column-widget-cm-div .cm-content',
  description: 'the CodeMirror field of the Formula pane of a calculated column in the context panel'});
element('formula hint', {selector: '.add-new-column-dialog-root .cm-hint-div:not(.cm-errort-div)',
  description: 'the line under the editor: the signature of the function at the caret, else how to start'});
element('formula error', {selector: '.cm-errort-div',
  description: 'the line under the editor that names what is wrong with the formula; empty while it is valid'});
element('column list viewer', {selector: '.add-new-column-columns-grid',
  description: 'the columns of the table, one row each ("text of cell N of __name"), on the left of the functions'});
element('preview grid viewer', {selector: '.ui-addnewcolumn-preview [name="viewer-Grid"]',
  description: 'the first rows of the columns the formula uses and the column it computes'});
element('functions list', {selector: '.add-new-column-dialog-root .grok-actions-browser-table',
  description: 'the functions on the right of the dialog, in the order the sort mode gives'});
element('functions panel', {selector: '.add-new-column-dialog-root .ui-widget-addnewcolumn-functions',
  description: 'the scrolling pane that holds the functions list'});
element('column type input', {selector: '[name="input-Add-New-Column---Type"]',
  description: 'the Type choice of the Add New Column dialog: auto (with the type the preview computed), a type, or plain text'});
element('functions sort icon', {selector: '.add-new-column-dialog-root [name="icon-sort-alt"]',
  description: 'the two arrows over the functions list: By name or By relevance'});
element('completion list', {selector: '.cm-tooltip-autocomplete',
  description: 'the autocomplete popup of the formula editor'});
element('signature tooltip', {selector: '.cm-tooltip-hover',
  description: 'what the formula editor shows over a function name under the pointer'});
element('resize corner of Add New Column dialog', {selector: '.add-new-column-dialog-root .d4-host-bottom-right-resizer',
  description: 'the bottom right corner handle of the dialog'});

kind('completion', {selector: '.cm-tooltip-autocomplete li', match: ['text'],
  description: 'an entry of the formula editor\'s autocomplete popup, by its label ("Round", "HEIGHT")'});
kind('function entry', {selector: '.grok-actions-browser-table tr', match: ['label'],
  labelSelector: '[data-entity-type="Func"] label', parts: {plus: '[name="icon-plus"]', name: '[data-entity-type="Func"]'},
  description: 'a row of the functions list, by the function name; "plus of X function entry" is the icon a hover reveals'});

/** The row of the column list that holds `name`, from the grid's own readings. */
async function columnCell(page: Page, target: ElementRef, name: string): Promise<string> {
  let names: string[] = [];
  await expect.poll(async () => {
    const values = await viewers.onViewer(page, target, (e) =>
      (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {}) as Record<string, unknown>;
    names = Object.entries(values).filter(([k]) => /^text of cell \d+ of __name$/.test(k)).map(([k, v]) => `${k}=${v}`);
    return names.some((n) => n.endsWith(`=${name}`));
  }, {timeout: pollMs(5000), message: `"${name}" in ${target.phrase}`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} shows no "${name}" row; it shows: ${names.map((n) => n.split('=')[1]).join(', ')}`);
  });
  return names.find((n) => n.endsWith(`=${name}`))!.split('=')[0].replace('text of ', '');
}

export const clickColumnName = When('user clicks on the {string} column in {widget}', async (page: Page, name: string, target: ElementRef) => {
  const box = await viewers.hitArea(page, target, await columnCell(page, target, name), true);
  const c = viewers.centerOf(box);
  await page.mouse.click(c.x, c.y);
}, {tier: 'ui', description: 'a click on the name of that column in the column list, found by the grid\'s "text of cell N of __name" readings'});

export const dragColumnName = When('user drags the {string} column of {widget} onto {element}', async (page: Page, name: string, source: ElementRef, target: ElementRef) =>
  dragAreaOntoElement(page, await columnCell(page, source, name), source, target),
{tier: 'ui', description: 'a pointer drag from the column\'s name in the column list to the element'});

/** The functions list as the user sees it: every visible row's function, top to bottom. */
function functionOrder(page: Page): Promise<{name: string; nqName: string}[]> {
  return page.evaluate(() => [...document.querySelectorAll('.add-new-column-dialog-root .grok-actions-browser-table tr [data-entity-type="Func"]')]
    .filter((e) => (e as HTMLElement).offsetParent !== null)
    .map((e) => ({name: (e.getAttribute('data-link') ?? '').replace(/^.*[./]/, ''), nqName: (e.getAttribute('data-link') ?? '').replace('/func/', '').replace('.', ':')})));
}

export const functionsStartWith = Then('the functions list should start with {string}', async (page: Page, names: string) => {
  const expected = names.split(',').map((s) => s.trim());
  let seen: string[] = [];
  await expect.poll(async () => (seen = (await functionOrder(page)).slice(0, expected.length).map((f) => f.name)).join(', '),
    {message: 'the first functions of the list'}).toBe(expected.join(', '));
}, {description: 'the names of the first rows, in order'});

/** What a function takes first, the way the dialog's type sorting reads it: the semantic type when
 * the parameter has one, else its type. */
export const functionsTakeFirst = Then('the first {int} functions of the functions list should take a {string} first',
  async (page: Page, count: number, type: string) => {
    let report = '';
    await expect.poll(async () => {
      const funcs = (await functionOrder(page)).slice(0, count).map((f) => f.nqName);
      const firsts: string[] = await page.evaluate((names) => names.map((n) => {
        const f = DG.Func.find({name: n.split(':').pop()}).find((x: any) => x.nqName.toLowerCase() === n.toLowerCase() || x.name === n.split(':').pop());
        const p = f?.inputs[0];
        return p ? (p.semType || p.propertyType) : 'nothing';
      }), funcs);
      const numeric = ['double', 'int', 'num', 'float', 'qnum', 'bigint', 'number'];
      const matches = (t: string) => type === 'number' ? numeric.includes(t) : t === type;
      report = funcs.map((f, i) => `${f}(${firsts[i]})`).join(', ');
      return funcs.length === count && firsts.every(matches);
    }, {message: `the first ${count} functions to take a ${type}`}).toBe(true).catch(() => {
      throw new Error(`the first ${count} functions of the list are ${report}`);
    });
  }, {description: 'the first parameter of each of the top rows: its semantic type if it has one, else its type ("number" for any numeric type)'});

export const functionsByName = Then('the functions list should be sorted by name', async (page: Page) => {
  let names: string[] = [];
  await expect.poll(async () => {
    names = (await functionOrder(page)).map((f) => f.name);
    const sorted = [...names].sort();
    return names.length > 20 && names.join(',') === sorted.join(',');
  }, {message: 'the functions list in alphabetical order'}).toBe(true).catch(() => {
    const i = names.findIndex((n, k) => k > 0 && names[k - 1] > n);
    throw new Error(`the functions list is not alphabetical (${names.length} rows): ${i < 0 ? names.slice(0, 5).join(', ') : `"${names[i - 1]}" comes before "${names[i]}"`}`);
  });
}, {description: 'every visible row by its function name, in character order (capitals first, as the platform sorts); more than 20 rows, so a filtered list does not pass'});

const rememberedOrders = new WeakMap<Page, string>();

export const rememberFunctions = When('user remembers the order of the functions list', async (page: Page) => {
  const order = (await functionOrder(page)).map((f) => f.nqName).join(',');
  expect(order, 'the functions list').not.toBe('');
  rememberedOrders.set(page, order);
}, {tier: 'api'});

export const functionsAsRemembered = Then('the functions list should be in the remembered order', async (page: Page) => {
  const kept = rememberedOrders.get(page);
  if (kept === undefined)
    throw new Error('no order was remembered — "user remembers the order of the functions list" comes first');
  expect((await functionOrder(page)).map((f) => f.nqName).join(','), 'the functions list against the remembered order').toBe(kept);
}, {description: 'every row in the same place; read once, after the clicks that could have moved them'});

export const functionsNotAsRemembered = Then('the functions list should not be in the remembered order', async (page: Page) => {
  const kept = rememberedOrders.get(page);
  if (kept === undefined)
    throw new Error('no order was remembered — "user remembers the order of the functions list" comes first');
  await expect.poll(async () => (await functionOrder(page)).map((f) => f.nqName).join(','),
    {message: 'the functions list against the remembered order'}).not.toBe(kept);
});

export const previewComputes = Then('the preview grid should show {string} computed as numbers', async (page: Page, column: string) => {
  let seen = '';
  await expect.poll(async () => {
    seen = await page.evaluate((c) => {
      const root = document.querySelector('.ui-addnewcolumn-preview [name="viewer-Grid"]');
      const grid = root ? DG.Widget.find(root as HTMLElement) : null;
      const col = grid?.dataFrame?.col(c);
      if (!col)
        return `no "${c}" column; it has ${grid?.dataFrame?.columns.names().join(', ')}`;
      const values = Array.from({length: col.length}, (_, i) => col.get(i));
      const numbers = values.filter((v) => typeof v === 'number' && Number.isFinite(v));
      return `${numbers.length} of ${values.length} (${col.type})`;
    }, column);
    const m = /^(\d+) of (\d+)/.exec(seen);
    return m !== null && Number(m[2]) > 0 && m[1] === m[2];
  }, {timeout: pollMs(15000), message: `the numbers of "${column}" in the preview`}).toBe(true).catch(() => {
    throw new Error(`the preview grid shows ${seen} finite numbers in "${column}"`);
  });
}, {description: 'the column the formula computes, over every preview row, all finite numbers'});

export const previewLacks = Then('the preview grid should not show {string}', async (page: Page, column: string) => {
  let names = '';
  await expect.poll(async () => {
    names = await page.evaluate(() => {
      const root = document.querySelector('.add-new-column-dialog-root .ui-addnewcolumn-preview [name="viewer-Grid"]');
      const grid = root ? DG.Widget.find(root as HTMLElement) : null;
      return grid?.dataFrame ? grid.dataFrame.columns.names().join('|') : '<no preview>';
    });
    return names !== '<no preview>' && !names.split('|').includes(column);
  }, {message: `"${column}" in the preview`}).toBe(true).catch(() => {
    throw new Error(`the preview grid shows the columns ${names}`);
  });
}, {description: 'the preview has been rebuilt without that column — what a preview left over from an earlier formula still shows'});

export const previewAbs = Then('the preview grid should show {string} as the absolute value of {string}', async (page: Page, column: string, source: string) => {
  let seen = '';
  await expect.poll(async () => {
    seen = await page.evaluate(([c, s]) => {
      const root = document.querySelector('.add-new-column-dialog-root .ui-addnewcolumn-preview [name="viewer-Grid"]');
      const df = root ? DG.Widget.find(root as HTMLElement)?.dataFrame : null;
      const a = df?.col(c); const b = df?.col(s);
      if (!a || !b)
        return `no "${a ? s : c}" column; it has ${df?.columns.names().join(', ')}`;
      for (let i = 0; i < df.rowCount; i++) {
        if (b.isNone(i))
          continue;
        if (!(Math.abs(a.get(i) - Math.abs(b.get(i))) <= 1e-4))
          return `row ${i + 1}: ${c} is ${a.get(i)}, ${s} is ${b.get(i)}`;
      }
      return df.rowCount > 0 ? '' : 'no rows';
    }, [column, source]);
    return seen;
  }, {timeout: pollMs(15000), message: `${column} against |${source}| in the preview`}).toBe('');
}, {description: 'every preview row, to four decimals'});

/** The runs of highlighted text in the formula editor, in document order: CodeMirror splits one
 * mark into several spans where another decoration (the bracket match at the caret) cuts it, so
 * adjacent highlighted text is joined back into the reference it is. */
function highlightedRuns(page: Page, target: ElementRef): Promise<{text: string; colors: string[]}[]> {
  return locate(page, target).then((loc) => loc.filter({visible: true}).first().evaluate((root) => {
    const runs: {text: string; colors: string[]}[] = [];
    let open = false;
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    for (let n = walker.nextNode(); n; n = walker.nextNode()) {
      const mark = n.parentElement?.closest('.cm-column-name');
      if (!mark) {
        open = false;
        continue;
      }
      if (!open)
        runs.push({text: '', colors: []});
      const run = runs[runs.length - 1];
      run.text += n.textContent;
      run.colors.push(getComputedStyle(n.parentElement!).color);
      open = true;
    }
    return runs;
  }));
}

export const highlightsExactly = Then('{element} should highlight the column references {string}', async (page: Page, target: ElementRef, list: string) => {
  const expected = list === '' ? [] : list.split(',').map((s) => s.trim());
  let seen: string[] = [];
  await expect.poll(async () => (seen = (await highlightedRuns(page, target)).map((r) => r.text)).join(' | '),
    {message: `the highlighted column references of ${target.phrase}`}).toBe(expected.join(' | '));
}, {description: 'every run of highlighted text, in order; "" for none — a reference cut into spans by the bracket match counts once'});

export const highlightsInColor = Then('every column reference of {element} should be drawn in the color of {string}', async (page: Page, target: ElementRef, cssVar: string) => {
  const want = await page.evaluate((v) => {
    const s = document.createElement('span');
    s.style.color = `var(${v})`;
    document.body.append(s);
    const c = getComputedStyle(s).color;
    s.remove();
    return c;
  }, cssVar);
  const runs = await highlightedRuns(page, target);
  expect(runs.length, `the highlighted column references of ${target.phrase}`).toBeGreaterThan(0);
  const off = runs.filter((r) => r.colors.some((c) => c !== want)).map((r) => `${r.text} (${[...new Set(r.colors)].join(', ')})`);
  expect(off, `references not drawn in ${cssVar} = ${want}`).toEqual([]);
}, {description: 'the computed color of every span of every highlighted reference against the design token, resolved in the page'});

export const highlightStandsOut = Then('every column reference of {element} should differ in color from the plain text of its line', async (page: Page, target: ElementRef) => {
  const report: string[] = await (await locate(page, target)).filter({visible: true}).first().evaluate((root) => {
    const out: string[] = [];
    for (const mark of Array.from(root.querySelectorAll('.cm-column-name'))) {
      const line = mark.closest('.cm-line');
      const plain = line ? getComputedStyle(line).color : '';
      const color = getComputedStyle(mark).color;
      if (!line || color === plain)
        out.push(`${mark.textContent} is ${color}, the plain text of its line is ${plain || 'not found'}`);
    }
    return root.querySelector('.cm-column-name') ? out : ['no highlighted column reference'];
  });
  expect(report, 'column references drawn like the plain text around them').toEqual([]);
}, {description: 'each highlighted span against the color its line gives unhighlighted text'});

const COMPLETION_INTERACTION_DELAY_MS = 75;

export const acceptCompletion = When('user accepts the highlighted completion with {key}', async (page: Page, key: string) => {
  const list = page.locator('.cm-tooltip-autocomplete').filter({visible: true});
  await expect(list, 'the autocomplete popup').toHaveCount(1);
  // CodeMirror's `interactionDelay` (@codemirror/autocomplete): keys pressed sooner do not go to the popup.
  await page.waitForTimeout(COMPLETION_INTERACTION_DELAY_MS);
  await (await locate(page, {phrase: 'formula editor'})).filter({visible: true}).first().press(key);
}, {tier: 'ui', description: 'the key pressed in the formula editor once the popup is shown and past its 75 ms interaction delay, as a person would; pressed once'});
