/* The Add New Column dialog of PowerPack (src/dialogs/add-new-column.ts) and the Formula pane it
   also hosts: the CodeMirror formula editor, its hint and error lines, the column list (a Dart grid
   whose rows are the table's columns, `__name` holding the name), the functions list (the platform's
   FunctionsWidget: one row per function, a plus icon revealed on hover, a sort icon), the preview
   grid, and the input history menu. */
import {Page} from '@playwright/test';
import {dataset, element, Given, kind, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, type ElementRef, expect, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const DG: any;
declare const grok: any;

dataset('SPGI', {path: 'System:DemoFiles/chem/SPGI.csv',
  description: '3624 molecules in "Structure" with the numeric, text and date columns of the SPGI demo (Id, Chemist, Species, Whole blood assay 1, Route Admin, Average Mass)'});

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
element('column name input', {selector: '[name="input-Add-New-Column---Name"]',
  description: 'the Name field of the Add New Column dialog; its placeholder follows the formula while it is empty'});
element('column type input', {selector: '[name="input-Add-New-Column---Type"]',
  description: 'the Type choice of the Add New Column dialog: auto (with the type the preview computed), a type, or plain text'});
element('functions sort icon', {selector: '.add-new-column-dialog-root [name="icon-sort-alt"]',
  description: 'the two arrows over the functions list: By name or By relevance'});
element('input history menu', {selector: '[name="input-history"]',
  description: 'the menu the history icon of a dialog opens: one entry per earlier run, the latest first'});
element('completion list', {selector: '.cm-tooltip-autocomplete',
  description: 'the autocomplete popup of the formula editor'});
element('signature tooltip', {selector: '.cm-tooltip-hover',
  description: 'what the formula editor shows over a function name under the pointer'});
element('resize corner of Add New Column dialog', {selector: '.add-new-column-dialog-root .d4-host-bottom-right-resizer',
  description: 'the bottom right corner handle of the dialog'});

kind('completion', {selector: '.cm-tooltip-autocomplete li', match: ['text'],
  description: 'an entry of the formula editor\'s autocomplete popup, by its label ("Round", "HEIGHT")'});
kind('project card', {selector: '.d4-gallery-card.entity-project', match: ['label'],
  labelSelector: '.grok-gallery-grid-item-title',
  description: 'a card of the Projects gallery (Browse > Dashboards), by the title it shows'});
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

async function dragTo(page: Page, from: {x: number; y: number}, target: ElementRef): Promise<void> {
  const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
  if (!box)
    throw new Error(`${target.phrase} has no box to drop onto`);
  const to = {x: box.x + Math.min(40, box.width / 2), y: box.y + Math.min(10, box.height / 2)};
  await page.mouse.move(from.x, from.y);
  await page.mouse.down();
  for (let i = 1; i <= 10; i++)
    await page.mouse.move(from.x + (to.x - from.x) * i / 10, from.y + (to.y - from.y) * i / 10);
  await page.mouse.up();
}

export const clickColumnName = When('user clicks on the {string} column in {widget}', async (page: Page, name: string, target: ElementRef) => {
  const box = await viewers.hitArea(page, target, await columnCell(page, target, name), true);
  const c = viewers.centerOf(box);
  await page.mouse.click(c.x, c.y);
}, {tier: 'ui', description: 'a click on the name of that column in the column list, found by the grid\'s "text of cell N of __name" readings'});

export const dragColumnName = When('user drags the {string} column of {widget} onto {element}', async (page: Page, name: string, source: ElementRef, target: ElementRef) =>
  dragTo(page, viewers.centerOf(await viewers.hitArea(page, source, await columnCell(page, source, name), true)), target),
{tier: 'ui', description: 'a pointer drag from the column\'s name in the column list to the element'});

export const dragAreaOnto = When('user drags the {string} area of {widget} onto {element}', async (page: Page, area: string, source: ElementRef, target: ElementRef) =>
  dragTo(page, viewers.centerOf(await viewers.hitArea(page, source, area, true)), target),
{tier: 'ui', description: 'a pointer drag from a hit area of a widget (a grid column header) to an element outside it'});

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

const rememberedBoxes = new WeakMap<Page, Map<string, {width: number; height: number}>>();

export const rememberSize = When('user remembers the size of {element}', async (page: Page, target: ElementRef) => {
  const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
  if (!box)
    throw new Error(`${target.phrase} has no box`);
  if (!rememberedBoxes.has(page))
    rememberedBoxes.set(page, new Map());
  rememberedBoxes.get(page)!.set(target.phrase, {width: box.width, height: box.height});
}, {tier: 'api', description: 'kept under the phrase, so several elements can be remembered at once'});

export const dragCornerBy = When('user drags {element} by {int} and {int} pixels', async (page: Page, target: ElementRef, dx: number, dy: number) => {
  const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
  if (!box)
    throw new Error(`${target.phrase} has no box to drag`);
  const from = {x: box.x + box.width / 2, y: box.y + box.height / 2};
  await page.mouse.move(from.x, from.y);
  await page.mouse.down();
  for (let i = 1; i <= 10; i++)
    await page.mouse.move(from.x + dx * i / 10, from.y + dy * i / 10);
  await page.mouse.up();
}, {tier: 'ui', description: 'a pointer drag of the element right and down by those amounts (negative: left and up)'});

export const sizeAgainstRemembered = Then('{element} should be {word} than remembered', async (page: Page, target: ElementRef, how: string) => {
  const kept = rememberedBoxes.get(page)?.get(target.phrase);
  if (!kept)
    throw new Error(`no size of ${target.phrase} was remembered — "user remembers the size of ${target.phrase}" comes first`);
  const checks: Record<string, (w: number, h: number) => boolean> = {
    larger: (w, h) => w > kept.width + 20 && h > kept.height + 20,
    smaller: (w, h) => w < kept.width - 20 && h < kept.height - 20,
    wider: (w) => w > kept.width + 20,
    taller: (_w, h) => h > kept.height + 20,
    narrower: (w) => w < kept.width - 20,
  };
  if (!checks[how])
    throw new Error(`"${how}": larger, smaller, wider, taller or narrower`);
  let now = {width: 0, height: 0};
  await expect.poll(async () => {
    const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
    now = box ?? now;
    return checks[how](now.width, now.height);
  }, {message: `${target.phrase} ${how} than remembered`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} is ${Math.round(now.width)}×${Math.round(now.height)}, it was ${Math.round(kept.width)}×${Math.round(kept.height)}`);
  });
}, {description: 'larger/smaller: both width and height by more than 20 px; wider/narrower/taller: that one dimension'});

export const keepsWidth = Then('{element} should keep its remembered width', async (page: Page, target: ElementRef) => {
  const kept = rememberedBoxes.get(page)?.get(target.phrase);
  if (!kept)
    throw new Error(`no size of ${target.phrase} was remembered`);
  const box = await (await locate(page, target)).filter({visible: true}).first().boundingBox();
  expect(Math.round(box?.width ?? -1), `the width of ${target.phrase} (it was ${Math.round(kept.width)})`).toBe(Math.round(kept.width));
}, {description: 'to the pixel, read once after the change that could have moved it'});

export const everyValueEquals = Then('every value of {string} column should equal {string} column plus {float}',
  async (page: Page, column: string, source: string, delta: number) => {
    const bad = await page.evaluate(([c, s, d]) => {
      const df = grok.shell.t;
      const a = df.col(c as string); const b = df.col(s as string);
      if (!a || !b)
        return `no "${a ? s : c}" column; the table has ${df.columns.names().join(', ')}`;
      for (let i = 0; i < df.rowCount; i++) {
        const x = a.get(i); const y = b.get(i);
        if (y === null || b.isNone(i)) {
          if (!a.isNone(i))
            return `row ${i + 1}: ${s} is empty and ${c} is ${x}`;
        }
        else if (!(Math.abs(x - (y + (d as number))) <= 1e-3))
          return `row ${i + 1}: ${c} is ${x}, ${s} is ${y}`;
      }
      return df.rowCount > 0 ? '' : 'the table has no rows';
    }, [column, source, delta] as [string, string, number]);
    expect(bad, `${column} against ${source} + ${delta}`).toBe('');
  }, {description: 'row by row to a thousandth (the columns are 32-bit floats), an empty source cell giving an empty value'});

export const everyValueLog = Then('every value of {string} column should be the decimal log of {string} column minus {float}',
  async (page: Page, column: string, source: string, delta: number) => {
    const bad = await page.evaluate(([c, s, d]) => {
      const df = grok.shell.t;
      const a = df.col(c as string); const b = df.col(s as string);
      if (!a || !b)
        return `no "${a ? s : c}" column; the table has ${df.columns.names().join(', ')}`;
      for (let i = 0; i < df.rowCount; i++) {
        const x = a.get(i); const y = b.get(i);
        if (b.isNone(i)) {
          if (!a.isNone(i))
            return `row ${i + 1}: ${s} is empty and ${c} is ${x}`;
        }
        else if (!(Math.abs(x - (Math.log10(y) - (d as number))) <= 1e-4))
          return `row ${i + 1}: ${c} is ${x}, ${s} is ${y}`;
      }
      return df.rowCount > 0 ? '' : 'the table has no rows';
    }, [column, source, delta] as [string, string, number]);
    expect(bad, `${column} against log10(${source}) - ${delta}`).toBe('');
  }, {description: 'row by row, to four decimals, an empty source cell giving an empty value'});

export const holdsFormula = Then('{element} should hold the formula {string}', async (page: Page, target: ElementRef, text: string) => {
  const loc = (await locate(page, target)).filter({visible: true});
  let seen: string[] = [];
  await expect.poll(async () => {
    seen = await loc.evaluateAll((els) => els.map((e) => e.textContent ?? ''));
    return seen.length === 1 && seen[0] === text;
  }, {message: `the text of ${target.phrase}`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} holds ${seen.length === 1 ? JSON.stringify(seen[0]) : `${seen.length} editors: ${JSON.stringify(seen)}`}, not ${JSON.stringify(text)}`);
  });
}, {description: 'the whole text of exactly one visible editor, character for character ("" for an empty one)'});

export const typeAtCaret = When('user types {string} at the caret', async (page: Page, text: string) => {
  await page.keyboard.type(text);
}, {tier: 'ui', description: 'keys into whatever holds the focus, where its caret is, without selecting what it holds first'});

export const hoverText = When('user hovers over the text {string} in {element}', async (page: Page, text: string, target: ElementRef) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  let box: {x: number; y: number; width: number; height: number} | null = null;
  await expect.poll(async () => (box = await loc.evaluate((root, t) => {
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    for (let n = walker.nextNode(); n; n = walker.nextNode()) {
      const i = n.textContent!.indexOf(t);
      if (i < 0)
        continue;
      const range = document.createRange();
      range.setStart(n, i);
      range.setEnd(n, i + t.length);
      const r = range.getBoundingClientRect();
      return {x: r.x, y: r.y, width: r.width, height: r.height};
    }
    return null;
  }, text)) !== null, {timeout: pollMs(5000), message: `"${text}" in ${target.phrase}`}).toBe(true);
  const b = box!;
  await page.mouse.move(b.x + b.width / 2 - 3, b.y + b.height / 2 + 3);
  await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2, {steps: 2});
}, {tier: 'ui', description: 'the pointer over the middle of the first place the element shows that text'});

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

export const noTablesOpen = Then('no table should be open', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => (grok.shell.tables ?? []).map((t: any) => t.name).join(', ')),
    {message: 'the tables open in the workspace'}).toBe('');
}, {description: 'the workspace holds no table, so what a reopen shows came from the server'});

export const fileInHome = Given('a copy of the {string} file is in the home folder as {string}', async (page: Page, source: string, name: string) => {
  const path: string = await page.evaluate(async ([src, n]) => {
    const project = grok.shell.user.project.name;
    const home = (await grok.dapi.connections.list()).find((c: any) => c.dataSource === 'Files' && c.nqName === `${project}:Home`);
    if (!home)
      throw new Error(`no home folder ${project}:Home on this stand`);
    const target = `${home.nqName}/${n}`;
    await grok.dapi.files.write(target, await grok.dapi.files.readAsBytes(src));
    return target;
  }, [source, name]);
  atFeatureEnd(page, async () => {
    await page.evaluate(async (p) => {
      if (await grok.dapi.files.exists(p))
        await grok.dapi.files.delete(p);
    }, path);
  });
}, {tier: 'api', description: 'the file written into the "My files" share of the current user, and deleted when the feature ends'});

const COMPLETION_INTERACTION_DELAY_MS = 75;

export const acceptCompletion = When('user accepts the highlighted completion with {key}', async (page: Page, key: string) => {
  const list = page.locator('.cm-tooltip-autocomplete').filter({visible: true});
  await expect(list, 'the autocomplete popup').toHaveCount(1);
  // CodeMirror's `interactionDelay` (@codemirror/autocomplete): keys pressed sooner do not go to the popup.
  await page.waitForTimeout(COMPLETION_INTERACTION_DELAY_MS);
  await (await locate(page, {phrase: 'formula editor'})).filter({visible: true}).first().press(key);
}, {tier: 'ui', description: 'the key pressed in the formula editor once the popup is shown and past its 75 ms interaction delay, as a person would; pressed once'});
