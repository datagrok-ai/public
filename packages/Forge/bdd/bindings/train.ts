/* How a Forge feature gets its data and its fixtures, the shown training, the train form's groups, readings noted and
   compared, and the gestures on a running training. */
import type {Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {el, type ElementRef, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import {catalogRowOf, modelArea} from './catalog-grid.js';
import {expectCreatedModel, guardIdentity, listOf, sweepAtFeatureEnd} from './grok-s.js';
import {slidePointer} from './pointer.js';
import {contextPanelTitle, readingText, type ResultsState, resultsState, settled, shownRun, shownSeed} from
  './results-state.js';

declare const grok: any;
declare const DG: any;

/* "Ready" alone is no proof of a training: a method waiting for Train is ready too, and the metrics grid is updated in
   place. A new result is a ready view whose shown run (data-forge-run-id) differs from the one remembered. */
const shownRuns = new WeakMap<Page, string | null>();

export const rememberShownTraining = When('user remembers the shown training', async (page: Page) => {
  await settled(page);
  shownRuns.set(page, await shownRun(page));
}, {tier: 'api', description: 'once the view has settled, the run whose result it shows (none while it shows a hint)'});

async function expectNewResult(page: Page, ms: number): Promise<void> {
  const before = shownRuns.get(page) ?? null;
  let last: ResultsState = {busy: null, invalid: null, run: null};
  await expect.poll(async () => {
    last = await resultsState(page);
    return last.busy === 'false' && last.invalid !== 'true' && last.run !== null && last.run !== before;
  }, {timeout: pollMs(ms), message: `a new training's result (the shown run was ${before ?? 'none'})`}).toBe(true);
  shownRuns.set(page, last.run);
}

export const newTrainingShown = Then('a new training\'s result should be shown', (page: Page) => expectNewResult(page, 15000),
  {description: 'the view is ready and shows the result of a run other than the one remembered; that run is remembered now'});

export const newTrainingShownWithin = Then('a new training\'s result should be shown within {int} seconds',
  (page: Page, seconds: number) => expectNewResult(page, seconds * 1000),
  {description: 'the same, with a budget of its own (SVM on thousands of rows)'});

export const sameTrainingShown = Then('the same training\'s result should be shown', async (page: Page) => {
  await settled(page);
  expect(await shownRun(page), 'the shown run, against the one remembered').toBe(shownRuns.get(page) ?? null);
}, {description: 'once settled, the shown run is still the one remembered: nothing new was trained and shown'});

export const clipboardHoldsSeed = Then('the clipboard should hold the shown seed', async (page: Page) => {
  const seed = await shownSeed(page);
  await expect.poll(() => gestures.readClipboard(page), {message: `the clipboard, the Performance panel shows ${seed}`})
    .toBe(seed);
}, {description: 'the clipboard text against the seed in the Performance panel'});

export const readingIsShownSeed = Then('the {string} reading of {widget} should be the shown seed',
  async (page: Page, name: string, target: ElementRef) => {
    const seed = await shownSeed(page);
    await expect.poll(() => readingText(page, target, name), {message: `"${name}" of ${target.phrase}, the seed shown ${seed}`})
      .toBe(seed);
  }, {description: 'a grid cell\'s text against the seed in the Performance panel'});

/* LIB-PROMOTE: the library compares one reading name across two viewers; the Runs grid names its columns after the
   metric ("Accuracy (validation)") and the Performance grid after the set ("Validation"). */
export const readingsEqual = Then('the {string} reading of {widget} should equal the {string} reading of {widget}',
  async (page: Page, name: string, target: ElementRef, otherName: string, other: ElementRef) => {
    // Both are read in every attempt: a grid drawing its columns reports a cell of one before the other.
    await expect.poll(async () => {
      const [value, expected] = await Promise.all([viewers.readingOf(page, target, name),
        viewers.readingOf(page, other, otherName)]);
      const isRead = !(value instanceof viewers.MissingReading) && !(expected instanceof viewers.MissingReading);
      return isRead && String(value) === String(expected) ? 'equal' : `${String(value)} against ${String(expected)}`;
    }, {message: `"${name}" of ${target.phrase} against "${otherName}" of ${other.phrase}`}).toBe('equal');
  }, {description: 'the two readings as text, both present, read together until they agree'});

export const readingHigher = Then('the {string} reading of {widget} should be higher than the {string} reading of {widget}',
  async (page: Page, name: string, target: ElementRef, otherName: string, other: ElementRef) => {
    const [a, b] = [await readingText(page, target, name), await readingText(page, other, otherName)];
    expect(Number(a) > Number(b), `"${name}" of ${target.phrase} (${a}) against "${otherName}" of ${other.phrase} (${b})`)
      .toBe(true);
  }, {description: 'the two readings as numbers, read once after the view is ready'});

export const catalogLists = Then('{widget} should list the Forge model {string}', async (page: Page, target: ElementRef,
  name: string) => {
  await catalogRowOf(page, target, name);
}, {description: 'a drawn cell of the name column shows the model\'s name'});

export const clickCatalogRow = When('user clicks on the Forge model {string} in {widget}', async (page: Page, name: string,
  target: ElementRef) => {
  const c = await modelArea(page, target, name, (row) => `cell ${row} of name`);
  await page.mouse.click(c.x, c.y);
}, {tier: 'ui', description: 'clicks the name cell of the row that shows the model'});

/** What the context panel shows of the Forge model [name]: 'shown', or what it is instead. */
async function panelOfModel(page: Page, name: string): Promise<string> {
  const [current, title] = await Promise.all([page.evaluate((n) => {
    const o: any = grok.shell.o;
    if (!(o instanceof DG.DomainRow) || o.typeName !== 'forge.model')
      return `the current object is ${String(o)}`;
    return o.values.name === n ? 'the model' : `the current object is the Forge model "${String(o.values.name)}"`;
  }, name), contextPanelTitle(page)]);
  // The title, not the panel's text: a table's panel names its models too (the pane Forge models).
  return current !== 'the model' ? current : title === name ? 'shown' : `the context panel's title is "${title}"`;
}

/* LIB-PROMOTE: "the context panel should show {string}" matches String(grok.shell.o), which for a domain row is not
   its name; checking the panel's title (or a `title` part of `context panel`) would cover a domain row. */
export const panelShowsModel = Then('the context panel should show the Forge model {string}', async (page: Page,
  name: string) => {
  await expect.poll(() => panelOfModel(page, name), {message: `the context panel and the Forge model "${name}"`})
    .toBe('shown');
}, {description: 'the current object is the forge.model row of that name, and the context panel\'s title is the name'});

/* The context panel is built again only for a model it does not show already (`showModel` in forge-app.ts; an
   application of the shown model is followed in place): a claim that needs a model's panel built anew (Applicable to
   after the open tables changed, a model the panel may already show) has another row shown first. */
export const reopenModel = When('user reopens the Forge model {string} in {widget}', async (page: Page, name: string,
  target: ElementRef) => {
  const other = await viewers.onViewer(page, target, (el, n) => {
    const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
    const names = Object.entries(values).filter(([k, v]) => /^text of cell \d+ of name$/.test(k) && v !== n)
      .map(([, v]) => String(v));
    return names.find((s) => s.startsWith('forge-test-')) ?? names[0] ?? null;
  }, name);
  if (other === null)
    throw new Error(`${target.phrase} draws no model other than "${name}"`);
  for (const model of [other, name]) {
    const c = await modelArea(page, target, model, (row) => `cell ${row} of name`);
    await page.mouse.click(c.x, c.y);
    await expect.poll(() => panelOfModel(page, model), {message: `the context panel and the Forge model "${model}"`})
      .toBe('shown');
  }
}, {tier: 'ui', description: 'clicks the name cell of another drawn model (a forge-test one first) and, once the context ' +
  'panel shows it, the model\'s; done when the panel shows the model'});

/* LIB-PROMOTE: the negative of "the context panel should show the Forge model". */
export const panelShowsNotModel = Then('the context panel should not show the Forge model {string}', async (page: Page,
  name: string) => {
  await expect.poll(async () => (await panelOfModel(page, name)).startsWith('the current object is'),
    {message: `the context panel left the Forge model "${name}"`}).toBe(true);
}, {description: 'polled: the current object is not the forge.model row of that name (another object, or nothing); ' +
  'read after the change that removed it (a delete and the catalog\'s reload)'});

/* LIB-PROMOTE: "{string} link in notification" takes the first balloon with the link; after two saves both balloons
   have one, and the claim is which model each opens. */
export const clickBalloonLink = When('user clicks on {string} link in the balloon of the Forge model {string}',
  async (page: Page, link: string, name: string) => {
    const balloon = page.locator('.d4-balloon', {hasText: `"${name}"`}).filter({visible: true}).first();
    await expect(balloon, `the balloon of "${name}"`).toBeVisible();
    await balloon.getByText(link, {exact: true}).first().click();
  }, {tier: 'ui', description: 'the link of the balloon whose text names the model in quotes'});

/** The shown panels of the results by row (the panels with one top), each row keyed by its titles in order: the
 * [width, height] of each panel by its title, without a value after it ("ROC Curve (AUC = 0.906)" is "ROC Curve"). */
type PanelRows = Record<string, Record<string, [number, number]>>;
const rememberedRows = new WeakMap<Page, PanelRows>();

function panelRows(page: Page): Promise<PanelRows> {
  return page.locator('.forge-results-panel').evaluateAll((all) => {
    const tops = new Map<number, [string, [number, number]][]>();
    for (const p of all) {
      const r = p.getBoundingClientRect();
      const title = p.querySelector(':scope > .forge-results-title')?.textContent?.trim().replace(/ \(.*\)$/, '');
      if (title && r.width > 0) {
        const top = Math.round(r.top);
        tops.set(top, [...tops.get(top) ?? [], [title, [Math.round(r.width), Math.round(r.height)]]]);
      }
    }
    return Object.fromEntries([...tops.values()].map((row) => [row.map(([t]) => t).join(', '), Object.fromEntries(row)]));
  });
}

export const rememberPanelSizes = When('user remembers the sizes of the results panels', async (page: Page) => {
  rememberedRows.set(page, await panelRows(page));
}, {tier: 'api', description: 'the width and height of every shown panel of the results, by its row and its title'});

/* A row of the results is laid out again when a panel of it comes or goes (Insights & Tips after a cutoff): the rows
   that keep their panels are claimed, and no remembered panel may go. */
export const panelSizesKept = Then('the results panels should keep their remembered sizes', async (page: Page) => {
  const kept = rememberedRows.get(page);
  if (kept === undefined)
    throw new Error('no panel sizes were remembered');
  const now = await panelRows(page);
  const shown = new Set(Object.values(now).flatMap((row) => Object.keys(row)));
  const gone = Object.values(kept).flatMap((row) => Object.keys(row)).filter((t) => !shown.has(t));
  const rows = Object.keys(kept).filter((k) => k in now);
  if (gone.length > 0 || rows.length === 0) {
    throw new Error(`the rows of the results were ${Object.keys(kept).join(' | ')}, now ${Object.keys(now).join(' | ')}` +
      (gone.length > 0 ? `: ${gone.join(', ')} gone` : ': no row kept its panels'));
  }
  expect(Object.fromEntries(rows.map((k) => [k, now[k]])), 'the rows of the results that kept their panels: [width, ' +
    'height] by title').toEqual(Object.fromEntries(rows.map((k) => [k, kept[k]])));
}, {description: 'every row of panels (the panels with one top) whose panels are the ones remembered, each panel to ' +
  'the pixel; no remembered panel gone, and at least one row kept; read once after the change that could have moved it'});

/** The captions of the shown groups of the train form, in order, or of the shown inputs of the group [group]. */
function shownCaptions(page: Page, group: string | null): Promise<string> {
  return page.locator('.forge-train-form').evaluateAll((forms, g) => {
    const form = forms.find((f) => f.getClientRects().length > 0);
    const isShown = (e: Element) => e.getClientRects().length > 0;
    const caption = (e: Element | null) => (e?.textContent ?? '').trim();
    const groups = Array.from(form?.querySelectorAll('.forge-group') ?? []).filter(isShown)
      .map((e) => ({e, caption: caption(e.querySelector('.forge-group-header label'))}));
    if (g === null)
      return groups.map((x) => x.caption).join(', ');
    const shown = groups.find((x) => x.caption === g)?.e;
    if (shown === undefined)
      return `(no shown "${g}" group)`;
    return Array.from(shown.querySelectorAll('.forge-group-body .ui-input-root')).filter(isShown)
      .map((r) => caption(r.querySelector('label.ui-input-label'))).join(', ');
  }, group);
}

export const groupsOfTrainForm = Then('the groups of train form should be {string}', async (page: Page, list: string) => {
  await expect.poll(() => shownCaptions(page, null), {message: 'the shown groups of the train form'})
    .toBe(listOf(list).join(', '));
}, {description: 'the captions of the groups of the Train view that are shown, in order (comma-separated)'});

/* LIB-PROMOTE: a form's visible inputs, by their captions; the library claims one input at a time. The selector needs
   generalising first: it is tied to Forge's `.forge-train-form` and `.forge-group`. */
export const groupShowsInputs = Then('the {string} group should show exactly the inputs {string}', async (page: Page,
  group: string, list: string) => {
  await expect.poll(() => shownCaptions(page, group), {message: `the shown inputs of the "${group}" group`})
    .toBe(listOf(list).join(', '));
}, {description: 'the captions of the shown inputs of the train form\'s group, in order (comma-separated); a hidden ' +
  'input does not count'});

const notedReadings = new WeakMap<Page, Map<string, string>>();
const noteKey = (name: string, target: ElementRef) => `${name}|${target.phrase}`;

function notedReading(page: Page, name: string, target: ElementRef): string {
  const noted = notedReadings.get(page)?.get(noteKey(name, target));
  if (noted === undefined)
    throw new Error(`"${name}" of ${target.phrase} was not noted: "user notes the … reading of …" comes first`);
  return noted;
}

/* LIB-PROMOTE: the library's "higher / lower than remembered" compare numeric readings only (viewers.ts
   `expectRememberedDirection`); a grid cell's reading is text ("0.882"). The direction pair over numeric text is what
   is missing; "as remembered" is the library's. */
export const noteReading = When('user notes the {string} reading of {widget}', async (page: Page, name: string,
  target: ElementRef) => {
  await expect.poll(() => readingText(page, target, name), {message: `"${name}" of ${target.phrase}`}).not.toMatch(/^\(/);
  const map = notedReadings.get(page) ?? new Map<string, string>();
  notedReadings.set(page, map);
  map.set(noteKey(name, target), await readingText(page, target, name));
}, {tier: 'api', description: 'the reading\'s text, once the viewer reports it, kept for this page'});

/* Polled: a cutoff re-cuts the shown result after the view's debounce, with no busy signal of its own. */
async function expectAgainstNoted(page: Page, name: string, target: ElementRef, isHigher: boolean): Promise<void> {
  const noted = notedReading(page, name, target);
  await expect.poll(async () => {
    const text = await readingText(page, target, name);
    return (isHigher ? Number(text) > Number(noted) : Number(text) < Number(noted)) ? 'holds' : text;
  }, {message: `"${name}" of ${target.phrase} ${isHigher ? 'above' : 'below'} the noted ${noted}`}).toBe('holds');
}

export const readingLowerThanNoted = Then('the {string} reading of {widget} should be lower than noted',
  (page: Page, name: string, target: ElementRef) => expectAgainstNoted(page, name, target, false),
  {description: 'polled: the reading and the noted text as numbers, the reading strictly lower'});

export const readingHigherThanNoted = Then('the {string} reading of {widget} should be higher than noted',
  (page: Page, name: string, target: ElementRef) => expectAgainstNoted(page, name, target, true),
  {description: 'polled: the reading and the noted text as numbers, the reading strictly higher'});

/** The title of the results panel [panel]. */
async function titleOf(page: Page, panel: ElementRef): Promise<Locator> {
  return (await locate(page, panel)).first().locator('.forge-results-title').first();
}

/** The number in "(AUC = …)" of a panel's [title], as written; null without one. */
async function titleAuc(title: Locator): Promise<string | null> {
  return /\(AUC = ([\d.]+)\)/.exec((await title.textContent()) ?? '')?.[1] ?? null;
}

export const aucInTitle = Then('the AUC in the title of {element} should equal the {string} reading of {widget}',
  async (page: Page, panel: ElementRef, name: string, target: ElementRef) => {
    const title = await titleOf(page, panel);
    await expect.poll(async () => {
      const [auc, value] = await Promise.all([titleAuc(title), readingText(page, target, name)]);
      return (auc ?? '(no AUC in the title)') === value ? 'equal' : `the title's ${auc ?? '(no AUC)'} against ${value}`;
    }, {message: `the AUC in the title of ${panel.phrase} against "${name}" of ${target.phrase}`}).toBe('equal');
  }, {description: 'the number in "(AUC = …)" of the panel\'s title against the reading\'s text (both with three ' +
  'decimals), read together until they agree'});

/** The area under the curve a scatter plot draws (its y over its x, the points sorted by x and closed at (0, 0) and
 * (1, 1), as a ROC curve), by the trapezoid rule. */
function curveArea(page: Page, target: ElementRef): Promise<number> {
  return viewers.onViewer(page, target, (el) => {
    const v = DG.Widget.find(el);
    const [x, y] = [v.dataFrame.col(v.props.xColumnName), v.dataFrame.col(v.props.yColumnName)];
    const points: [number, number][] = [[0, 0], [1, 1]];
    for (let i = 0; i < v.dataFrame.rowCount; i++) {
      if (!x.isNone(i) && !y.isNone(i))
        points.push([x.get(i), y.get(i)]);
    }
    points.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    let area = 0;
    for (let i = 1; i < points.length; i++)
      area += (points[i][0] - points[i - 1][0]) * (points[i][1] + points[i - 1][1]) / 2;
    return area;
  });
}

export const curveAreaIsAuc = Then('the area under the curve of {widget} should equal the AUC in the title of {element} ' +
  'within {float}', async (page: Page, target: ElementRef, panel: ElementRef, within: number) => {
  const title = await titleOf(page, panel);
  await expect.poll(async () => {
    const [text, area] = await Promise.all([titleAuc(title), curveArea(page, target)]);
    if (text === null)
      return '(no AUC in the title)';
    const auc = Number(text);
    return Math.abs(area - auc) <= within ? 'equal' : `the curve's area ${area.toFixed(3)} against the title's ${auc}`;
  }, {message: `the area under the curve of ${target.phrase} against the AUC in the title of ${panel.phrase}`})
    .toBe('equal');
}, {description: 'the trapezoid area of the points the plot draws (its x and y columns, closed at (0, 0) and (1, 1)) ' +
  'against the number in "(AUC = …)" of the panel\'s title, polled until within the tolerance'});

/* LIB-PROMOTE: a viewer's current row is no reading of the scatter plot's status (a grid has "current row"). */
export const currentRowOf = Then('the current row of {widget} should be {int}', async (page: Page, target: ElementRef,
  row: number) => {
  await expect.poll(() => viewers.onViewer(page, target, (el) => DG.Widget.find(el).dataFrame.currentRowIdx + 1),
    {message: `the current row of ${target.phrase} (1-based)`}).toBe(row);
}, {description: 'the 1-based current row of the frame the viewer shows'});

/** Records every frame the platform's semantic type detection reports from now on, once per document: a table opened
 * after it is waited for by awaitTableReady, whenever its detection lands. */
function watchDetection(page: Page): Promise<void> {
  return page.evaluate(() => {
    const w = window as any;
    if (w.__forgeDetected !== undefined)
      return;
    w.__forgeDetected = [];
    grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      w.__forgeDetected = [...w.__forgeDetected.slice(-19), a?.args?.dataFrame?.dart];
    });
  });
}

/* A table opened in the page is ready, as the library's opener waits for, once the platform has detected its semantic
   types (watchDetection, from before the open) and its view is the current one with its own grid shown. */
async function awaitTableReady(page: Page, name: string): Promise<void> {
  await expect.poll(() => viewers.evaluate(page, (n) => {
    const w = window as any;
    const t = w.__bdd.tableNamed(n);
    if (!w.__forgeDetected.includes(t.dart))
      return 'its semantic types are not detected';
    const tv = grok.shell.tv;
    return tv?.dataFrame?.dart !== t.dart ? `the current view is ${String(tv?.name)}` :
      tv.root.querySelector('[name="viewer-Grid"]')?.getClientRects().length > 0 ? 'ready' : 'its grid is not shown';
  }, name), {timeout: pollMs(60000), message: `the table "${name}" open`}).toBe('ready');
}

/* The step ends when the table's view has its grid and the table is the current one: a menu picked before the grid is
   built left a scatter plot legend's icon over the next feature's menu. LIB-PROMOTE: the library's dataset opener reads
   files, and the manual's console line opens a generated table; with `openTable` taking a generated table as its
   source and the harness finding viewer chrome left after `closeAll` (the legend icon), this step goes. */
export const demoTableOpen = Given('a demo table {string} of {int} rows is open', async (page: Page, name: string,
  rows: number) => {
  await page.evaluate(([n, count]) => {
    const t = grok.data.demo.demog(count);
    t.name = n;
    grok.shell.addTableView(t);
  }, [name, rows] as [string, number]);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
  await expect.poll(() => page.evaluate(() => String(grok.shell.t?.name ?? '')), {message: 'the current table'})
    .toBe(name);
}, {tier: 'api', description: 'grok.data.demo.demog(rows) named and added as a table view; done when it is the ' +
  'current table'});

/* Where a scenario says "no balloon", an info balloon fails too (the library's floor counts errors and warnings).
   LIB-PROMOTE: the harness keeping the Background's balloons and errors for the first scenario (or the generator
   adding the floor) lets the Background floors go. */
export const noBalloonAtAll = Then('no balloon of any colour should have been shown', async (page: Page) => {
  const shown = await viewers.takeBalloons(page);
  expect(shown.map((b) => `${b.type}: ${b.message}`), 'balloons since the last check').toEqual([]);
}, {description: 'the platform\'s balloons of every type since the previous balloon check; checking clears them'});

/* LIB-PROMOTE: the library's task-bar watch records entry texts, not entries (a text update is recorded again), so it
   cannot count how many times a job ran. This one keeps every entry element added from now on. */
export const countTaskBarEntries = Given('user counts the task bar entries', (page: Page) => page.evaluate(() => {
  const w = window as any;
  w.__forgeTaskBar?.observer.disconnect();
  // The shell's one task bar: its entries stay under it, alone or in its tooltip of several tasks.
  const bar = document.querySelector('.d4-task-bar');
  if (bar === null)
    throw new Error('no task bar on the page');
  const entries = new Map<Element, string[]>();
  const observer = new MutationObserver(() => {
    for (const e of Array.from(bar.querySelectorAll('.d4-task-bar-entry'))) {
      const text = (e.textContent ?? '').trim();
      const texts = entries.get(e);
      if (texts === undefined)
        entries.set(e, [text]);
      else if (!texts.includes(text))
        texts.push(text);
    }
  });
  observer.observe(bar, {childList: true, subtree: true, characterData: true});
  w.__forgeTaskBar = {observer, entries};
}), {tier: 'ui', description: 'keeps every task bar entry added from now until the page resets, with the texts it showed'});

interface CountedEntries {texts: string[]; isShown: boolean; count: number}

/** The task bar entries counted since "user counts the task bar entries" that showed [text]: every text they showed,
 * whether one is still on the page (marked `data-forge-entry` for a gesture on it), and how many there were. */
function entryTexts(page: Page, text: string): Promise<CountedEntries> {
  return page.evaluate((t) => {
    const entries: Map<Element, string[]> | undefined = (window as any).__forgeTaskBar?.entries;
    if (entries === undefined)
      throw new Error('no task bar entries are counted: "user counts the task bar entries" first');
    const shown = Array.from(entries).filter(([, texts]) => texts.some((s) => s.includes(t)));
    const live = shown.find(([e]) => e.isConnected)?.[0];
    for (const e of Array.from(document.querySelectorAll('[data-forge-entry]')))
      e.removeAttribute('data-forge-entry');
    live?.setAttribute('data-forge-entry', 'true');
    return {texts: shown.flatMap(([, texts]) => texts), isShown: live !== undefined, count: shown.length};
  }, text);
}

/** Waits for a counted entry that showed [text] to be on the page; it is then marked `data-forge-entry`. */
async function awaitEntryShown(page: Page, text: string): Promise<void> {
  await expect.poll(async () => (await entryTexts(page, text)).isShown ? 'shown' : 'not shown',
    {timeout: pollMs(15000), message: `a task bar entry that showed "${text}" since "user counts the task bar entries"`})
    .toBe('shown');
}

export const taskBarEntryCount = Then('the task bar should have shown {int} entry/entries {string}', async (page: Page, n: number,
  text: string) => {
  await expect.poll(async () => (await entryTexts(page, text)).count,
    {message: `task bar entries showing "${text}" since "user counts the task bar entries"`}).toBe(n);
}, {description: 'distinct entries, each a job the task bar showed'});

/* LIB-PROMOTE: the library's "the task bar should have finished" looks for the text, and a running job replaces it with
   its progress ("Rows 1-5098 of 5098"), so it passes while the job still runs. This one follows the entry itself. */
export const taskBarEntryEnded = Then('the task bar entry that showed {string} should have ended', async (page: Page,
  text: string) => {
  await expect.poll(async () => {
    const {texts, isShown} = await entryTexts(page, text);
    return texts.length === 0 ? 'never shown' : isShown ? 'still shown' : 'ended';
  }, {message: `the task bar entry that showed "${text}" since "user counts the task bar entries"`,
    timeout: pollMs(120000)}).toBe('ended');
}, {description: 'an entry counted since "user counts the task bar entries" showed the text and every such entry has ' +
  'left the page, whatever text it shows by then: the job ended; on a 120 s budget'});

/* LIB-PROMOTE: a job's progress texts ("Rows 1-1 of 6497"), which the library's watch records only as an entry is
   added. The positive is polled; the negative is read once the entry has left the page, so no later text is missed. */
export const taskBarEntryShowedText = Then('the task bar entry that showed {string} should have shown {string}',
  async (page: Page, text: string, shown: string) => {
    await expect.poll(async () => (await entryTexts(page, text)).texts.some((s) => s.includes(shown)),
      {message: `the texts of the task bar entry that showed "${text}", one containing "${shown}"`}).toBe(true);
  }, {description: 'a text the entry counted since "user counts the task bar entries" showed at some point'});

export const taskBarEntryNotShowedText = Then('the task bar entry that showed {string} should not have shown {string}',
  async (page: Page, text: string, shown: string) => {
    await expect.poll(async () => {
      const read = await entryTexts(page, text);
      return read.count === 0 ? 'never shown' : read.isShown ? 'still shown' :
        read.texts.filter((s) => s.includes(shown)).join(' | ') || 'ended without it';
    }, {message: `the texts of the task bar entry that showed "${text}", once it has ended`}).toBe('ended without it');
  }, {description: 'once the entry has ended, none of the texts it showed contains the text'});

/* LIB-PROMOTE: the library's dataset opener reads a csv into a fresh frame, which has no creation script, so Save would
   not offer Reference. This one opens the demo file as Browse > Files does (OpenServerFile, unprocessed: the platform
   records the creation script) and names the table, so its runs are the feature's own. The step ends when the table
   is ready (awaitTableReady), as the library's opener does. */
export const openDemoFile = Given('user opens the demo file {string} as {string}', async (page: Page, file: string, name: string) => {
  await watchDetection(page);
  await page.evaluate(async ([path, tableName]) => {
    const call = DG.Func.byName('OpenServerFile').prepare({fullPath: path});
    await call.call(false, undefined, {processed: false});
    const tables: unknown = call.getOutputParamValue();
    const table: any = Array.isArray(tables) ? tables[0] : null;
    if (!(table instanceof DG.DataFrame))
      throw new Error(`OpenServerFile gave no table for ${path}`);
    table.name = tableName;
    const view = Array.from(grok.shell.tableViews as any[]).find((v: any) => v.dataFrame?.dart === table.dart);
    if (view === undefined)
      grok.shell.addTableView(table);
    else
      grok.shell.v = view;
  }, [`System:DemoFiles/${file}`, name] as [string, string]);
  await awaitTableReady(page, name);
}, {tier: 'api', description: 'OpenServerFile of System:DemoFiles/<file>, unprocessed (the table keeps its creation ' +
  'script, so Save offers Reference), named; done when the platform has detected its semantic types and its view is ' +
  'current with its grid shown'});

const FIXTURE_FIELDS = ['target', 'features', 'method', 'storage', 'tags', 'description'];

/* A saved model for a feature that is not about training (M1b's catalog and apply features): Forge:trainModel in the
   page with the view's defaults unless a row says otherwise. The model and the runs of the table are swept at feature
   end, and the sweep read back. */
export const trainThroughApiWith = Given('a Forge model {string} is trained on the current table through the API with:',
  async (page: Page, name: string, rows: string[][]) => {
    const fields: Record<string, string> = {};
    for (const [field, value] of rows) {
      if (!FIXTURE_FIELDS.includes(field))
        throw new Error(`no fixture field "${field}"; the fields are: ${FIXTURE_FIELDS.join(', ')}`);
      fields[field] = value;
    }
    await guardIdentity(page);
    const table = await page.evaluate(() => grok.shell.t ? String(grok.shell.t.name) : null);
    if (table === null)
      throw new Error('no current table to train on');
    sweepAtFeatureEnd(page, [name], [table]);
    const features = fields.features === undefined ? null : listOf(fields.features);
    await page.evaluate(async ([modelName, f, list]) => {
      const t = grok.shell.t;
      await grok.functions.call('Forge:trainModel', {table: t, target: f.target ?? t.columns.byIndex(t.columns.length - 1).name,
        features: list, method: f.method ?? null, name: modelName, storage: f.storage ?? 'none', tags: f.tags ?? null,
        description: f.description ?? null, showProgress: false});
    }, [name, fields, features] as [string, Record<string, string | undefined>, string[] | null]);
    await expectCreatedModel(name, `the Forge model "${name}" trained on ${table}`);
  }, {tier: 'api', description: 'Forge:trainModel on the current table, saved as the name, with | field | value | rows: ' +
  'target (default the last column), features (comma-separated; default the API\'s), method, storage (default none), ' +
  'tags, description; deleted with Forge:deleteModel at feature end, and the runs of the table swept'});

/* LIB-PROMOTE: a journey's scenario inherits the balloons of the ones before it (the harness removes them between
   features only); a claim about "the" balloon's link needs the earlier ones gone. */
export const balloonsClosed = Given('no balloon is shown', async (page: Page) => {
  await page.evaluate(() => DG.Balloon.closeAll());
  await expect(page.locator('.d4-balloon').filter({visible: true}), 'balloons after closing them').toHaveCount(0);
}, {tier: 'api', description: 'DG.Balloon.closeAll(), as the shell reset between features does'});

/* LIB-PROMOTE: an area of a canvas below the visible part of a scrolling panel (the Runs grid under the rows of
   results) cannot take a pointer gesture; a person scrolls to it first. */
export const scrollIntoView = When('user scrolls {element} into view', async (page: Page, target: ElementRef) => {
  await (await locate(page, target)).first().scrollIntoViewIfNeeded();
}, {tier: 'ui', description: 'scrolls the nearest scrolling ancestors until the element is in view'});

/* LIB-PROMOTE: every cancellable task-bar progress could use it; the entry's cross is the user's cancel. A running
   job replaces the entry's text with its progress ("Rows 1-1 of 6497"), so the entry is the one that showed the text
   since "user counts the task bar entries". */
export const cancelTaskBarEntry = When('user cancels the {string} task bar entry', async (page: Page, text: string) => {
  await awaitEntryShown(page, text);
  await page.locator('.d4-task-bar-entry[data-forge-entry] .fa-times').first().click();
}, {tier: 'ui', description: 'clicks the cross of the task bar entry that showed the text, counted since "user counts ' +
  'the task bar entries" and still on the page'});

/* LIB-PROMOTE: the library's hover lands on the element in one move and waits for its own mouseenter; a disabled
   button (ui.setDisabled) takes no pointer events, and the overlay that carries its reason shows the tooltip only
   after moves over it, as a person's pointer makes (probed on the Train view, 2026-10-19). */
export const slidePointerOnto = When('user slides the pointer onto {element}', async (page: Page, target: ElementRef) => {
  const box = await (await locate(page, target)).first().boundingBox();
  if (box === null)
    throw new Error(`${target.phrase}: no box to slide onto`);
  const x = box.x + box.width / 2;
  await slidePointer(page, {x, y: box.y + box.height / 2}, {x, y: Math.max(0, box.y - 40)});
}, {tier: 'ui', description: 'from 40 px above the element to its centre in five moves'});

/* A change while a training runs (scenario 5.3): the first value is typed and committed, the second goes in once the
   training of the first is in the task bar (waited for: a training queued behind another starts later). */
export const changeDuringTraining = When('user sets {string} input to {string} and then to {string} while the training runs',
  async (page: Page, caption: string, first: string, second: string) => {
    const input = el(`"${caption}" input`);
    await gestures.typeInto(page, input, first, true);
    await awaitEntryShown(page, 'Training');
    await gestures.typeInto(page, input, second, true);
  }, {tier: 'ui', description: 'types and commits the first value, and the second while a "Training …" task bar entry ' +
    'counted since "user counts the task bar entries" is shown'});
