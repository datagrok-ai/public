/* How a Forge feature gets its data and its fixtures, the shown training, and the gestures on a running training. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, type ElementRef, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import {catalogRowOf, modelArea} from './catalog-grid.js';
import {expectCreatedModel, guardIdentity, listOf, sweepModels, sweepRuns} from './grok-s.js';
import {readingText, type ResultsState, resultsState, settled, shownRun, shownSeed} from './results-state.js';

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
function panelOfModel(page: Page, name: string): Promise<string> {
  return page.evaluate((n) => {
    const o: any = grok.shell.o;
    if (!(o instanceof DG.DomainRow) || o.typeName !== 'forge.model')
      return `the current object is ${String(o)}`;
    if (o.values.name !== n)
      return `the current object is the Forge model "${String(o.values.name)}"`;
    const panel = document.querySelector('.grok-prop-panel');
    return panel !== null && (panel.textContent ?? '').includes(n) ? 'shown' : 'the context panel does not show its name';
  }, name);
}

/* LIB-PROMOTE: "the context panel should show {string}" matches String(grok.shell.o), which for a domain row is not
   its name. */
export const panelShowsModel = Then('the context panel should show the Forge model {string}', async (page: Page,
  name: string) => {
  await expect.poll(() => panelOfModel(page, name), {message: `the context panel and the Forge model "${name}"`})
    .toBe('shown');
}, {description: 'the current object is the forge.model row of that name, and the context panel shows the name'});

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

type Sizes = Record<string, [number, number]>;
const rememberedSizes = new WeakMap<Page, Sizes>();

function panelSizes(page: Page): Promise<Sizes> {
  return page.locator('.forge-results-panel').evaluateAll((all) => {
    const sizes: Record<string, [number, number]> = {};
    for (const p of all) {
      const r = p.getBoundingClientRect();
      const title = p.querySelector(':scope > .forge-results-title')?.textContent?.trim();
      if (title && r.width > 0)
        sizes[title] = [Math.round(r.width), Math.round(r.height)];
    }
    return sizes;
  });
}

export const rememberPanelSizes = When('user remembers the sizes of the results panels', async (page: Page) => {
  rememberedSizes.set(page, await panelSizes(page));
}, {tier: 'api', description: 'the width and height of every shown panel of the results, by its title'});

export const panelSizesKept = Then('the results panels should keep their remembered sizes', async (page: Page) => {
  const kept = rememberedSizes.get(page);
  if (kept === undefined)
    throw new Error('no panel sizes were remembered');
  expect(await panelSizes(page), 'the panels of the results, by title: [width, height]').toEqual(kept);
}, {description: 'every panel, to the pixel, read once after the change that could have moved it'});

/* LIB-PROMOTE: the library's balloon floor counts errors and warnings only; where a scenario says "no balloon", an
   info balloon is a failure too. */
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
   records the creation script) and names the table, so its runs are the feature's own. The step ends when the
   platform's semantic type detection reports the frame, as the library's opener does. */
export const openDemoFile = Given('user opens the demo file {string} as {string}', async (page: Page, file: string, name: string) => {
  await page.evaluate(async ([path, tableName, timeout]) => {
    // One subscription from before the call: the detection can report the frame before the call returns it.
    const detected: unknown[] = [];
    const waiting: (() => void)[] = [];
    let table: any = null;
    const sub = grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      detected.push(a?.args?.dataFrame?.dart);
      if (table !== null && detected.includes(table.dart)) {
        for (const done of waiting.splice(0))
          done();
      }
    });
    try {
      const call = DG.Func.byName('OpenServerFile').prepare({fullPath: path});
      await call.call(false, undefined, {processed: false});
      const tables: unknown = call.getOutputParamValue();
      table = Array.isArray(tables) ? tables[0] : null;
      if (!(table instanceof DG.DataFrame))
        throw new Error(`OpenServerFile gave no table for ${path}`);
      table.name = tableName;
      const view = Array.from(grok.shell.tableViews as any[]).find((v: any) => v.dataFrame?.dart === table.dart);
      if (view === undefined)
        grok.shell.addTableView(table);
      else
        grok.shell.v = view;
      await new Promise<void>((resolve, reject) => {
        const timer = setTimeout(() => reject(new Error(`${path}: semantic types were not detected`)), timeout);
        waiting.push(() => {
          clearTimeout(timer);
          resolve();
        });
        if (detected.includes(table.dart)) {
          for (const done of waiting.splice(0))
            done();
        }
      });
    }
    finally {
      sub.unsubscribe();
    }
  }, [`System:DemoFiles/${file}`, name, pollMs(60000)] as [string, string, number]);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
}, {tier: 'api', description: 'OpenServerFile of System:DemoFiles/<file>, unprocessed (the table keeps its creation ' +
  'script, so Save offers Reference), named; done when the platform has detected its semantic types'});

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
    atFeatureEnd(page, async () => {
      await sweepModels(page, [name]);
      await sweepRuns(page, [table]);
    });
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
