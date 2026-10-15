/* How a Forge feature gets its data and its fixtures, the shown training, and the gestures on a running training. */
import type {Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, type ElementRef, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import {expectModelCount, guardIdentity, sweepModels, sweepRuns} from './grok-s.js';
import {type ResultsState, resultsState, settled, shownRun, shownSeed} from './results-state.js';

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

/** A reading as text; an absent one as the library's message, which no value equals. */
async function readingText(page: Page, target: ElementRef, name: string): Promise<string> {
  const value = await viewers.readingOf(page, target, name);
  return value instanceof viewers.MissingReading ? `(${value.text})` : String(value);
}

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
    const expected = await readingText(page, other, otherName);
    await expect.poll(() => readingText(page, target, name),
      {message: `"${name}" of ${target.phrase} against "${otherName}" of ${other.phrase} (${expected})`}).toBe(expected);
  }, {description: 'the two readings as text'});

export const readingHigher = Then('the {string} reading of {widget} should be higher than the {string} reading of {widget}',
  async (page: Page, name: string, target: ElementRef, otherName: string, other: ElementRef) => {
    const [a, b] = [await readingText(page, target, name), await readingText(page, other, otherName)];
    expect(Number(a) > Number(b), `"${name}" of ${target.phrase} (${a}) against "${otherName}" of ${other.phrase} (${b})`)
      .toBe(true);
  }, {description: 'the two readings as numbers, read once after the view is ready'});

/** The 1-based row of the catalog grid [target] whose name cell shows [name], or the names it shows. */
function catalogRow(page: Page, target: ElementRef, name: string): Promise<number | string> {
  return viewers.onViewer(page, target, (el, n) => {
    const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
    const names: string[] = [];
    for (const [key, value] of Object.entries(values)) {
      const row = /^text of cell (\d+) of name$/.exec(key);
      if (row === null)
        continue;
      if (value === n)
        return Number(row[1]);
      names.push(String(value));
    }
    return names.join(', ') || 'no rows';
  }, name);
}

/** Waits for the row of [name] in [target]; its 1-based number. */
async function catalogRowOf(page: Page, target: ElementRef, name: string): Promise<number> {
  let row: number | string = '';
  await expect.poll(async () => typeof (row = await catalogRow(page, target, name)) === 'number',
    {message: `a row "${name}" in ${target.phrase}; it shows ${row}`}).toBe(true);
  return Number(row);
}

export const catalogLists = Then('{widget} should list the Forge model {string}', async (page: Page, target: ElementRef,
  name: string) => {
  await catalogRowOf(page, target, name);
}, {description: 'a drawn cell of the name column shows the model\'s name'});

export const clickCatalogRow = When('user clicks on the Forge model {string} in {widget}', async (page: Page, name: string,
  target: ElementRef) => {
  const row = await catalogRowOf(page, target, name);
  const c = viewers.centerOf(await viewers.hitArea(page, target, `cell ${row} of name`, true));
  await page.mouse.click(c.x, c.y);
}, {tier: 'ui', description: 'clicks the name cell of the row that shows the model'});

/* LIB-PROMOTE: "the context panel should show {string}" matches String(grok.shell.o), which for a domain row is not
   its name. */
export const panelShowsModel = Then('the context panel should show the Forge model {string}', async (page: Page,
  name: string) => {
  await expect.poll(() => page.evaluate((n) => {
    const o: any = grok.shell.o;
    if (!(o instanceof DG.DomainRow) || o.typeName !== 'forge.model')
      return `the current object is ${String(o)}`;
    if (o.values.name !== n)
      return `the current object is the Forge model "${String(o.values.name)}"`;
    const panel = document.querySelector('.grok-prop-panel');
    return panel !== null && (panel.textContent ?? '').includes(n) ? 'shown' : 'the context panel does not show its name';
  }, name), {message: `the context panel and the Forge model "${name}"`}).toBe('shown');
}, {description: 'the current object is the forge.model row of that name, and the context panel shows the name'});

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

export const taskBarEntryCount = Then('the task bar should have shown {int} entry/entries {string}', async (page: Page, n: number,
  text: string) => {
  await expect.poll(() => page.evaluate((t) => {
    const entries: Map<Element, string[]> = (window as any).__forgeTaskBar?.entries ?? new Map();
    return Array.from(entries.values()).filter((texts) => texts.some((s) => s.includes(t))).length;
  }, text), {message: `task bar entries showing "${text}" since "user counts the task bar entries"`}).toBe(n);
}, {description: 'distinct entries, each a job the task bar showed'});

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

/* A saved model for a feature that is not about training (M1b's catalog and apply features): Forge:trainModel in the
   page, the view's defaults, the table's last column as the target (the view's default too), no data stored. The model
   and the runs of the table are swept at feature end, and the sweep read back. */
export const trainThroughApi = Given('a Forge model {string} is trained on the current table through the API',
  async (page: Page, name: string) => {
    await guardIdentity(page);
    const table = await page.evaluate(() => grok.shell.t ? String(grok.shell.t.name) : null);
    if (table === null)
      throw new Error('no current table to train on');
    atFeatureEnd(page, async () => {
      await sweepModels(page, [name]);
      await sweepRuns(page, table);
    });
    await page.evaluate(async (modelName) => {
      const t = grok.shell.t;
      await grok.functions.call('Forge:trainModel', {table: t, target: t.columns.byIndex(t.columns.length - 1).name,
        name: modelName, storage: 'none', showProgress: false});
    }, name);
    await expectModelCount(name, 1, `the Forge model "${name}" trained on ${table}`);
  }, {tier: 'api', description: 'Forge:trainModel {table, target: the last column, name, storage: none}; deleted with ' +
  'Forge:deleteModel at feature end, and the runs of the table swept'});

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

/* LIB-PROMOTE: a wait for a task bar entry by its text. The entry appears when the work starts; it is waited for (a
   training queued behind another starts later), never a fixed time. */
async function taskBarEntry(page: Page, text: string, message: string): Promise<Locator> {
  const entry = page.locator('.d4-task-bar-entry', {hasText: text}).first();
  await expect(entry, message).toBeVisible({timeout: pollMs(15000)});
  return entry;
}

/* LIB-PROMOTE: every cancellable task-bar progress could use it; the entry's cross is the user's cancel. */
export const cancelTaskBarEntry = When('user cancels the {string} task bar entry', async (page: Page, text: string) => {
  await (await taskBarEntry(page, text, `a task bar entry "${text}"`)).locator('.fa-times').first().click();
}, {tier: 'ui', description: 'clicks the cross of the task bar entry that shows the text'});

/* A change while a training runs (scenario 5.3): the first value is typed and committed, the second goes in once the
   training of the first is in the task bar. */
export const changeDuringTraining = When('user sets {string} input to {string} and then to {string} while the training runs',
  async (page: Page, caption: string, first: string, second: string) => {
    const input = el(`"${caption}" input`);
    await gestures.typeInto(page, input, first, true);
    await taskBarEntry(page, 'Training', 'a training in the task bar after the first value');
    await gestures.typeInto(page, input, second, true);
  }, {tier: 'ui', description: 'types and commits the first value, and the second while a "Training …" task bar entry is ' +
    'shown'});
