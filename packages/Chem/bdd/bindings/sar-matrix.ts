/* The readiness barriers of the SAR Matrix viewer: building the matrices of a few thousand compounds
   runs past the timeout a reading step waits. The barrier is the viewer's own readings, whose
   "computing" holds from the moment a rebuild is queued. Both are waits, so a guide leaves them out. Then a reading
   held against what the Summary shows, both read in one poll, since the series names are numbered by the build. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {type ElementRef, escapeRegExp, expect, locate, pollMs, silent, viewers} from '@datagrok-libraries/bdd/runtime';

type Run = {matrices: number, computing: boolean, problem: string, failed: string};

async function waitForRun(page: Page, target: ElementRef, done: (run: Run) => boolean, what: string): Promise<Run> {
  let seen: Run = {matrices: 0, computing: true, problem: '', failed: ''};
  await expect.poll(async () => {
    seen = await viewers.onViewer(page, target, (e) => {
      const values = (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {};
      return {matrices: Number(values['matrices'] ?? 0), computing: values['computing'] === true,
        problem: String(values['problem'] ?? ''), failed: String(values['failed'] ?? '')};
    });
    // A failed build is final: the viewer queues nothing after it until a setting changes.
    return done(seen) || (!seen.computing && seen.failed !== '');
  }, {timeout: pollMs(300000), intervals: [1000], message: what}).toBe(true).catch(() => {
    throw new Error(`the viewer did not finish ${what} in five minutes; it reads ${JSON.stringify(seen)}`);
  });
  expect(seen.failed, 'the viewer\'s "failed" reading: the build threw').toBe('');
  return seen;
}

export const matricesBuilt = Then('{widget} should have built its matrices', async (page: Page, target: ElementRef) => {
  silent(page);
  const run = await waitForRun(page, target, (r) => (!r.computing && r.matrices > 0) || r.problem !== '',
    'building its matrices');
  expect(run.problem, 'the viewer\'s "problem" reading').toBe('');
}, {description: 'polls the viewer\'s "matrices" and "computing" readings for as long as a run takes (up to five minutes); a run the viewer refuses, or one that failed, fails the step at once with its message'});

export const analysisFinished = Then('{widget} should have finished building', async (page: Page, target: ElementRef) => {
  silent(page);
  await waitForRun(page, target, (r) => !r.computing, 'its analysis');
}, {description: 'waits (up to five minutes) until the viewer has no run under way or queued, whatever it found: the barrier before claiming an empty or refused analysis; a build that failed fails it with the viewer\'s "failed" message'});

/** The visible text of an element, '' while it has none to read. */
async function shownText(page: Page, element: ElementRef): Promise<string> {
  const loc = (await locate(page, element)).filter({visible: true}).first();
  return ((await loc.textContent({timeout: 1000}).catch(() => null)) ?? '').trim();
}

export const readingIsTextOf = Then('the {string} reading of {widget} should be the text of {element}',
  async (page: Page, name: string, target: ElementRef, element: ElementRef) => {
    await expect.poll(async () => {
      const reading = await viewers.readingOf(page, target, name);
      if (reading instanceof viewers.MissingReading)
        return String(reading);
      const text = await shownText(page, element);
      return text !== '' && text === String(reading) ? 'same' : `the reading is "${String(reading)}", ${element.phrase} shows "${text}"`;
    }, {message: `"${name}" reading of ${target.phrase} against the text of ${element.phrase}`}).toBe('same');
  }, {description: 'a reading equal to the text an element shows, both read in the same poll: the claim holds whatever the names are numbered'});

export const textHoldsReading = Then('{element} should contain the {string} reading of {widget}',
  async (page: Page, element: ElementRef, name: string, target: ElementRef) => {
    await expect.poll(async () => {
      const reading = await viewers.readingOf(page, target, name);
      if (reading instanceof viewers.MissingReading)
        return String(reading);
      const text = await shownText(page, element);
      const word = new RegExp(`(^|\\W)${escapeRegExp(String(reading))}(?!\\w)`);
      return String(reading) !== '' && word.test(text) ? 'contained' : `${element.phrase} shows "${text}", the reading is "${String(reading)}"`;
    }, {message: `the text of ${element.phrase} against the "${name}" reading of ${target.phrase}`}).toBe('contained');
  }, {description: 'the reading as a whole word of what the element shows ("Series 5" is not in "Series 52_1 · Core 1")'});
