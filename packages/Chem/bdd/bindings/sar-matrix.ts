/* The readiness barrier of the SAR Matrix viewer: building the matrices of a few thousand compounds
   runs past the timeout a reading step waits. The barrier is the viewer's own readings. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

export const matricesBuilt = Then('{widget} should have built its matrices', async (page: Page, target: ElementRef) => {
  let seen = {matrices: 0, computing: true, problem: ''};
  await expect.poll(async () => {
    seen = await viewers.onViewer(page, target, (e) => {
      const values = (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {};
      return {matrices: Number(values['matrices'] ?? 0), computing: values['computing'] === true,
        problem: String(values['problem'] ?? '')};
    });
    return (!seen.computing && seen.matrices > 0) || seen.problem !== '';
  }, {timeout: pollMs(300000), intervals: [1000], message: 'the matrices of the analysis'}).toBe(true).catch(() => {
    throw new Error(`the viewer built no matrix in five minutes; it reads ${JSON.stringify(seen)}`);
  });
  expect(seen.problem, 'the viewer\'s "problem" reading').toBe('');
}, {description: 'polls the viewer\'s "matrices" and "computing" readings for as long as a run takes (up to five minutes); a run the viewer refuses fails the step'});
