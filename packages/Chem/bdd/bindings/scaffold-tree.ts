/* The readiness barrier of the Scaffold Tree viewer: generating a tree over a molecule column runs
   for a minute or more on a hundred molecules, past the timeout a reading step waits. The barrier is
   the viewer's own "nodes" reading, not a wait of a fixed length. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

export const treeBuilt = Then('{widget} should have finished building its tree', async (page: Page, target: ElementRef) => {
  let seen = -1;
  await expect.poll(async () => {
    const r = await viewers.onViewer(page, target, (e) => {
      const values = (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {};
      return {nodes: Number(values['nodes'] ?? -1), message: String(values['message'] ?? '')};
    });
    seen = r.nodes;
    return r.nodes > 0 || /failed/i.test(r.message);
  }, {timeout: pollMs(300000), intervals: [1000], message: 'the nodes of the scaffold tree'}).toBe(true).catch(() => {
    throw new Error(`the scaffold tree built no node in five minutes; its "nodes" reading is ${seen}`);
  });
}, {description: 'polls the viewer\'s "nodes" reading for as long as a generation takes (up to five minutes)'});

export const mmpReady = Then('{widget} should have finished its analysis', async (page: Page, target: ElementRef) => {
  let seen = -1;
  await expect.poll(async () => {
    seen = await viewers.onViewer(page, target, (e) =>
      Number((window as any).__bdd.viewerOf(e).getWidgetStatus()?.values?.['substitutions'] ?? -1));
    return seen > 0;
  }, {timeout: pollMs(300000), intervals: [1000], message: 'the substitutions of the analysis'}).toBe(true).catch(() => {
    throw new Error(`the analysis found no substitution in five minutes; its "substitutions" reading is ${seen}`);
  });
}, {description: 'polls the viewer\'s "substitutions" reading for as long as a matched-pairs run takes (up to five minutes)'});
