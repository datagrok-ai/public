import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {historyHooked, openStage, placeAt, placeOf, routed, samePlace, showStage, startRouter, withRowOf}
  from '../stages';
import type {Place} from '../stages';
import {hubView, openBrought, demoBrought} from '../hub-view';
import {monitorView} from '../monitor-view';
import {PLAN_TAG} from '../plan-file';
import {exampleStudy} from './example-study';
import {PackageFunctions} from '../package';
import {modal, modalOf} from '../ui';
import {appIsOpen} from '../app-open';

const pmxViews = () => [...grok.shell.views].filter((v) => placeOf(v));
const here = () => placeAt(location.pathname + location.search)?.place;

/** The page on screen is `place`, its address too, and no page of the app is open twice. */
async function on(place: Place, what: string): Promise<void> {
  // (the state named when the wait ends, not as it starts)
  await awaitCheck(() => placeOf(grok.shell.v) === place && here() === place, '', 30000).catch(() => {
    throw new Error(`${what}: not on ${place} (shown ${placeOf(grok.shell.v)} ${grok.shell.v?.name}, address ` +
      `${location.pathname}${location.search}; views ${[...grok.shell.views].map((v) => v.name + ':' +
        (v.root.dataset.pmxStage ?? ''))}; dialogs ${DG.Dialog.getOpenDialogs().map((d) => d.title)})`);
  });
  await routed();
  const places = pmxViews().map((v) => placeOf(v));
  expect(new Set(places).size, places.length, `${what}: a page open twice (${places.join()})`);
}

category('Browser history (one view)', () => {
  test('an address names its page and its saved row', async () => {
    expect(JSON.stringify(placeAt('/apps/PmxStudio')), JSON.stringify({place: 'home'}));
    expect(placeAt('/apps/PmxStudio/plan?plan=X')!.row, 'plan=X');
    expect(placeAt('/apps/PmxStudio/monitor?study=Y&browse=apps')!.row, 'study=Y', 'the platform\'s own: no part');
    expect(placeAt('/apps/pmxstudio/analyze?study=Y')!.place, 'analyze', 'any case');
    expect(placeAt('/apps/PmxStudio/compare')!.place, 'compare');
    expect(placeAt('/'), null, 'Datagrok\'s home: not the app\'s');
    expect(placeAt('/apps/PmxStudioX'), null);
    expect(samePlace('/apps/PmxStudio/plan?plan=X', '/apps/PmxStudio/plan?plan=X&browse=apps'), true);
    expect(samePlace('/apps/PmxStudio/plan?plan=X', '/apps/PmxStudio/plan?plan=Z'), false);
    expect(samePlace('/apps/PmxStudio/monitor?study=Y', '/apps/PmxStudio/analyze?study=Y'), false);
  });

  // (QA 0.1.3 technician F14, analyst #8: the platform's push lost ?study=, or kept the page before's)
  test('the platform\'s address of a page names that page\'s saved row, or none', async () => {
    const page = (stage: string, row?: string) => ({root: {dataset: {pmxStage: stage, pmxAdded: 'true',
      ...row ? {pmxRowOf: row} : {}}}}) as unknown as DG.ViewBase;
    expect(withRowOf('/apps/PmxStudio/monitor?browse=apps', page('monitor', 'study=Y')),
      '/apps/PmxStudio/monitor?browse=apps&study=Y', 'its row added, the platform\'s own kept');
    expect(withRowOf('/apps/PmxStudio/analyze?study=CT26', page('analyze')), '/apps/PmxStudio/analyze',
      'a study not saved carries no other\'s');
    expect(withRowOf('/', page('plan', 'plan=P')), '/apps/PmxStudio/plan?plan=P', 'Datagrok\'s home: the page\'s own');
    expect(withRowOf('/apps/PmxStudio/plan?plan=Q', page('plan', 'plan=P')), '/apps/PmxStudio/plan?plan=P');
  });

  test('Back and Forward walk the stages inside the app', async () => {
    startRouter();
    const length0 = history.length;
    try {
      await showStage('home', await hubView());
      await on('home', 'Overview');
      const plan = (await exampleStudy()).getTag(PLAN_TAG)!;
      await openStage('plan', plan);
      await on('plan', 'Plan');
      await showStage('monitor', await monitorView(await exampleStudy()));
      await on('monitor', 'Monitor');
      await openStage('analyze');
      await on('analyze', 'Analyze');
      // (the browser keeps at most 50 entries: past them, Back below checks each move instead)
      expect(history.length >= Math.min(50, length0 + 3), true, 'one entry per move');
      for (const place of ['monitor', 'plan', 'home'] as const) {
        history.back();
        await on(place, `Back to ${place}`);
      }
      for (const place of ['plan', 'monitor', 'analyze'] as const) {
        history.forward();
        await on(place, `Forward to ${place}`);
      }
      // the platform's app function call for the last Forward comes late (a server call). Onto a page of the app
      // routeOpen drops it as stale; but with every page closed before it lands (this test's teardown, at once) it
      // is the app opened at that address, and Analyze would come back over the next test's Overview: let it land
      await delay(3000);
      await routed();
    } finally {
      pmxViews().forEach((v) => v.close());
    }
  }, {timeout: 180000});

  test('Back closes a dialog of the page left, as Cancel (QA 0.1.3: Compare with other studies stayed open)',
    async () => {
      startRouter();
      try {
        await showStage('home', await hubView());
        await on('home', 'Overview');
        await openStage('plan');
        await on('plan', 'Plan');
        const d = modal('test-route-dialog', 'A dialog')!;
        let closed = 0;
        d.onClose.subscribe(() => closed++);
        d.show();
        expect(modalOf('test-route-dialog') === d, true);
        history.back();
        await on('home', 'Back to the Overview');
        await awaitCheck(() => closed > 0, 'the dialog stayed open over the Overview', 5000);
        // (not expect(x, undefined): the test library's expected defaults to true, an explicit undefined too)
        expect(modalOf('test-route-dialog') === undefined, true, 'the dialog is still registered');
      } finally {
        pmxViews().forEach((v) => v.close());
      }
    }, {timeout: 120000});

  test('a study opened from the Overview replaces it: no second view', async () => {
    startRouter();
    try {
      await showStage('home', await hubView());
      await on('home', 'Overview');
      await openBrought(await demoBrought('monitor'), 'monitor');
      await on('monitor', 'the study opened');
      await delay(300);
      expect(pmxViews().some((v) => placeOf(v) === 'home'), false, 'the Overview stayed open');
      // (every page of the app is named PMX Studio: the study is named in its header)
      await awaitCheck(() => grok.shell.v?.name === 'PMX Studio', `the study's page named ${grok.shell.v?.name}`, 5000);
      history.back();
      await on('home', 'Back to the Overview');
    } finally {
      pmxViews().forEach((v) => v.close());
    }
  }, {timeout: 120000});

  test('the history hook and the app\'s listeners live only while a page of the app is open', async () => {
    // (a page shown and closed first: an earlier test may have left the app counted open with no page ever added)
    await showStage('home', await hubView());
    pmxViews().forEach((v) => v.close());
    await awaitCheck(() => !appIsOpen(), 'the app still counted open with no page of it', 5000);
    try {
      await showStage('home', await hubView());
      await on('home', 'Overview');
      expect(appIsOpen(), true);
      // (the hook may sit under the package bundle's own, put over it: the tests' bundle is a second copy of the app)
      expect(historyHooked() && history.pushState !== null, true, 'no hook while a page is open');
    } finally {
      pmxViews().forEach((v) => v.close());
    }
    // (the browser's own back on top, or, with two copies of the app, a hook that only passes through: never null)
    await awaitCheck(() => !appIsOpen() && !historyHooked() && typeof history.pushState === 'function',
      'the hook stayed on after the last page closed', 5000);
  }, {timeout: 60000});

  // (what the platform does on Back onto another stage: it runs the app function again with the path)
  test('the app run again with a page open moves that page: no second view', async () => {
    startRouter();
    try {
      await showStage('home', await hubView());
      await on('home', 'Overview');
      expect(await PackageFunctions.pmxStudioApp('/plan'), null, 'a view answered: docked as a second one');
      await on('plan', 'the app run on /plan');
      expect(pmxViews().length, 1);
    } finally {
      pmxViews().forEach((v) => v.close());
    }
  }, {timeout: 120000});
});
