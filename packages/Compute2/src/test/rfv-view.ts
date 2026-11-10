import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {after, category, test, expect, waitUntil} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {provideDgViewService} from '@datagrok-libraries/webcomponents-vue';
import type {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RichFunctionView} from '../components/RFV/RichFunctionView';
import {setVueAppOptions} from '../utils';
import {awaitWebComponents, closeViewsAfter} from './utils';
import {buttons, choose, openKeyboardWorkflow, press, type} from './keyboard';
import {ribbonIcons} from './rfv';
import {clickRibbon, closedPanel, openWizard, panel} from './tree-wizard';

type Track = <T extends DG.ViewBase>(view: T) => T;

const OUTPUT_TAB_KEY = 'opened_output_tab_Compute2:TestRfvTabs';

async function openEditor(track: Track, call: DG.FuncCall) {
  await awaitWebComponents();
  const view = track(await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as DG.ViewBase);
  grok.shell.addView(view);
  await waitUntil(() => buttons(view.root, 'Run').some((el) => !el.classList.contains('d4-disabled')),
    'Run stays disabled', 15000);
  return view;
}

const openTabs = (track: Track, rows = 2) => openEditor(track, DG.Func.byName('Compute2:TestRfvTabs').prepare({rows}));

const rowsInput = (view: DG.ViewBase) => view.root.querySelector<HTMLInputElement>('input.ui-input-editor')!;

// runs and waits for the outputs; Run is gone after a run until an input changes
async function run(view: DG.ViewBase) {
  buttons(view.root, 'Run')[0].click();
  await waitUntil(() => buttons(view.root, 'Run').length === 0, 'The run did not finish', 10000);
}

async function rerunWith(view: DG.ViewBase, rows: number) {
  type(rowsInput(view), `${rows}`);
  await waitUntil(() => buttons(view.root, 'Run').some((el) => !el.classList.contains('d4-disabled')),
    'Run stays disabled', 5000);
  await run(view);
}

const dock = (view: DG.ViewBase) => view.root.querySelector('dock-spawn-ts')!.shadowRoot!;

const selectedTabs = (view: DG.ViewBase) =>
  Array.from(dock(view).querySelectorAll('.dockspan-tab-handle-selected')).map((el) => el.textContent?.trim());

// dock-spawn selects a tab on mouse down and drags it until the mouse is up
function chooseTab(view: DG.ViewBase, title: string) {
  const handle = Array.from(dock(view).querySelectorAll<HTMLElement>('.dockspan-tab-handle'))
    .find((el) => el.textContent?.trim() === title);
  if (!handle)
    throw new Error(`No tab "${title}"`);
  handle.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true}));
  window.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
}

const outputTabs = (view: DG.ViewBase) => ['first', 'middle', 'last'].filter((title) => !!panel(view, title));

const RAN: FuncCallStateInfo = {
  isRunning: false, isRunnable: true, isOutputOutdated: false, runError: undefined, pendingDependencies: [],
};

const ranTabsCall = async () => {
  const call = DG.Func.byName('Compute2:TestRfvTabs').prepare({rows: 2});
  await call.call();
  return call;
};

// a workflow step or a run opened from history gets a new call in the same RFV
async function mountRfv(track: Track, first: DG.FuncCall) {
  await awaitWebComponents();
  const view = track(DG.View.create());
  view.root.classList.remove('ui-panel', 'ui-box');
  grok.shell.addView(view);
  const call = Vue.shallowRef(first);
  const app = Vue.createApp({
    setup: () => () => Vue.h(RichFunctionView, {funcCall: call.value, callState: RAN, standalone: true}),
  });
  setVueAppOptions(app);
  provideDgViewService(app, view);
  app.mount(view.root);
  return {view, setCall: (next: DG.FuncCall) => call.value = next, unmount: () => app.unmount()};
}

// the ribbon renders menu items as elements, which Menu.find does not match by text
const menuItems = (view: DG.ViewBase, group: string) => {
  const menu = view.ribbonMenu.find(group);
  return menu ? Array.from(menu.root.querySelectorAll('.d4-menu-item')).map((el) => el.textContent?.trim()) : null;
};

category('RFV: tabs', () => {
  const track = closeViewsAfter();
  after(async () => sessionStorage.removeItem(OUTPUT_TAB_KEY));

  test('A run shows the last output tab, or the one chosen before for the function', async () => {
    sessionStorage.removeItem(OUTPUT_TAB_KEY);
    const view = await openTabs(track);
    await run(view);
    await waitUntil(() => selectedTabs(view).includes('last'),
      () => `last should be shown after the run, is ${selectedTabs(view)}`, 5000);
    chooseTab(view, 'middle');
    await waitUntil(() => sessionStorage.getItem(OUTPUT_TAB_KEY) === 'middle', 'The chosen tab was not kept', 5000);

    const again = await openTabs(track);
    await run(again);
    await waitUntil(() => selectedTabs(again).includes('middle'),
      () => `middle should be shown in the next view, is ${selectedTabs(again)}`, 5000);
  });

  test('A closed output tab stays closed after a rerun', async () => {
    const view = await openTabs(track);
    await run(view);
    await waitUntil(() => outputTabs(view).length === 3, () => `Three output tabs expected, got ${outputTabs(view)}`,
      5000);
    await closedPanel(view, 'first');
    await rerunWith(view, 3);
    expectDeepEqual(outputTabs(view), ['middle', 'last']);
  });

  test('An empty output gets its tab once a rerun fills it', async () => {
    const view = await openTabs(track, 0);
    await run(view);
    await waitUntil(() => outputTabs(view).length > 0, 'No output tabs after the run', 5000);
    expectDeepEqual(outputTabs(view), ['first', 'middle'], {prefix: 'An empty last gets no tab'});
    await rerunWith(view, 2);
    await waitUntil(() => !!panel(view, 'last'), 'A filled last should get its tab', 5000);
  });

  test('A new call in the same view brings closed tabs back and shows the chosen tab', async () => {
    sessionStorage.removeItem(OUTPUT_TAB_KEY);
    const rfv = await mountRfv(track, await ranTabsCall());
    try {
      await waitUntil(() => outputTabs(rfv.view).length === 3 && selectedTabs(rfv.view).includes('last'),
        () => `Three output tabs with last shown expected, got ${outputTabs(rfv.view)} / ${selectedTabs(rfv.view)}`,
        10000);
      await closedPanel(rfv.view, 'first');
      // as chosen in another view of the function
      sessionStorage.setItem(OUTPUT_TAB_KEY, 'middle');
      rfv.setCall(await ranTabsCall());
      await waitUntil(() => outputTabs(rfv.view).length === 3, 'The closed tab should come back', 5000);
      await waitUntil(() => selectedTabs(rfv.view).includes('middle'),
        () => `middle should be shown, is ${selectedTabs(rfv.view)}`, 5000);
    } finally {
      rfv.unmount();
    }
  });
});

category('RFV: panels', () => {
  const track = closeViewsAfter();
  after(async () => sessionStorage.removeItem(OUTPUT_TAB_KEY));

  test('The help icon opens and closes help, and closing the panel resets it', async () => {
    const view = await openTabs(track);
    const helpIcon = () => ribbonIcons(view, 'question')[0] as HTMLElement;
    expect(!!panel(view, 'Help'), false, 'Help starts closed');
    helpIcon().click();
    await waitUntil(() => !!panel(view, 'Help'), 'The icon should open help', 5000);
    expect(helpIcon().getAttribute('aria-label'), 'Close help panel');
    await closedPanel(view, 'Help');
    await waitUntil(() => helpIcon().getAttribute('aria-label') === 'Open help panel',
      'Closing the panel should reset the icon', 5000);
    helpIcon().click();
    await waitUntil(() => !!panel(view, 'Help'), 'The icon should open help again', 5000);
  });

  test('A closed form stays closed when the function runs', async () => {
    const view = await openTabs(track);
    await waitUntil(() => !!panel(view, 'Inputs')?.getClientRects().length, 'The form is not laid out', 5000);
    await closedPanel(view, 'Inputs');
    // the view takes the keys pressed with the focus lost to the page
    (document.activeElement as HTMLElement | null)?.blur();
    expect(press(document.body, 'Enter', {ctrl: true}), false, 'Ctrl+Enter should run');
    await waitUntil(() => outputTabs(view).length === 3, 'The run did not show the outputs', 10000);
    expect(!!panel(view, 'Inputs'), false, 'The run reopened the form');
  });

  test('In a workflow step the runs panel switches between opening and read-only loading', async () => {
    const view = await openWizard(track);
    await choose(view, 'first 3');
    clickRibbon(view, 'history');
    await waitUntil(() => !!panel(view, 'Runs'), 'The history icon should open the runs', 5000);
    clickRibbon(view, 'lock');
    await waitUntil(() => !!panel(view, 'Runs (read-only load)'), 'The lock icon should switch to read-only', 5000);
    clickRibbon(view, 'lock');
    await waitUntil(() => !panel(view, 'Runs (read-only load)') && !panel(view, 'Runs'),
      'The lock icon should close the read-only runs', 5000);
  });

  test('The input toggle of a step flips and is kept for the next view of the same function', async () => {
    const toggleOf = (view: DG.ViewBase) => Array.from(view.root.querySelectorAll<HTMLElement>('.rfv2-input-toggle'))
      .find((el) => el.offsetParent !== null);
    const label = (view: DG.ViewBase) => toggleOf(view)?.textContent?.trim() ?? '';
    const flipped = (text: string) =>
      text.startsWith('Show') ? text.replace('Show', 'Hide') : text.replace('Hide', 'Show');

    // item 1 gets its linked input once fixed mul sets it
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    type(view.root.querySelector<HTMLInputElement>('[dock-spawn-title="Inputs"] input.ui-input-editor')!, '7');
    await choose(view, 'item 1');
    await waitUntil(() => !!toggleOf(view), 'No input toggle', 5000);
    const initial = label(view);
    toggleOf(view)!.click();
    // the toggle state lives for the session, so leave it as it was
    let current = view;
    try {
      await waitUntil(() => label(view) === flipped(initial), () => `The toggle should flip, shows "${label(view)}"`,
        5000);
      const again = await openKeyboardWorkflow(track);
      await choose(again, 'fixed mul');
      type(again.root.querySelector<HTMLInputElement>('[dock-spawn-title="Inputs"] input.ui-input-editor')!, '7');
      await choose(again, 'item 1');
      await waitUntil(() => !!toggleOf(again), 'No input toggle', 5000);
      current = again;
      expect(label(again), flipped(initial), 'The flip was not kept');
    } finally {
      if (label(current) !== initial) {
        toggleOf(current)?.click();
        await waitUntil(() => label(current) === initial, 'The toggle did not flip back', 5000);
      }
    }
  });
});

category('RFV: exports', () => {
  const track = closeViewsAfter();

  test('Exports show after a run: the default report, the custom ones and the project', async () => {
    const view = await openEditor(track, DG.Func.byName('Compute2:TestCustomExportModel').prepare({a: 1}));
    expect(ribbonIcons(view, 'arrow-to-bottom').length, 0, 'No report before the run');
    expect(menuItems(view, 'Step exports') === null, true, 'No exports before the run');
    await run(view);
    await waitUntil(() => ribbonIcons(view, 'arrow-to-bottom').length === 1, 'The report icon should show', 5000);
    await waitUntil(() => menuItems(view, 'Step exports') !== null, 'The exports should show', 5000);
    expectDeepEqual(menuItems(view, 'Step exports'), ['Default Excel', 'rec', 'Save as project...']);
  });

  test('Fitting and sensitivity analysis icons follow the function features', async () => {
    const labels = (view: DG.ViewBase) => view.getRibbonPanels().flat()
      .flatMap((el) => Array.from(el.querySelectorAll('[aria-label]')).map((icon) => icon.getAttribute('aria-label')));
    const analyses = ['Fit inputs', 'Run sensitivity analysis'];
    const withFeatures = await openEditor(track, DG.Func.byName('Compute2:ObjectCooling2').prepare());
    expectDeepEqual(analyses.filter((label) => labels(withFeatures).includes(label)), analyses);
    const without = await openTabs(track);
    expectDeepEqual(analyses.filter((label) => labels(without).includes(label)), []);
  });
});
