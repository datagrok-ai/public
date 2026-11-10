import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, test, expect, waitUntil} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {awaitWebComponents, closeViewsAfter, deleteSavedWorkflow} from './utils';
import {
  balloons, buttons, choose, clearBalloons, closeDialog, dialogTitles, labels, openKeyboardWorkflow, press, row,
  rowLabel, rows, type,
} from './keyboard';
import {ribbonIcons, savedRunId} from './rfv';
import {uniqueName} from './run-storage';
import {SELECTED_STEP_BACKGROUND} from '../components/TreeWizard/tree-utils';

const currentStep = (view: DG.ViewBase) =>
  new URLSearchParams((view.path ?? '').split('?')[1] ?? '').get('currentStep');

const chosen = (view: DG.ViewBase) =>
  rows(view).filter((el) => el.querySelector(`[style*="${SELECTED_STEP_BACKGROUND}"]`)).map(rowLabel);

async function chooseRow(view: DG.ViewBase, label: string) {
  const el = row(view, label);
  el.focus();
  press(el, 'Enter');
  await waitUntil(() => chosen(view).join() === label, `${label} should be chosen, is ${chosen(view)}`, 5000);
}

const panel = (view: DG.ViewBase, title: string) => view.root.querySelector(`[dock-spawn-title="${title}"]`);

// the workflow's own dock is the first, the open step has its own inside
function closePanel(view: DG.ViewBase, title: string) {
  const shadow = view.root.querySelector('dock-spawn-ts')!.shadowRoot!;
  const bar = Array.from(shadow.querySelectorAll('.panel-titlebar'))
    .find((el) => el.querySelector('.panel-titlebar-text')?.textContent?.trim() === title);
  if (bar) {
    bar.querySelector('.panel-titlebar-button-close')!.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    return;
  }
  const tab = Array.from(shadow.querySelectorAll('.dockspan-tab-handle'))
    .find((el) => el.textContent?.trim() === title);
  if (!tab)
    throw new Error(`No panel "${title}" in the workflow dock`);
  tab.querySelector<HTMLElement>('.dockspan-tab-handle-close-button')!.click();
}

async function closedPanel(view: DG.ViewBase, title: string) {
  closePanel(view, title);
  await waitUntil(() => !panel(view, title), `The ${title} panel should close`, 5000);
}

const WIZARD = ['Wizard title', 'first 3', 'group title', 'inner', 'review', 'last'];

async function openWizard(track: <T extends DG.ViewBase>(view: T) => T, resolve?: (data?: any) => void) {
  await awaitWebComponents();
  const call = DG.Func.byName('Compute2:MockTreeWizardPipeline').prepare();
  if (resolve)
    call.aux.resolve = resolve;
  const view = track(await grok.functions.call('Compute2:TreeWizardEditor', {call}) as DG.ViewBase);
  grok.shell.addView(view);
  await waitUntil(() => labels(view).join() === WIZARD.join(),
    () => `The tree should show ${WIZARD}, shows ${labels(view)}`, 15000);
  return view;
}

const clickRibbon = (view: DG.ViewBase, icon: string) => (ribbonIcons(view, icon)[0] as HTMLElement).click();

// the compositor overlay covers the workflow while it is busy, so a user can't press its buttons then
const idle = (view: DG.ViewBase) => waitUntil(
  () => !view.root.querySelector('.c2-compositor-overlay--active'), 'The workflow stays busy', 10000);

const comparisonNames = () => Array.from(document.querySelectorAll<HTMLInputElement>(
  '.d4-dialog .c2-comparison-add-names input')).map((el) => el.value);

const saveRunName = () => document.querySelector<HTMLInputElement>('.d4-dialog .c2-runs-field input')?.value;

category('TreeWizard: routing', () => {
  const track = closeViewsAfter();

  test('The chosen step is kept in the view path', async () => {
    const view = await openKeyboardWorkflow(track);
    expect(chosen(view).join(), 'Keyboard');
    await waitUntil(() => currentStep(view) === '0', `currentStep should be "0", is ${currentStep(view)}`, 5000);
    await choose(view, 'item 2');
    await waitUntil(() => currentStep(view) === '0 2 1', `currentStep should be "0 2 1", is ${currentStep(view)}`,
      5000);
    expect(chosen(view).join(), 'item 2');
    await chooseRow(view, 'Keyboard');
    await waitUntil(() => currentStep(view) === '0', `currentStep should be "0", is ${currentStep(view)}`, 5000);
  });

  test('A new workflow opens on its root even when the URL has the step of another view', async () => {
    const other = await openWizard(track);
    await choose(other, 'inner');
    await waitUntil(() => window.location.search.includes('currentStep=0+1+0'),
      `The URL should have the step of inner, is ${window.location.search}`, 5000);
    const view = await openWizard(track);
    await waitUntil(() => currentStep(view) === '0', `currentStep should be "0", is ${currentStep(view)}`, 5000);
    expect(chosen(view).join(), 'Wizard title');
  });

  test('Removing the chosen step chooses the step that took its place', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'item 2');
    press(row(view, 'item 2'), 'Delete');
    closeDialog('ok');
    await waitUntil(() => !labels(view).includes('item 2'), 'item 2 should be removed', 5000);
    await waitUntil(() => chosen(view).join() === 'item 3', `item 3 should be chosen, is ${chosen(view)}`, 5000);
    expect(currentStep(view), '0 2 1');
  });

  test('A run opened by id shows its name and id until it is edited', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:MockPipeline1');
    const model = func.friendlyName;
    const view = track(await grok.functions.call('Compute2:TreeWizardEditor', {call: func.prepare()}) as DG.ViewBase);
    await waitUntil(() => ribbonIcons(view, 'save').length > 0, 'Expected a save icon', 15000);
    (ribbonIcons(view, 'save')[0] as HTMLElement).click();
    await waitUntil(() => document.querySelector('.d4-dialog .c2-runs-field input') !== null,
      'Save dialog not shown', 5000);
    const name = uniqueName('Opened by id');
    const input = document.querySelector<HTMLInputElement>('.d4-dialog .c2-runs-field input')!;
    input.value = name;
    input.dispatchEvent(new Event('input'));
    await Vue.nextTick();
    closeDialog('save');
    await waitUntil(() => savedRunId(view) != null, 'URL should carry the saved run id', 10000);
    const id = savedRunId(view)!;
    try {
      const call = func.prepare();
      call.aux.initialRunId = id;
      const opened = track(await grok.functions.call('Compute2:TreeWizardEditor', {call}) as DG.ViewBase);
      grok.shell.addView(opened);
      await waitUntil(() => opened.name === `${model} - ${name}` && savedRunId(opened) === id,
        () => `The run should open as "${model} - ${name}" with its id, is "${opened.name}" at ${opened.path}`, 15000);
      clickRibbon(opened, 'exchange');
      await waitUntil(() => comparisonNames().length === 1, 'No add to comparison dialog', 5000);
      expect(comparisonNames()[0], name, 'An unedited run should be compared as the saved run');
      closeDialog('cancel');

      await choose(opened, 'step1');
      const inputs = Array.from(panel(opened, 'Inputs')!.querySelectorAll<HTMLInputElement>('input.ui-input-editor'))
        .filter((el) => el.offsetParent !== null);
      type(inputs[0], `${Number(inputs[0].value || 0) + 1}`);
      await waitUntil(() => opened.name === `${model} - ${name} (edited)` && savedRunId(opened) === null,
        () => `An edit should mark the name and drop the id, is "${opened.name}" at ${opened.path}`, 5000);
      clickRibbon(opened, 'exchange');
      await waitUntil(() => comparisonNames().length === 1, 'No add to comparison dialog', 5000);
      expect(comparisonNames()[0].startsWith('Pipeline 1 title (current, '), true,
        `An edited run should be compared as a snapshot, is "${comparisonNames()[0]}"`);
      closeDialog('cancel');
    } finally {
      await deleteSavedWorkflow(id);
    }
  });
});

category('TreeWizard: panels', () => {
  const track = closeViewsAfter();

  test('A closed step panel shows again when a step is chosen in the tree', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    await closedPanel(view, 'Step review');
    await chooseRow(view, 'fixed mul');
    await waitUntil(() => !!panel(view, 'Step review'), 'Choosing the same step should show its panel', 5000);
    await closedPanel(view, 'Step review');
    await chooseRow(view, 'fixed add');
    await waitUntil(() => !!panel(view, 'Step review'), 'Choosing another step should show its panel', 5000);
  });

  test('A closed workflow panel shows again when a workflow is chosen in the tree', async () => {
    const view = await openKeyboardWorkflow(track);
    await waitUntil(() => !!panel(view, 'Step sequence review'), 'The root should show its panel', 5000);
    await closedPanel(view, 'Step sequence review');
    await chooseRow(view, 'Keyboard');
    await waitUntil(() => !!panel(view, 'Step sequence review'), 'Choosing the root should show its panel', 5000);
  });

  test('Closing the inspector panel resets its ribbon toggle', async () => {
    const view = await openKeyboardWorkflow(track);
    const toggle = () => ribbonIcons(view, 'bug')[0] as HTMLElement;
    expect(!!panel(view, 'Inspector'), false, 'The inspector starts hidden');
    toggle().click();
    await waitUntil(() => !!panel(view, 'Inspector'), 'The toggle should show the inspector', 5000);
    await closedPanel(view, 'Inspector');
    toggle().click();
    await waitUntil(() => !!panel(view, 'Inspector'), 'After its panel closed the toggle should show it again', 5000);
  });

  test('Back chooses the previous step and moves the focus into it', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    buttons(view.root, 'Back')[0].click();
    await waitUntil(() => chosen(view).join() === 'fixed add', `fixed add should be chosen, is ${chosen(view)}`, 5000);
    await waitUntil(() => view.root.contains(document.activeElement) && !document.activeElement!.closest('.mtl-tree'),
      `The focus should be in the opened step, is on ${document.activeElement?.tagName}`, 10000);
  });
});

category('TreeWizard: confirmations', () => {
  const track = closeViewsAfter();

  test('Confirm data asks first, then returns the result and closes the view', async () => {
    const returned: any[] = [];
    const view = await openWizard(track, (data) => returned.push(data));
    clickRibbon(view, 'check');
    expect(dialogTitles().includes('Action confirmation'), true, 'No return confirmation');
    closeDialog('cancel');
    expect(returned.length, 0, 'Cancel returned');
    clickRibbon(view, 'check');
    closeDialog('ok');
    // the editor resolves once more, without a result, when its view is closed
    await waitUntil(() => returned.length === 2, () => `The view should return and close, got ${returned.length}`,
      5000);
    expectDeepEqual(returned, [{a: 3}, undefined]);
  });

  test('An action asks first only when it has a confirmation message', async () => {
    const view = await openWizard(track);
    await choose(view, 'first 3');
    buttons(view.root, 'Double a')[0].click();
    expect(dialogTitles().includes('Action confirmation'), false, 'Double a has no confirmation message');
    await waitUntil(() => labels(view).includes('first 6'), () => `Double a should run, shows ${labels(view)}`, 5000);
    await idle(view);
    buttons(view.root, 'Reset a')[0].click();
    expect(dialogTitles().includes('Action confirmation'), true, 'Reset a should ask first');
    closeDialog('cancel');
    // actions run in order, so a reset run by Cancel would show as first 2
    buttons(view.root, 'Double a')[0].click();
    await waitUntil(() => labels(view).includes('first 12'), () => `Cancel should not reset, shows ${labels(view)}`,
      5000);
    await idle(view);
    buttons(view.root, 'Reset a')[0].click();
    closeDialog('ok');
    await waitUntil(() => labels(view).includes('first 1'), () => `OK should reset, shows ${labels(view)}`, 5000);
  });
});

category('TreeWizard: saving and comparing', () => {
  const track = closeViewsAfter();

  test('The workflow is added to a comparison under its title', async () => {
    const view = await openWizard(track);
    clickRibbon(view, 'exchange');
    await waitUntil(() => comparisonNames().length === 1, 'No add to comparison dialog', 5000);
    try {
      expect(comparisonNames()[0].startsWith('Wizard title (current, '), true, `Named "${comparisonNames()[0]}"`);
    } finally {
      closeDialog('cancel');
    }
  });

  test('Saving a nested workflow starts from its title', async () => {
    const view = await openWizard(track);
    await chooseRow(view, 'group title');
    row(view, 'group title').querySelector('.fa-save')!.closest<HTMLElement>('dg-icon-fa')!.click();
    await waitUntil(() => saveRunName() !== undefined, 'No Save run dialog', 5000);
    try {
      expect(saveRunName(), 'group title');
    } finally {
      closeDialog('cancel');
    }
  });

  test('Saving a step to its history starts from its title', async () => {
    const view = await openWizard(track);
    await choose(view, 'first 3');
    press(row(view, 'first 3'), 'Enter', {ctrl: true});
    await waitUntil(() => ribbonIcons(view, 'cloud-upload-alt').length > 0,
      'The run step should offer saving to its history', 10000);
    clickRibbon(view, 'cloud-upload-alt');
    await waitUntil(() => saveRunName() !== undefined, 'No Save run dialog', 5000);
    try {
      expect(saveRunName(), 'first 3');
    } finally {
      closeDialog('cancel');
    }
  });
});

category('TreeWizard: navigation', () => {
  const track = closeViewsAfter();

  test('Next on a workflow goes into its first step, on an action step to the next step', async () => {
    const view = await openWizard(track);
    await chooseRow(view, 'group title');
    buttons(view.root, 'Next')[0].click();
    await waitUntil(() => chosen(view).join() === 'inner', () => `inner should be chosen, is ${chosen(view)}`, 5000);
    await chooseRow(view, 'review');
    buttons(view.root, 'Next')[0].click();
    await waitUntil(() => chosen(view).join() === 'last', () => `last should be chosen, is ${chosen(view)}`, 5000);
  });

  test('Next on the last step says what to do next', async () => {
    const view = await openWizard(track);
    await choose(view, 'last');
    press(row(view, 'last'), 'Enter', {ctrl: true});
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'last should run and offer Next', 10000);
    clearBalloons();
    buttons(view.root, 'Next')[0].click();
    expect(balloons().some((text) => text.includes('Now you can')), true, 'No hint after the last step');
    expect(chosen(view).join(), 'last');
  });
});

category('TreeWizard: export', () => {
  const track = closeViewsAfter();

  test('Export shows once every step has run, with the default and the custom exports', async () => {
    const view = await openWizard(track);
    const exportMenu = () => view.ribbonMenu.find('Export');
    // the ribbon renders menu items as elements, which Menu.find does not match by text
    const exportItems = () => Array.from(exportMenu().root.querySelectorAll<HTMLElement>('.d4-menu-item'))
      .map((el) => el.textContent?.trim());
    expect(exportMenu() == null, true, 'Export should wait for the steps to run');
    clickRibbon(view, 'forward');
    await waitUntil(() => exportMenu() != null, 'Export should show once every step has run', 10000);
    expectDeepEqual(exportItems(), ['Default Excel', 'Record export']);
  });
});
