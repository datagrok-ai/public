import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, test, expect, waitUntil} from '@datagrok-libraries/test/src/test';
import {closeViewsAfter, awaitWebComponents} from './utils';
import {RunDetailsForm} from '../components/RunStorage/RunDetailsForm';

// Keyboard tests need a view in the shell: focus and the dock's panel placement only work on attached elements.

type Keys = {alt?: boolean, ctrl?: boolean, shift?: boolean};

export const press = (el: Element, key: string, keys: Keys = {}) => el.dispatchEvent(new KeyboardEvent('keydown', {
  key, altKey: !!keys.alt, ctrlKey: !!keys.ctrl, shiftKey: !!keys.shift, bubbles: true, cancelable: true,
}));

export const rowLabel = (row: Element) => row.querySelector('.mtl-ml')?.textContent?.trim() ?? '';
export const rows = (view: DG.ViewBase) => Array.from(view.root.querySelectorAll<HTMLElement>('.mtl-tree .tree-node'));
export const labels = (view: DG.ViewBase) => rows(view).map(rowLabel);
const levels = (view: DG.ViewBase) => rows(view).map((row) => row.getAttribute('aria-level'));
export const row = (view: DG.ViewBase, label: string) => {
  const found = rows(view).find((el) => rowLabel(el) === label);
  if (!found)
    throw new Error(`No tree row "${label}" in ${JSON.stringify(labels(view))}`);
  return found;
};
export const chosen = (view: DG.ViewBase) =>
  rows(view).filter((el) => el.querySelector('[aria-current="step"]')).map(rowLabel);
const focusedLabel = () => {
  const el = document.activeElement as HTMLElement | null;
  return el?.classList.contains('tree-node') ? rowLabel(el) : null;
};

export const dialogTitles = () =>
  Array.from(document.querySelectorAll('.d4-dialog .d4-dialog-header')).map((el) => el.textContent?.trim());
export const closeDialog = (button: string) => {
  const found = Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog button'))
    .find((el) => el.textContent?.trim().toLowerCase() === button);
  if (!found)
    throw new Error(`No dialog button "${button}"`);
  found.click();
};

export const balloons = () =>
  Array.from(document.querySelectorAll('.d4-balloon')).map((el) => el.textContent?.trim() ?? '');
export const clearBalloons = () => document.querySelectorAll('.d4-balloon').forEach((el) => el.remove());

const describe = (el: Element | null) => el ? `${el.tagName}.${el.className}` : 'nothing';

// the ribbon icons that take the focus, in order
const ribbonStops = (view: DG.ViewBase) => view.getRibbonPanels().flat()
  .flatMap((panel) => Array.from(panel.querySelectorAll<HTMLElement>('[tabindex="0"]')));

const panelOf = (el: Element | null) =>
  el?.closest('[dock-spawn-title]')?.getAttribute('dock-spawn-title') ?? describe(el);

// presses Tab from the start until the focus is back on it: the panels it went through, sorted
function tabRound(start: HTMLElement) {
  start.focus();
  const panels = new Set<string>();
  for (let i = 0; i < 200; i++) {
    press(document.activeElement!, 'Tab');
    panels.add(panelOf(document.activeElement));
    if (document.activeElement === start)
      return JSON.stringify([...panels].sort());
  }
  throw new Error(`Tab never came back to ${describe(start)}, went through ${JSON.stringify([...panels])}`);
}

function expectTabEnters(from: HTMLElement, panels: string[]) {
  from.focus();
  press(from, 'Tab');
  expect(panels.includes(panelOf(document.activeElement)), true,
    `Tab from ${describe(from)} should go into ${panels.join(', ')}, is on ${describe(document.activeElement)}`);
}

const treeStop = (view: DG.ViewBase) => rows(view).find((el) => el.tabIndex === 0)!;

export const type = (input: HTMLInputElement, value: string) => {
  input.focus();
  input.value = value;
  input.dispatchEvent(new Event('input', {bubbles: true}));
  input.dispatchEvent(new Event('change', {bubbles: true}));
};

const nextFrame = () => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));

export const buttons = (root: HTMLElement, label: string) =>
  Array.from(root.querySelectorAll<HTMLButtonElement>('button'))
    .filter((el) => el.textContent?.trim().toLowerCase() === label.toLowerCase());

const INITIAL = [
  'Keyboard', 'fixed add', 'fixed mul', 'list', 'item 1', 'item 2', 'item 3', 'pinned', 'slow', 'group', 'grouped step',
];
const listed = (view: DG.ViewBase) => labels(view).slice(4, 8).join();

export async function openKeyboardWorkflow(track: <T extends DG.ViewBase>(view: T) => T) {
  await awaitWebComponents();
  const call = DG.Func.byName('Compute2:MockKeyboardPipeline').prepare();
  const view = track(await grok.functions.call('Compute2:TreeWizardEditor', {call}) as unknown as DG.ViewBase);
  grok.shell.addView(view);
  await waitUntil(() => JSON.stringify(labels(view)) === JSON.stringify(INITIAL),
    `Tree should show ${JSON.stringify(INITIAL)}`, 15000);
  return view;
}

// waits for the form of the newly opened step, not the one it replaces
export async function choose(view: DG.ViewBase, label: string) {
  const inputsPanel = () => view.root.querySelector('[dock-spawn-title="Inputs"]');
  const before = inputsPanel();
  const el = row(view, label);
  el.focus();
  press(el, 'Enter');
  await waitUntil(() => {
    const panel = inputsPanel();
    return !!panel && panel !== before &&
      Array.from(panel.querySelectorAll('input')).some((el) => el.offsetParent !== null);
  }, `The step "${label}" did not open`, 10000);
}

// the disabled button must be what the pointer hits, once the dock has laid out the view
async function expectHoverable(button: HTMLElement) {
  let hit: Element | null = null;
  await waitUntil(() => {
    button.scrollIntoView({block: 'center'});
    const box = button.getBoundingClientRect();
    hit = document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2);
    return !!hit && button.contains(hit);
  }, () => `The pointer should hit the disabled button, hits ${describe(hit)}`, 3000);
}

// from the open fixed mul: the link sets item 1's b from its a, then b, behind the linked inputs toggle, is set apart
async function makeItem1Inconsistent(view: DG.ViewBase) {
  const visible = <T extends HTMLElement>(selector: string) =>
    Array.from(view.root.querySelectorAll<T>(selector)).filter((el) => el.offsetParent !== null);
  const stepInputs = () => visible<HTMLInputElement>('[dock-spawn-title="Inputs"] input.ui-input-editor');

  await waitUntil(() => stepInputs().length > 0, 'No fixed mul inputs', 5000);
  type(stepInputs()[0], '7');
  await choose(view, 'item 1');
  await waitUntil(() => visible('.rfv2-input-toggle').length > 0, 'No linked inputs toggle', 5000);
  // the toggle stays open for the session, so an earlier test may have left b shown
  if (stepInputs().length < 2)
    visible('.rfv2-input-toggle')[0].click();
  await waitUntil(() => stepInputs().length === 2 && stepInputs()[1].value === '7', 'b should show 7 from the link',
    5000);
  type(stepInputs()[1], '5');
  await waitUntil(() => buttons(view.root, 'Update').length > 0, 'Changing the linked input b should offer Update',
    5000);
  return stepInputs()[1];
}

// he-tree's own moves show on the next render and the driver applies moves in order, so the keys moved nothing when the
// tree is the same a frame later and a real move pressed after them is the only change
async function expectNoMove(view: DG.ViewBase, pressKeys: () => void) {
  const levelsBefore = levels(view);
  pressKeys();
  await nextFrame();
  expect(JSON.stringify(labels(view)), JSON.stringify(INITIAL), 'The tree order changed');
  expect(JSON.stringify(levels(view)), JSON.stringify(levelsBefore), 'The nesting changed');
  press(row(view, 'item 1'), 'ArrowDown', {alt: true});
  const moved = [...INITIAL.slice(0, 4), 'item 2', 'item 1', ...INITIAL.slice(6)];
  await waitUntil(() => JSON.stringify(labels(view)) === JSON.stringify(moved),
    () => `Only item 1 should move, got ${JSON.stringify(labels(view))}`, 5000);
}

category('TreeWizard: keyboard moves', () => {
  const track = closeViewsAfter();

  test('Alt+Down moves a list item in the workflow and keeps the focus on it', async () => {
    const view = await openKeyboardWorkflow(track);
    const item = row(view, 'item 1');
    item.focus();
    press(item, 'ArrowDown', {alt: true});
    await waitUntil(() => listed(view) === 'item 2,item 1,item 3,pinned',
      `item 1 should move below item 2, got ${JSON.stringify(labels(view))}`, 5000);
    await waitUntil(() => focusedLabel() === 'item 1', `The focus should stay on item 1, is on ${focusedLabel()}`, 5000);

    press(row(view, 'item 1'), 'ArrowUp', {alt: true});
    await waitUntil(() => listed(view) === 'item 1,item 2,item 3,pinned',
      `item 1 should move back up, got ${JSON.stringify(labels(view))}`, 5000);
    await waitUntil(() => focusedLabel() === 'item 1', `The focus should stay on item 1, is on ${focusedLabel()}`, 5000);
  });

  // he-tree moves a step on its own Alt+arrows, ignoring the drag rules and the driver
  test('Alt+arrows move no step that cannot be dragged, and never change the nesting', async () => {
    const view = await openKeyboardWorkflow(track);
    await expectNoMove(view, () => {
      press(row(view, 'fixed add'), 'ArrowDown', {alt: true});
      press(row(view, 'list'), 'ArrowUp', {alt: true});
      press(row(view, 'pinned'), 'ArrowUp', {alt: true});
      press(row(view, 'item 2'), 'ArrowLeft', {alt: true});
      press(row(view, 'item 2'), 'ArrowRight', {alt: true});
    });
  });
});

category('TreeWizard: collapsed workflows', () => {
  const track = closeViewsAfter();

  // he-tree rebuilds its rows from every new state, and the tree unmounts while its panel is hidden
  test('A collapsed workflow stays collapsed when the tree rebuilds or its panel is hidden and shown', async () => {
    const view = await openKeyboardWorkflow(track);
    const group = row(view, 'group');
    group.focus();
    press(group, 'ArrowLeft');
    await waitUntil(() => !labels(view).includes('grouped step'), 'group did not collapse', 5000);

    press(row(view, 'item 1'), 'ArrowDown', {alt: true});
    await waitUntil(() => listed(view) === 'item 2,item 1,item 3,pinned', 'item 1 did not move', 5000);
    expect(labels(view).includes('grouped step'), false, 'The rebuild opened group');

    const treeToggle = () => ribbonStops(view).find((el) => el.classList.contains('fa-folder-tree'))!;
    treeToggle().click();
    await waitUntil(() => rows(view).length === 0, 'The tree did not hide', 5000);
    treeToggle().click();
    await waitUntil(() => labels(view).includes('group'), 'The tree did not come back', 5000);
    expect(labels(view).includes('grouped step'), false, 'Hiding and showing the tree opened group');
  });
});

category('TreeWizard: keyboard while busy', () => {
  const track = closeViewsAfter();

  test('A move pressed during a run changes nothing', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'slow');
    press(row(view, 'slow'), 'Enter', {ctrl: true});
    await waitUntil(() => view.root.querySelector('.c2-compositor-overlay--active') !== null, 'The run did not start',
      5000);
    press(row(view, 'item 1'), 'ArrowDown', {alt: true});
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'The run did not finish', 15000);
    await expectNoMove(view, () => {});
  });
});

category('TreeWizard: keyboard removal', () => {
  const track = closeViewsAfter();

  test('Delete then OK removes the item and the focus moves to the next row', async () => {
    const view = await openKeyboardWorkflow(track);
    const item = row(view, 'item 2');
    item.focus();
    press(item, 'Delete');
    // the confirmation shows at once, which the check for steps that cannot be removed relies on
    expect(dialogTitles().includes('Remove confirmation'), true, 'No remove confirmation');
    closeDialog('ok');
    await waitUntil(() => !labels(view).includes('item 2'), 'item 2 should be removed', 5000);
    await waitUntil(() => focusedLabel() === 'item 3',
      `The focus should move to the next row, is on ${focusedLabel() ?? document.activeElement?.tagName}`, 5000);
  });

  test('Delete on a step that cannot be removed does nothing', async () => {
    const view = await openKeyboardWorkflow(track);
    press(row(view, 'pinned'), 'Delete');
    press(row(view, 'fixed add'), 'Delete');
    expect(dialogTitles().includes('Remove confirmation'), false, 'No confirmation for fixed or pinned steps');
  });
});

category('TreeWizard: keyboard shortcuts', () => {
  const track = closeViewsAfter();

  test('Ctrl+Enter on an invalid step warns and does not run', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed add');
    clearBalloons();
    press(row(view, 'fixed add'), 'Enter', {ctrl: true});
    await waitUntil(() => balloons().some((text) => text.includes('Fix the invalid inputs')), 'No warning', 5000);
    expect(buttons(view.root, 'Next').length, 0, 'The step must not run');
  });

  test('Ctrl+Shift+Enter presses Update of an inconsistent step, whose confirmation ignores Tab', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true, shift: true});
    expect(dialogTitles().includes('Update confirmation'), false, 'fixed mul has nothing to update');

    const b = await makeItem1Inconsistent(view);
    press(b, 'Enter', {ctrl: true, shift: true});
    // shown at once, so the check above that nothing showed holds
    expect(dialogTitles().includes('Update confirmation'), true, 'No update confirmation');
    const dialog = Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog'))
      .find((el) => el.querySelector('.d4-dialog-header')?.textContent?.trim() === 'Update confirmation')!;
    // Enter confirms and Esc cancels wherever the focus is, so Tab moves it nowhere
    const button = dialog.querySelector<HTMLElement>('button')!;
    button.focus();
    for (const shift of [false, true]) {
      expect(press(button, 'Tab', {shift}), false, `${shift ? 'Shift+' : ''}Tab should be swallowed`);
      expect(document.activeElement === button, true, `${shift ? 'Shift+' : ''}Tab moved the focus`);
    }
    closeDialog('cancel');
  });

  test('Ctrl+Enter on a focused row runs that step, not the open one', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed add');
    clearBalloons();
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true});
    // pressing the open invalid step would warn at once
    expect(balloons().some((text) => text.includes('Fix the invalid inputs')), false,
      'The open invalid step was pressed');
    await choose(view, 'fixed mul');
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'fixed mul should have run', 10000);
  });

  test('Ctrl+Shift+Enter on a focused row updates that step, or the substeps of a workflow', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    await makeItem1Inconsistent(view);
    await choose(view, 'fixed mul');
    // item 1 waits for fixed mul until it runs, and the workflow row offers its update only then
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true});
    await waitUntil(() => !!row(view, 'list').querySelector('[aria-label="Update substeps with consistent values"]'),
      'The list row should offer to update its substeps', 10000);
    for (const label of ['item 1', 'list']) {
      press(row(view, label), 'Enter', {ctrl: true, shift: true});
      await waitUntil(() => dialogTitles().includes('Update confirmation'), `No update confirmation for ${label}`, 5000);
      closeDialog('cancel');
      await waitUntil(() => !dialogTitles().includes('Update confirmation'), 'The confirmation stays open', 5000);
    }
  });

  test('Ctrl+S saves with the focus anywhere in the view or lost, and never opens the browser save', async () => {
    const view = await openKeyboardWorkflow(track);
    const saveDialogs = () => dialogTitles().filter((title) => title === 'Save run').length;
    press(row(view, 'fixed mul'), 's', {ctrl: true});
    await waitUntil(() => saveDialogs() === 1, 'No Save run dialog', 5000);
    closeDialog('cancel');
    await waitUntil(() => saveDialogs() === 0, 'The Save run dialog stays open', 5000);

    const outside = ui.button('outside', () => {});
    document.body.append(outside);
    try {
      outside.focus();
      expect(press(outside, 's', {ctrl: true}), false, 'From outside the view the browser save would open');
    } finally {
      outside.remove();
    }
    // the focus lost to the page, as after a panel closes; the press from outside, on the same path, saved nothing
    expect(document.activeElement === document.body, true, 'The focus should be lost');
    expect(press(document.body, 's', {ctrl: true}), false, 'With the focus lost the browser save would open');
    await waitUntil(() => saveDialogs() > 0, 'No Save run dialog with the focus lost', 5000);
    expect(saveDialogs(), 1, 'The press from outside the view saved too');
    closeDialog('cancel');
  });

  test('Next moves the focus into the opened step', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true});
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'The step should run and offer Next', 10000);
    buttons(view.root, 'Next')[0].click();
    await waitUntil(() => view.root.contains(document.activeElement) && !document.activeElement!.closest('.mtl-tree'),
      `The focus should be in the opened step, is on ${document.activeElement?.tagName}`, 10000);
  });

  test('Alt+PgDn and Alt+PgUp go to the next and previous step, even before a run or with the focus lost', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed add');
    const expectChosen = (label: string) =>
      waitUntil(() => chosen(view).join() === label, () => `${label} should be chosen, is ${chosen(view)}`, 5000);
    press(row(view, 'fixed add'), 'PageDown', {alt: true});
    await expectChosen('fixed mul');
    // the list holds steps, so Next skips it
    press(document.activeElement!, 'PageDown', {alt: true});
    await expectChosen('item 1');
    (document.activeElement as HTMLElement).blur();
    press(document.body, 'PageUp', {alt: true});
    await expectChosen('fixed mul');
  });
});

category('TreeWizard: keyboard focus', () => {
  const track = closeViewsAfter();

  test('Back and Next show the keyboard focus', async () => {
    const view = await openKeyboardWorkflow(track);
    await choose(view, 'fixed mul');
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true});
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'The step should run and offer Next', 10000);
    for (const label of ['Back', 'Next']) {
      const button = buttons(view.root, label)[0];
      button.focus();
      expect(getComputedStyle(button).outlineStyle !== 'none', true, `${label} shows no focus outline`);
    }

    const dialog = ui.dialog('Outside Compute2').show();
    try {
      const button = dialog.root.querySelector('button')!;
      button.focus();
      expect(getComputedStyle(button).outlineStyle, 'none', 'Buttons outside Compute2 keep the platform style');
    } finally {
      dialog.close();
    }
  });

  test('Tab goes round the steps tree and the open step only', async () => {
    const view = await openKeyboardWorkflow(track);
    await waitUntil(() => buttons(view.root, 'Next').some((el) => el.getClientRects().length > 0),
      'The workflow page is not laid out', 5000);
    expect(tabRound(treeStop(view)), JSON.stringify(['Actions', 'Steps']), 'Workflow page');

    await choose(view, 'fixed mul');
    press(row(view, 'fixed mul'), 'Enter', {ctrl: true});
    await waitUntil(() => buttons(view.root, 'Next').length > 0, 'The step should run and offer Next', 10000);
    expect(tabRound(treeStop(view)), JSON.stringify(['Inputs', 'Steps']), 'Step');
    expectTabEnters(ribbonStops(view)[0], ['Steps', 'Inputs']);
  });

  test('The plus of a list takes the focus and opens the steps to add', async () => {
    const view = await openKeyboardWorkflow(track);
    const list = row(view, 'list');
    list.focus();
    press(list, 'Enter');
    const plusOf = () => row(view, 'list').querySelector<HTMLElement>('dg-combo-popup');
    await waitUntil(() => !!plusOf(), 'The chosen list should show its plus', 5000);
    const plus = plusOf()!;
    expect(plus.tabIndex, 0, 'The plus should take the focus');
    plus.focus();
    press(plus, 'Enter');
    await waitUntil(() => !!document.activeElement?.closest('[dock-spawn-title="Steps to add"]'),
      `The focus should be in the steps to add, is on ${document.activeElement?.className}`, 10000);
    const stepsToAdd = document.activeElement!.closest('[dock-spawn-title]')!;
    const firstAdd = stepsToAdd.querySelector<HTMLElement>('[tabindex="0"]')!;
    expect(tabRound(firstAdd), JSON.stringify(['Actions', 'Steps', 'Steps to add']), 'Tab round');
  });
});

category('RFV: keyboard shortcuts', () => {
  const track = closeViewsAfter();

  async function openSimpleInputs() {
    await awaitWebComponents();
    const call = DG.Func.byName('Compute2:SimpleInputs2').prepare();
    const view = track(await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as unknown as DG.ViewBase);
    grok.shell.addView(view);
    await waitUntil(() => view.root.querySelectorAll('input.ui-input-editor').length >= 3, 'Form not rendered', 15000);
    return view;
  }

  const inputs = (view: DG.ViewBase) =>
    Array.from(view.root.querySelectorAll<HTMLInputElement>('input.ui-input-editor'));

  test('Ctrl+Enter on invalid inputs warns, then runs once they are valid', async () => {
    const view = await openSimpleInputs();
    clearBalloons();
    press(inputs(view)[0], 'Enter', {ctrl: true});
    await waitUntil(() => balloons().some((text) => text.includes('Fix the invalid inputs')), 'No warning', 5000);

    const [a, b, c] = inputs(view);
    type(a, '1');
    type(b, '2');
    type(c, 'x');
    await waitUntil(() => buttons(view.root, 'Run')[0]?.getAttribute('aria-disabled') !== 'true',
      'Run stays disabled', 5000);
    press(c, 'Enter', {ctrl: true});
    await waitUntil(() => view.root.querySelector('[dock-spawn-title="Output"]') !== null, 'No outputs after the run',
      10000);
  });

  test('Ctrl+S warns like Save before a run and opens the save dialog after it, even with the focus lost', async () => {
    const view = await openSimpleInputs();
    clearBalloons();
    press(inputs(view)[0], 's', {ctrl: true});
    await waitUntil(() => balloons().some((text) => text.includes('Results are outdated')), 'No warning', 5000);
    expect(dialogTitles().includes('Save run'), false, 'No dialog before a run');

    const [a, b, c] = inputs(view);
    type(a, '1');
    type(b, '2');
    type(c, 'x');
    await waitUntil(() => buttons(view.root, 'Run')[0]?.getAttribute('aria-disabled') !== 'true',
      'Run stays disabled', 5000);
    press(c, 'Enter', {ctrl: true});
    await waitUntil(() => view.root.querySelector('[dock-spawn-title="Output"]') !== null, 'No outputs after the run',
      10000);
    // with the focus lost to the page, as after a panel closes
    c.blur();
    expect(press(document.body, 's', {ctrl: true}), false, 'The browser save would open');
    await waitUntil(() => dialogTitles().includes('Save run'), 'No Save run dialog', 5000);
    closeDialog('cancel');
  });

  test('The disabled Run button explains why on hover', async () => {
    const view = await openSimpleInputs();
    const run = buttons(view.root, 'Run')[0];
    expect(run.getAttribute('aria-disabled') === 'true', true, 'Run should be disabled');
    await expectHoverable(run);
    run.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    run.dispatchEvent(new MouseEvent('mouseenter'));
    await waitUntil(() => Array.from(document.querySelectorAll('.d4-tooltip'))
      .some((el) => (el.textContent ?? '').includes('invalid')), 'No tooltip with the reason', 3000);
  });

  test('Tab goes round the form only, not the ribbon or an output grid', async () => {
    await awaitWebComponents();
    const call = DG.Script.create(`//name: RfvTabGrid
//language: javascript
//input: int n = 3
//output: dataframe table {viewer: Grid()}
table = DG.DataFrame.fromCsv('x\\n1\\n2\\n3');`).prepare({n: 3});
    const view = track(await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as unknown as DG.ViewBase);
    grok.shell.addView(view);
    const visibleGrid = () => Array.from(view.root.querySelectorAll<HTMLElement>('dg-viewer canvas'))
      .find((el) => el.getClientRects().length > 0);
    await waitUntil(() => buttons(view.root, 'Run').some((el) => el.getAttribute('aria-disabled') !== 'true'),
      'Run stays disabled', 15000);
    buttons(view.root, 'Run')[0].click();
    await waitUntil(() => !!visibleGrid(), 'No output grid', 10000);

    expect(tabRound(inputs(view)[0]), JSON.stringify(['Inputs']), 'Tab round');
    expectTabEnters(visibleGrid()!, ['Inputs']);
    expectTabEnters(ribbonStops(view)[0], ['Inputs']);
  });

  test('Tab in the run details of the view stays in their fields', async () => {
    const view = await openSimpleInputs();
    const host = ui.div();
    view.root.append(host);
    const app = Vue.createApp(RunDetailsForm, {run: DG.Func.byName('Compute2:SimpleInputs2').prepare()});
    app.mount(host);
    try {
      const fields = Array.from(host.querySelectorAll<HTMLElement>('.ui-input-editor'));
      // Vue skips an event stamped in the millisecond its listener was attached
      await new Promise((resolve) => setTimeout(resolve, 2));
      fields[0].focus();
      press(fields[0], 'Tab');
      expect(document.activeElement === fields[1], true,
        `Tab from Name should go to Description, is on ${panelOf(document.activeElement)}`);
    } finally {
      app.unmount();
      host.remove();
    }
  });
});
