import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray, awaitCheck, waitUntil} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {OUTPUT_OUTDATED_PATH} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  addRunToComparison, addToComparison, compareWithOriginal, confirmOutdated, openComparisonViews, openRunComparison,
} from '../components/RunComparison/comparison-views';
import {ComparisonAddition, ComparisonService} from '../components/RunComparison/comparison-service';
import {entryFromLiveCall} from '../components/RunComparison/entry-extraction';
import {deleteSavedRun, deleteSavedWorkflow} from './utils';
import {DAY_CALCULATION, withRunStorage} from './run-storage';
import {runWorkflow, saveWorkflowRun} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';

const snapshot = (name: string) => entryFromLiveCall(
  DG.Script.create('//name: rcViewsProbe\n//language: javascript\n//output: int x\nx = 1;').prepare(),
  {isOutputOutdated: false}, name);

// picks an option of the Add to comparison dialog and confirms it
async function pickInDialog(option: (value: string) => boolean) {
  const radios = () => Array.from(document.querySelectorAll<HTMLInputElement>('.d4-dialog input[type="radio"]'));
  await waitUntil(() => radios().length > 0, 'Add to comparison dialog not shown', 3000);
  radios().find((radio) => option(radio.dataset.value ?? ''))!.click();
  Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog button'))
    .find((b) => b.textContent?.trim() === 'OK')!.click();
}

// the model picker's header, and its list once opened
const modelPicker = (root: HTMLElement) => {
  const header = Array.from(root.querySelectorAll<HTMLElement>('div'))
    .find((div) => div.lastElementChild?.textContent === '\u25BE')!;
  return {
    selected: () => header.firstElementChild?.textContent,
    open: () => header.click(),
    items: () => Array.from(root.querySelector('input[placeholder="Filter models..."]')?.parentElement?.children ?? [])
      .filter((el) => el.tagName === 'DIV').map((el) => el.textContent),
  };
};

const modelLabel = (f: DG.Func) => f.friendlyName ?? f.name;

category('RunComparison: service', () => {
  test('the service passes additions and entry names to the connected view', async () => {
    const service = new ComparisonService();
    const received: ComparisonAddition[] = [];
    service.connect({add: async (addition) => {
      received.push(addition);
    }, entryNames: () => ['A', 'B']});
    await service.add({runIds: ['r1'], selectModel: true});
    expect(received.length, 1);
    expectArray(received[0].runIds ?? [], ['r1']);
    expect(received[0].selectModel, true);
    expectArray(service.entryNames, ['A', 'B']);
  });

  test('a view takes runs, entries and models through its service', async () => {
    const model = DG.Func.byName('Compute2:SimpleInputs2');
    const other = DG.Func.byName('Compute2:ObjectCooling2');
    const handle = openRunComparison({entries: [snapshot('First')], model});
    try {
      await waitUntil(() => handle.service.entryNames.length === 1, 'start-up entries not added', 3000);
      const picker = modelPicker(handle.view.root);
      expect(picker.selected(), modelLabel(model));

      // a run that can't be loaded is reported, the rest is still added
      await handle.service.add({runIds: ['rc-run-that-does-not-exist'], entries: [snapshot('Second')], model: other});
      expectArray(handle.service.entryNames, ['First', 'Second']);
      expect(picker.selected(), modelLabel(model));
      picker.open();
      await waitUntil(() => picker.items().includes(modelLabel(other)), 'added model not in the list', 3000);
      expect(picker.items().includes(modelLabel(model)), true);
    } finally {
      handle.view.close();
    }
  });
});

category('RunComparison: views', () => {
  test('the run picker lists runs of all models and adds a selection of several models', () =>
    withRunStorage(async (f) => {
      const stamp = Date.now();
      await f.saveDayRun(`RC day ${stamp}`);
      const result = await runWorkflow('Bioreactorworkflowdemo:BioreactorWorkflow', '1.0');
      const workflow = await saveWorkflowRun('Bioreactorworkflowdemo:BioreactorWorkflow', '1.0', result.state,
        {title: `RC workflow ${stamp}`});
      const handle = openRunComparison({model: DG.Func.byName(DAY_CALCULATION)});
      const root = handle.view.root;
      try {
        const picker = modelPicker(root);
        picker.open();
        await awaitCheck(() => picker.items().includes('All models'), 'no All models option', 3000);
        Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-model-dropdown .c2-runs-row'))
          .find((row) => row.textContent === 'All models')!.click();
        const row = (name: string) => root.querySelector<HTMLElement>(`.c2-runs-row-name[title="${name}"]`)
          ?.closest<HTMLElement>('.c2-runs-row');
        await awaitCheck(() => !!row(`RC day ${stamp}`) && !!row(`RC workflow ${stamp}`),
          'runs of both models should be listed', 15000);
        row(`RC day ${stamp}`)!.dispatchEvent(new MouseEvent('click', {bubbles: true, detail: 1}));
        row(`RC workflow ${stamp}`)!.dispatchEvent(new MouseEvent('click', {bubbles: true, detail: 1, ctrlKey: true}));
        const add = () => Array.from(root.querySelectorAll<HTMLButtonElement>('button'))
          .find((button) => button.textContent?.startsWith('Add selected runs'));
        await awaitCheck(() => add()?.textContent === 'Add selected runs (2)', 'both runs should be selected', 3000);
        add()!.click();
        await awaitCheck(() => handle.service.entryNames.length === 2, 'both runs should be compared', 15000);
        expectArray([...handle.service.entryNames].sort(), [`RC day ${stamp}`, `RC workflow ${stamp}`]);
      } finally {
        handle.view.close();
        await deleteSavedWorkflow(workflow.id);
      }
    }));

  test('a snapshot without outdated outputs compares without asking', async () => {
    let called = false;
    confirmOutdated(snapshot('Current'), () => called = true);
    expect(called, true);
  });

  test('adding opens a view when none is open, then adds to the picked one or a new one', async () => {
    openComparisonViews.value.forEach((item) => item.view.close());
    await waitUntil(() => openComparisonViews.value.length === 0, 'views from earlier tests are still open', 3000);
    try {
      addToComparison(snapshot('First'));
      expect(openComparisonViews.value.length, 1);
      const [first] = openComparisonViews.value;
      await waitUntil(() => first.service.entryNames.length === 1, 'the snapshot did not reach the new view', 3000);

      addToComparison(snapshot('Second'));
      await pickInDialog((value) => value.startsWith('1. '));
      await waitUntil(() => first.service.entryNames.length === 2, 'the snapshot did not reach the picked view', 3000);
      expectArray(first.service.entryNames, ['First', 'Second']);

      addToComparison(snapshot('Third'));
      await pickInDialog((value) => value === 'New comparison');
      expect(openComparisonViews.value.length, 2);
      const second = openComparisonViews.value[1];
      await waitUntil(() => second.service.entryNames.length === 1, 'the snapshot did not reach the new view', 3000);
      expectArray(second.service.entryNames, ['Third']);
      expectArray(first.service.entryNames, ['First', 'Second']);

      // a run that can't be loaded is reported, the snapshot still arrives
      compareWithOriginal(snapshot('Edited'), 'rc-run-that-does-not-exist');
      expect(openComparisonViews.value.length, 3);
      const third = openComparisonViews.value[2];
      await waitUntil(() => third.service.entryNames.length === 1,
        'the snapshot did not reach the original view', 3000);
      expectArray(third.service.entryNames, ['Edited']);
    } finally {
      document.querySelectorAll('.d4-dialog').forEach((el) => el.remove());
      openComparisonViews.value.forEach((item) => item.view.close());
    }
    await waitUntil(() => openComparisonViews.value.length === 0, 'closed views are still tracked', 3000);
  });

  test('a saved run is added without asking, once per view', async () => {
    openComparisonViews.value.forEach((item) => item.view.close());
    await awaitCheck(() => openComparisonViews.value.length === 0, 'views from earlier tests are still open', 3000);
    const call = DG.Func.byName('Compute2:SimpleInputs2').prepare({a: 1, b: 2, c: 'x'});
    call.options['title'] = 'rcViewsSaved';
    call.options[OUTPUT_OUTDATED_PATH] = serialize(true);
    const saved = await historyUtils.saveRun(call, {newId: true});
    try {
      addRunToComparison(saved.id, 'rcViewsSaved');
      expect(openComparisonViews.value.length, 1);
      const [view] = openComparisonViews.value;
      await awaitCheck(() => view.service.entryNames.length === 1, 'the saved run did not reach the new view', 3000);
      expectArray(view.service.entryNames, ['rcViewsSaved']);

      await view.service.add({runIds: [saved.id]});
      expectArray(view.service.entryNames, ['rcViewsSaved']);
    } finally {
      openComparisonViews.value.forEach((item) => item.view.close());
      await deleteSavedRun(saved.id);
    }
  });
});
