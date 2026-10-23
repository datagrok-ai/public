import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray, awaitCheck} from '@datagrok-libraries/test/src/test';
import {
  addToComparison, compareWithOriginal, confirmOutdated, openComparisonViews, openRunComparison,
} from '../components/RunComparison/comparison-views';
import {ComparisonAddition, ComparisonService} from '../components/RunComparison/comparison-service';
import {entryFromLiveCall} from '../components/RunComparison/entry-extraction';

const snapshot = (name: string) => entryFromLiveCall(
  DG.Script.create('//name: rcViewsProbe\n//language: javascript\n//output: int x\nx = 1;').prepare(), false, name);

// picks an option of the Add to comparison dialog and confirms it
async function pickInDialog(option: (value: string) => boolean) {
  const radios = () => Array.from(document.querySelectorAll<HTMLInputElement>('.d4-dialog input[type="radio"]'));
  await awaitCheck(() => radios().length > 0, 'Add to comparison dialog not shown', 3000);
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
      await awaitCheck(() => handle.service.entryNames.length === 1, 'start-up entries not added', 3000);
      const picker = modelPicker(handle.view.root);
      expect(picker.selected(), modelLabel(model));

      // a run that can't be loaded is reported, the rest is still added
      await handle.service.add({runIds: ['rc-run-that-does-not-exist'], entries: [snapshot('Second')], model: other});
      expectArray(handle.service.entryNames, ['First', 'Second']);
      expect(picker.selected(), modelLabel(model));
      picker.open();
      await awaitCheck(() => picker.items().includes(modelLabel(other)), 'added model not in the list', 3000);
      expect(picker.items().includes(modelLabel(model)), true);
    } finally {
      handle.view.close();
    }
  });
});

category('RunComparison: views', () => {
  test('a snapshot without outdated outputs compares without asking', async () => {
    let called = false;
    confirmOutdated(snapshot('Current'), () => called = true);
    expect(called, true);
  });

  test('adding opens a view when none is open, then adds to the picked one or a new one', async () => {
    openComparisonViews.value.forEach((item) => item.view.close());
    await awaitCheck(() => openComparisonViews.value.length === 0, 'views from earlier tests are still open', 3000);
    try {
      addToComparison(snapshot('First'));
      expect(openComparisonViews.value.length, 1);
      const [first] = openComparisonViews.value;
      await awaitCheck(() => first.service.entryNames.length === 1, 'the snapshot did not reach the new view', 3000);

      addToComparison(snapshot('Second'));
      await pickInDialog((value) => value.startsWith('1. '));
      await awaitCheck(() => first.service.entryNames.length === 2, 'the snapshot did not reach the picked view', 3000);
      expectArray(first.service.entryNames, ['First', 'Second']);

      addToComparison(snapshot('Third'));
      await pickInDialog((value) => value === 'New comparison');
      expect(openComparisonViews.value.length, 2);
      const second = openComparisonViews.value[1];
      await awaitCheck(() => second.service.entryNames.length === 1, 'the snapshot did not reach the new view', 3000);
      expectArray(second.service.entryNames, ['Third']);
      expectArray(first.service.entryNames, ['First', 'Second']);

      // a run that can't be loaded is reported, the snapshot still arrives
      compareWithOriginal(snapshot('Edited'), 'rc-run-that-does-not-exist');
      expect(openComparisonViews.value.length, 3);
      const third = openComparisonViews.value[2];
      await awaitCheck(() => third.service.entryNames.length === 1,
        'the snapshot did not reach the original view', 3000);
      expectArray(third.service.entryNames, ['Edited']);
    } finally {
      document.querySelectorAll('.d4-dialog').forEach((el) => el.remove());
      openComparisonViews.value.forEach((item) => item.view.close());
    }
    await awaitCheck(() => openComparisonViews.value.length === 0, 'closed views are still tracked', 3000);
  });
});
