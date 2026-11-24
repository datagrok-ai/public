import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {CustomFunctionView} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import {historyInput} from '@datagrok-libraries/compute-utils/function-views/src/history-input';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {DAY_CALCULATION, RunStorageFixture, spaceRuns, withRunStorage} from './run-storage';
import {RunsChange, listMyModelRuns, runsChanged$} from '../components/RunStorage/run-storage';
import {closeView, deleteSavedRun, findButton} from './utils';

const TEST_CUSTOM_VIEW = 'Compute2:TestCustomView';

class CountingView extends CustomFunctionView {
  counts = {beforeSave: 0, afterSave: 0, afterLoad: 0};

  buildIO(): HTMLElement {
    return ui.div();
  }

  async onBeforeSaveRun() {
    this.counts.beforeSave++;
  }

  async onAfterSaveRun() {
    this.counts.afterSave++;
  }

  async onAfterLoadRun() {
    this.counts.afterLoad++;
  }
}

class ContentView extends CountingView {
  content?: HTMLElement;

  constructor(private makeContent: () => HTMLElement) {
    super(DAY_CALCULATION);
  }

  buildIO() {
    this.content = this.makeContent();
    return this.content;
  }
}

async function readyView(f: RunStorageFixture) {
  const view = new CountingView(DAY_CALCULATION);
  await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
  view.linkFunccall(await f.dayRun());
  return view;
}

const dialogsNow = (selector: string) => new Set(document.querySelectorAll(selector));

/** The Compute2 dialog matching the selector that opened since `before`. */
async function openedDialog(selector: string, before: Set<Element>): Promise<HTMLElement> {
  let found: HTMLElement | undefined;
  await awaitCheck(() => {
    found = [...document.querySelectorAll<HTMLElement>(selector)].find((el) => !before.has(el));
    return !!found;
  }, `${selector} should open`, 10000);
  return found!;
}

/** Selects the run named `name` in the open run browser. */
async function selectListedRun(browser: HTMLElement, name: string) {
  const row = () => [...browser.querySelectorAll<HTMLElement>('.c2-runs-row')]
    .find((el) => el.querySelector('.c2-runs-row-name')?.textContent === name);
  await awaitCheck(() => !!row(), 'the saved run should be listed', 15000);
  row()!.click();
}

/** Opens the run named `name` in the run browser that opens after `before`. */
async function openInBrowser(name: string, before: Set<Element>) {
  const browser = await openedDialog('.c2-runs-browser', before);
  await selectListedRun(browser, name);
  findButton(browser, 'open').click();
}

const typeInto = (field: HTMLInputElement | HTMLTextAreaElement, text: string) => {
  field.value = text;
  field.dispatchEvent(new Event('input'));
};

// the id of a run that was never saved
const missingRunId = () => {
  const call = DG.Func.byName(DAY_CALCULATION).prepare({});
  call.newId();
  return call.id;
};

const ribbonIcons = (view: CustomFunctionView) => view.buildRibbonPanels().flat();
const hasIcon = (view: CustomFunctionView, name: string) =>
  ribbonIcons(view).some((icon) => icon.classList.contains(`fa-${name}`));
const exportIcon = (view: CustomFunctionView) =>
  ribbonIcons(view).find((icon) => icon.classList.contains('fa-arrow-to-bottom'));

async function newView() {
  const view = new CountingView(DAY_CALCULATION);
  await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
  return view;
}

const openEditor = (): Promise<CustomFunctionView> =>
  grok.functions.call('Compute2:CustomFunctionViewEditor', {call: DG.Func.byName(TEST_CUSTOM_VIEW).prepare({})});

/** Opens the custom view editor as the page does when its link has the run id. */
async function openFromLink(runId: string) {
  const handled = globalThis.initialURLHandled;
  // the start URL can't be set otherwise; the own property hides the prototype getter until deleted
  Object.defineProperty(grok.shell, 'startUri', {get: () => `${location.origin}/?id=${runId}`, configurable: true});
  globalThis.initialURLHandled = false;
  try {
    return await openEditor();
  } finally {
    delete (grok.shell as any).startUri;
    globalThis.initialURLHandled = handled;
  }
}

/** Records the names of the files downloaded during the action instead of downloading them. */
async function withDownloads(action: (downloads: string[]) => Promise<void>) {
  const download = DG.Utils.download;
  const downloads: string[] = [];
  DG.Utils.download = (fileName: string) => {
    downloads.push(fileName);
  };
  try {
    await action(downloads);
  } finally {
    DG.Utils.download = download;
  }
}

category('Custom function view', () => {
  test('Save asks in the Save dialog, then saves the run with its details', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const name = `CFV run ${Date.now()}`;
    const before = dialogsNow('.c2-runs-dialog');
    const saving = view.saveRunWithDialog(view.funcCall!);
    const dialog = await openedDialog('.c2-runs-dialog', before);
    const [nameField, descriptionField] = dialog.querySelectorAll<HTMLInputElement>('.c2-runs-form .ui-input-editor');
    typeInto(nameField, name);
    typeInto(descriptionField, 'From the view');
    const saveButton = dialog.querySelector<HTMLElement>('.ui-btn-ok')!;
    await awaitCheck(() => saveButton.getAttribute('aria-disabled') !== 'true', 'the save should be allowed', 10000);
    saveButton.click();
    const saved = await saving;
    f.track(saved!.id);

    expect(view.counts.beforeSave, 1);
    expect(view.counts.afterSave, 1);
    expect(view.funcCall!.id, saved!.id);
    const loaded = await historyUtils.loadRun(saved!.id);
    expect(getRunName(loaded), name);
    expect(loaded.options['description'], 'From the view');
  }));

  test('Save into a space links the run there', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const space = await f.space();
    const saved = await view.saveRun(view.funcCall!, {name: `CFV space run ${Date.now()}`}, {spaceId: space.id});
    f.track(saved.id);

    const linked = await spaceRuns(space, DG.Func.byName(DAY_CALCULATION));
    expect(JSON.stringify(linked.map((run) => run.id)), JSON.stringify([saved.id]));
  }));

  test('A save to a space tells open run browsers', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const space = await f.space();
    const changes: RunsChange[] = [];
    const sub = runsChanged$.subscribe((change) => changes.push(change));
    try {
      const saved = await view.saveRun(view.funcCall!, {name: `CFV notified ${Date.now()}`}, {spaceId: space.id});
      f.track(saved.id);
      await awaitCheck(() => changes.length > 0, 'the save should reach run browsers', 3000);
      expect(JSON.stringify(changes), JSON.stringify([{spaceIds: [space.id]}]));
    } finally {
      sub.unsubscribe();
    }
  }));

  test('A save that fails only to link into the space keeps the run', () =>
    withRunStorage(async (f) => {
      const view = await readyView(f);
      const space = await f.space();
      await grok.dapi.spaces.delete(space);
      const saved = await view.saveRun(view.funcCall!, {name: `CFV unlinked ${Date.now()}`}, {spaceId: space.id});
      f.track(saved.id);

      expect(view.funcCall!.id, saved.id);
      expect((await historyUtils.loadRun(saved.id)).id, saved.id);
    }));

  test('A cancelled Save saves nothing', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const before = view.funcCall;
    const dialogs = dialogsNow('.c2-runs-dialog');
    const saving = view.saveRunWithDialog(view.funcCall!);
    findButton(await openedDialog('.c2-runs-dialog', dialogs), 'cancel').click();
    expect(await saving, null);
    expect(view.counts.beforeSave + view.counts.afterSave, 0);
    expect(view.funcCall === before, true);
  }));

  test('Open loads the run picked in the run browser', () => withRunStorage(async (f) => {
    const name = `CFV open ${Date.now()}`;
    const saved = await f.saveDayRun(name);
    const view = await readyView(f);
    const browsers = dialogsNow('.c2-runs-browser');
    const opening = view.openRunWithDialog();
    await openInBrowser(name, browsers);
    const opened = await opening;

    expect(opened!.id, saved.id);
    expect(view.funcCall!.id, saved.id);
    expect(view.funcCall!.inputs['solutionAdded'], 0.5);
    expect(view.counts.afterLoad, 1);
  }));

  test('Open and Save dialogs list runs without the not-completed mark', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const name = `CFV no mark ${Date.now()}`;
    f.track((await view.saveRun(view.funcCall!, {name})).id);
    const listedRow = async (dialog: HTMLElement) => {
      let row: HTMLElement | undefined;
      await awaitCheck(() => !!(row = [...dialog.querySelectorAll<HTMLElement>('.c2-runs-row')]
        .find((el) => el.querySelector('.c2-runs-row-name')?.textContent === name)),
      'the saved run should be listed', 15000);
      return row!;
    };

    const browsers = dialogsNow('.c2-runs-browser');
    const opening = view.openRunWithDialog();
    const browser = await openedDialog('.c2-runs-browser', browsers);
    const browserMark = (await listedRow(browser)).querySelector('.c2-runs-not-completed');
    findButton(browser, 'cancel').click();
    await opening;

    const dialogs = dialogsNow('.c2-runs-dialog');
    const saving = view.saveRunWithDialog(view.funcCall!);
    const saveDialog = await openedDialog('.c2-runs-dialog', dialogs);
    const saveMark = (await listedRow(saveDialog)).querySelector('.c2-runs-not-completed');
    findButton(saveDialog, 'cancel').click();
    await saving;

    expect(browserMark, null);
    expect(saveMark, null);
  }));

  test('A cancelled Open keeps the current run', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const before = view.funcCall;
    const browsers = dialogsNow('.c2-runs-browser');
    const opening = view.openRunWithDialog();
    findButton(await openedDialog('.c2-runs-browser', browsers), 'cancel').click();
    expect(await opening, null);
    expect(view.funcCall === before, true);
    expect(view.counts.afterLoad, 0);
  }));

  test('An Open that fails keeps the current run', () => withRunStorage(async (f) => {
    const name = `CFV open fails ${Date.now()}`;
    const saved = await f.saveDayRun(name);
    const view = await readyView(f);
    const before = view.funcCall;
    const browsers = dialogsNow('.c2-runs-browser');
    const opening = view.openRunWithDialog();
    const browser = await openedDialog('.c2-runs-browser', browsers);
    await selectListedRun(browser, name);
    await deleteSavedRun(saved.id);
    findButton(browser, 'open').click();
    expect(await opening, null);
    expect(view.funcCall === before, true);
    expect(view.counts.afterLoad, 0);
  }));

  test('With Compute2 the view has History and Save buttons', async () => {
    const view = await newView();
    expect(hasIcon(view, 'history') && hasIcon(view, 'save'), true);
  });

  test('Load by id links the run and calls the load hook once', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(`CFV load ${Date.now()}`);
    const view = await readyView(f);
    await view.loadRun(saved.id);
    expect(view.counts.afterLoad, 1);
    expect(view.funcCall!.id, saved.id);
  }));

  test('A new view starts with the default inputs', async () => {
    const view = new CountingView(DAY_CALCULATION);
    await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
    expect(view.funcCall!.inputs['incomingVolume'], 1);
    expect(view.funcCall!.inputs['maximumGrowthRate'], 0.4);
  });

  test('A deep link opens its run in the editor', () => withRunStorage(async (f) => {
    const view = await openEditor();
    view.funcCall!.inputs['a'] = 2;
    view.funcCall!.inputs['b'] = 3;
    await view.funcCall!.call();
    const saved = await view.saveRun(view.funcCall!, {name: `CFV deep link ${Date.now()}`});
    f.track(saved.id);

    const opened = await openFromLink(saved.id);
    expect(opened.funcCall!.id, saved.id);
    expect(opened.funcCall!.inputs['a'], 2);
  }));

  test('A run saved from a custom view is a Compute2 run that previews in the view', () => withRunStorage(async (f) => {
    const view = await openEditor();
    view.funcCall!.inputs['a'] = 2;
    view.funcCall!.inputs['b'] = 3;
    await view.funcCall!.call();
    const saved = await view.saveRun(view.funcCall!, {name: `CFV handler ${Date.now()}`});
    f.track(saved.id);

    const run = await grok.dapi.functions.calls.allPackageVersions().include('session.user,func.package').find(saved.id);
    const handler = DG.ObjectHandler.forEntity(run);
    expect(handler?.type, 'ComputeRun');
    const preview = await handler!.renderPreview(run) as unknown as CustomFunctionView;
    expect(preview.funcCall!.id, saved.id);
    expect(preview.funcCall!.inputs['a'], 2);
  }));

  test('A deep link to a run that fails to load opens a new call', async () => {
    const missing = missingRunId();
    const view = await openFromLink(missing);
    expect(!!view.funcCall && view.funcCall.id !== missing, true);
  });

  test('The Export icon downloads the first export', () => withDownloads(async (downloads) => {
    const view = await newView();
    view.exports = [
      {name: 'First', export: async () => DG.FileInfo.fromString('first.txt', '1')},
      {name: 'Second', export: async () => DG.FileInfo.fromString('second.txt', '2')},
    ];
    exportIcon(view)!.click();
    await awaitCheck(() => downloads.length > 0, 'the export should download', 3000);
    expect(JSON.stringify(downloads), JSON.stringify(['first.txt']));
  }));

  test('Export turned off greys the icon and downloads nothing', () => withDownloads(async (downloads) => {
    const view = await newView();
    view.exports = [{name: 'First', export: async () => DG.FileInfo.fromString('first.txt', '1')}];
    view.exportEnabled = false;
    const icon = exportIcon(view)!;
    expect(icon.classList.contains('d4-disabled'), true);
    icon.click();
    expect(downloads.length, 0);
    view.exportEnabled = true;
    expect(exportIcon(view)!.classList.contains('d4-disabled'), false);
  }));

  const contents = [['div', () => ui.div()], ['box', () => ui.splitH([ui.box(), ui.box()])]] as const;
  for (const [kind, makeContent] of contents) {
    test(`A ${kind} content fills the view`, async () => {
      const view = new ContentView(makeContent);
      await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
      grok.shell.addView(view);
      try {
        const root = view.root.getBoundingClientRect();
        const content = view.content!.getBoundingClientRect();
        expect(Math.round(content.width), Math.round(root.width));
        expect(Math.round(content.height), Math.round(root.height));
      } finally {
        closeView(view);
      }
    });
  }

  test('Without exports there is no Export icon', async () => {
    const view = await newView();
    expect(exportIcon(view) == null, true);
  });

  test('The editor opens even when its first run fails to load', async () => {
    const editorCall = DG.Func.byName(TEST_CUSTOM_VIEW).prepare({});
    const missing = missingRunId();
    editorCall.aux.initialRunId = missing;
    const view: CustomFunctionView = await grok.functions.call('Compute2:CustomFunctionViewEditor',
      {call: editorCall});
    expect(!!view.funcCall && view.funcCall.id !== missing, true);
  });

  test('A run saved from a custom view is listed under it and opens in it again', () => withRunStorage(async (f) => {
    const editorCall = DG.Func.byName(TEST_CUSTOM_VIEW).prepare({});
    const view: CustomFunctionView = await grok.functions.call('Compute2:CustomFunctionViewEditor',
      {call: editorCall});
    const space = await f.space();
    let opened: DG.ViewBase | undefined;
    try {
      view.funcCall!.inputs['a'] = 2;
      view.funcCall!.inputs['b'] = 3;
      await view.funcCall!.call();
      const saved = await view.saveRun(view.funcCall!, {name: `CFV reopen ${Date.now()}`}, {spaceId: space.id});
      f.track(saved.id);
      const customView = DG.Func.byName(TEST_CUSTOM_VIEW);
      expect((await listMyModelRuns(customView)).some((run) => run.id === saved.id), true);
      expect((await spaceRuns(space, customView)).some((run) => run.id === saved.id), true);

      await grok.functions.call('Compute2:OpenWorkflowRun', {id: saved.id});
      const isReopened = (v: DG.ViewBase) => (v as CustomFunctionView).funcCall?.id === saved.id;
      await awaitCheck(() => [...grok.shell.views].some(isReopened), 'the run should open in its custom view', 15000);
      opened = [...grok.shell.views].find(isReopened);
      expect((opened as CustomFunctionView).funcCall!.inputs['a'], 2);
    } finally {
      if (opened)
        closeView(opened);
    }
  }));
});

const icon = (input: DG.InputBase, name: string) => input.root.querySelector<HTMLElement>(`.fa-${name}`)!;

category('Custom function view: history input', () => {
  test('Picks a run in the run browser', () => withRunStorage(async (f) => {
    const name = `History input ${Date.now()}`;
    const saved = await f.save(await f.dayRun({called: true}), {metadata: {name}});
    const changes: (DG.FuncCall | null)[] = [];
    const input = historyInput('Previous run', {
      model: DAY_CALCULATION, nullable: true, onValueChanged: (run) => changes.push(run),
    });
    expect(input.value, null);
    expect(icon(input, 'times').style.display, 'none');

    const browsers = dialogsNow('.c2-runs-browser');
    icon(input, 'search').click();
    await openInBrowser(name, browsers);
    await awaitCheck(() => input.value !== null, 'the opened run should be the value', 15000);
    expect(input.value!.id, saved.id);
    expect(input.value!.outputs['dailyProfile'] instanceof DG.DataFrame, true, 'the run comes with outputs');
    expect(input.stringValue.startsWith(name), true);
    expect(changes.length, 1);

    icon(input, 'times').click();
    expect(input.value, null);
    expect(input.stringValue, '');
    expect(changes.length, 2);
  }));

  test('A cancelled pick keeps the value', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(`History input cancel ${Date.now()}`);
    const run = await historyUtils.loadRun(saved.id);
    const input = historyInput('Previous run', {model: DAY_CALCULATION, value: run});
    const browsers = dialogsNow('.c2-runs-browser');
    icon(input, 'search').click();
    findButton(await openedDialog('.c2-runs-browser', browsers), 'cancel').click();
    await awaitCheck(() => dialogsNow('.c2-runs-browser').size === browsers.size, 'the browser should close', 5000);
    expect(input.value === run, true);
  }));
});
