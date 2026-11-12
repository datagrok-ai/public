import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  CustomFunctionView, RunStorageHost,
} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {attachRunStorage} from '../components/RunStorage/custom-function-view-runs';
import {DAY_CALCULATION, RunStorageFixture, spaceRuns, withRunStorage} from './run-storage';
import {RunsChange, listMyModelRuns, runsChanged$} from '../components/RunStorage/run-storage';
import {closeView} from './utils';

const TEST_CUSTOM_VIEW = 'Compute2:TestCustomView';

class CountingView extends CustomFunctionView {
  counts = {beforeSave: 0, afterSave: 0, afterLoad: 0};

  buildIO() {
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

async function readyView(f: RunStorageFixture) {
  const view = new CountingView(DAY_CALCULATION);
  await view.isReady.pipe(filter((ready) => ready), take(1)).toPromise();
  view.linkFunccall(await f.dayRun());
  return view;
}

type SaveResult = Awaited<ReturnType<RunStorageHost['saveDialog']>>;

/** A host whose dialogs resolve fixed answers right away and record what they were asked. */
const fixedHost = (save: SaveResult, open: DG.FuncCall | null = null) => {
  const asked = {save: [] as DG.FuncCall[], open: [] as DG.Func[]};
  const host: RunStorageHost = {
    saveDialog: async (call) => {
      asked.save.push(call);
      return save;
    },
    openDialog: async (func) => {
      asked.open.push(func);
      return open;
    },
  };
  return {host, asked};
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

category('Custom function view', () => {
  test('Save asks the host, then saves the run with its details', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const name = `CFV run ${Date.now()}`;
    const {host, asked} = fixedHost({metadata: {name, description: 'From the view', tags: ['cfv']}, target: {}});
    view.runStorage = host;
    const saved = await view.saveRunWithDialog(view.funcCall!);
    f.track(saved!.id);

    expect(asked.save.length, 1);
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
    view.runStorage = fixedHost({metadata: {name: `CFV space run ${Date.now()}`}, target: {spaceId: space.id}}).host;
    const saved = await view.saveRunWithDialog(view.funcCall!);
    f.track(saved!.id);

    const linked = await spaceRuns(space, DG.Func.byName(DAY_CALCULATION));
    expect(JSON.stringify(linked.map((run) => run.id)), JSON.stringify([saved!.id]));
  }));

  test('A save to a space tells open run browsers', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const space = await f.space();
    view.runStorage = fixedHost({metadata: {name: `CFV notified ${Date.now()}`}, target: {spaceId: space.id}}).host;
    const changes: RunsChange[] = [];
    const sub = runsChanged$.subscribe((change) => changes.push(change));
    try {
      const saved = await view.saveRunWithDialog(view.funcCall!);
      f.track(saved!.id);
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
      view.runStorage = fixedHost({metadata: {name: `CFV unlinked ${Date.now()}`}, target: {spaceId: space.id}}).host;
      const saved = await view.saveRunWithDialog(view.funcCall!);
      f.track(saved!.id);

      expect(view.funcCall!.id, saved!.id);
      expect((await historyUtils.loadRun(saved!.id)).id, saved!.id);
    }));

  test('A cancelled Save saves nothing', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    view.runStorage = fixedHost(null).host;
    const before = view.funcCall;
    expect(await view.saveRunWithDialog(view.funcCall!), null);
    expect(view.counts.beforeSave + view.counts.afterSave, 0);
    expect(view.funcCall === before, true);
  }));

  test('Open loads the run picked in the host dialog', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(`CFV open ${Date.now()}`);
    const view = await readyView(f);
    const {host, asked} = fixedHost(null, await historyUtils.loadRun(saved.id));
    view.runStorage = host;
    const opened = await view.openRunWithDialog();

    expect(asked.open.length, 1);
    expect(asked.open[0].nqName, DG.Func.byName(DAY_CALCULATION).nqName);
    expect(opened!.id, saved.id);
    expect(view.funcCall!.id, saved.id);
    expect(view.funcCall!.inputs['solutionAdded'], 0.5);
    expect(view.counts.afterLoad, 1);
  }));

  test('A cancelled Open keeps the current run', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    view.runStorage = fixedHost(null, null).host;
    const before = view.funcCall;
    expect(await view.openRunWithDialog(), null);
    expect(view.funcCall === before, true);
    expect(view.counts.afterLoad, 0);
  }));

  test('An Open that fails keeps the current run', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    view.runStorage = {...fixedHost(null).host, openDialog: () => historyUtils.loadRun(missingRunId())};
    const before = view.funcCall;
    expect(await view.openRunWithDialog(), null);
    expect(view.funcCall === before, true);
    expect(view.counts.afterLoad, 0);
  }));

  test('Without a host there are no History and Save buttons', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    expect(hasIcon(view, 'history') || hasIcon(view, 'save'), false);
    expect(await view.saveRunWithDialog(view.funcCall!), null);
    expect(await view.openRunWithDialog(), null);
    expect(view.counts.beforeSave + view.counts.afterLoad, 0);
    view.runStorage = fixedHost(null).host;
    expect(hasIcon(view, 'history') && hasIcon(view, 'save'), true);
  }));

  test('Load by id links the run and calls the load hook once', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(`CFV load ${Date.now()}`);
    const view = await readyView(f);
    await view.loadRun(saved.id);
    expect(view.counts.afterLoad, 1);
    expect(view.funcCall!.id, saved.id);
  }));

  test('The editor setup gives the view both dialogs and takes them back', () => withRunStorage(async (f) => {
    const view = await readyView(f);
    const detach = attachRunStorage(view);
    expect(typeof view.runStorage?.saveDialog, 'function');
    expect(typeof view.runStorage?.openDialog, 'function');
    expect(hasIcon(view, 'history') && hasIcon(view, 'save'), true);
    detach();
    expect(view.runStorage == null, true);
    expect(hasIcon(view, 'history') || hasIcon(view, 'save'), false);
  }));

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
