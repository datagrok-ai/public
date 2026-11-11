import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {unzipSync} from 'fflate';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {HeadlessResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {loadWorkflowConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {
  OUTPUT_OUTDATED_PATH, readRunCompleted,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {
  DEFAULT_EXPORT, callExports, exportSavedCall, exportSavedWorkflow, workflowExports,
} from '../run-exports';
import {CUSTOM_VIEW_OPTION} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ExportRun, exportRuns, modelExports, modelKind, zipExports} from '../components/RunStorage/bulk-export';
import {SpaceNode, makeIndex, unlistedSpaces} from '../components/RunStorage/space-index';
import {BulkExportDialog} from '../components/RunStorage/BulkExportDialog';
import {RunBrowser} from '../components/RunStorage/RunBrowser';
import {uniqueName, withRunStorage} from './run-storage';
import {deleteSavedWorkflow} from './utils';

const EXPORT_MODEL = 'Compute2:TestCustomExportModel';
const CUSTOM_VIEW = 'Compute2:TestCustomView';

function mount(component: Vue.Component, props: Record<string, any>) {
  const root = ui.div();
  document.body.append(root);
  const app = Vue.createApp(component, props);
  app.mount(root);
  return {root, unmount: () => {
    app.unmount();
    root.remove();
  }};
}

const text = (file: DG.FileInfo) => new TextDecoder().decode(file.data);

const callOf = (nqName: string, params: Record<string, any>, outdated: boolean) => {
  const call = DG.Func.byName(nqName).prepare(params);
  call.options[OUTPUT_OUTDATED_PATH] = String(outdated);
  return call;
};

async function saveWorkflow(nqName: string): Promise<DG.FuncCall> {
  const run = async (func: string, params: Record<string, any>) => {
    const call = DG.Func.byName(func).prepare(params);
    await call.call();
    return call.getOutputParamValue();
  };
  const result: HeadlessResult = await run('Compute2:RunWorkflow', {nqName, version: '1.0'});
  return run('Compute2:SaveWorkflowRun', {nqName, version: '1.0', state: result.state, metadata: {title: 'bulk'}});
}

category('Run storage: bulk export', () => {
  test('Export lists follow the menus, the default first', async () => {
    const model = DG.Func.byName(EXPORT_MODEL);
    expectDeepEqual(callExports(model).map(({id}) => id), [DEFAULT_EXPORT, 'rec', 'none', 'hook']);
    expectDeepEqual(callExports(model, true).map(({id}) => id), ['rec', 'none', 'hook']);
    expectDeepEqual(workflowExports({disableDefaultExport: true,
      customExports: [{id: 'a', friendlyName: 'A', handler: async () => DG.FileInfo.fromString('a.txt', '')}]}),
    [{id: 'a', name: 'A'}]);

    const wizard = DG.Func.byName('Compute2:MockTreeWizardPipeline');
    expect(modelKind(wizard, false), 'workflow');
    expect(modelKind(wizard, true), 'call');
    expectDeepEqual((await modelExports(wizard, 'workflow')).map(({name}) => name), ['Default Excel', 'Record export']);
    // a compact view exports its step
    expectDeepEqual(await modelExports(DG.Func.byName('Compute2:MockSingleStepPipeline'), 'workflow'),
      callExports(DG.Func.byName('Compute2:ObjectCooling2')));
  });

  test('The spaces that couldn\'t be listed are found below a space, each once', async () => {
    const node = (id: string, children: string[], error?: string): [string, SpaceNode] =>
      [id, {space: {id} as DG.Project, children, runs: [], error}];
    const index = makeIndex(['a', 'd'], new Map([
      node('a', ['b', 'c']), node('b', ['c'], 'no access'), node('c', [], 'timeout'), node('d', [], 'other'),
    ]));
    expectDeepEqual(unlistedSpaces(index, 'a').map(({space}) => space.id), ['b', 'c']);
    expectDeepEqual(unlistedSpaces(index, 'c').map(({space}) => space.id), ['c']);
  });

  test('The zip keeps the folders and numbers a name taken in a folder', async () => {
    const file = (name: string) => ({folders: [], file: DG.FileInfo.fromString(name, name)});
    const zip = unzipSync(zipExports([
      file('a.xlsx'), file('A.xlsx'), file('a.xlsx'), file('b'),
      {folders: ['Sub: 1', 'Sub 2'], file: DG.FileInfo.fromString('a.xlsx', 'nested')},
      {folders: ['..'], file: DG.FileInfo.fromString('a.xlsx', 'up')},
    ]));
    expectDeepEqual(Object.keys(zip).sort(),
      ['A (2).xlsx', 'Sub_ 1/Sub 2/a.xlsx', '_/a.xlsx', 'a (3).xlsx', 'a.xlsx', 'b'].sort());
    expect(new TextDecoder().decode(zip['Sub_ 1/Sub 2/a.xlsx']), 'nested');
  });

  test('The dialog lists the models with their default export, with and without subspaces', async () => {
    const run = (nqName: string, spaceId: string, outdated = false): ExportRun =>
      ({run: callOf(nqName, {}, outdated), model: DG.Func.byName(nqName), spaceId, folders: []});
    const runs = [
      run(EXPORT_MODEL, 'top'), run(EXPORT_MODEL, 'top', true), run(EXPORT_MODEL, 'sub'),
      run('Compute2:TestAdd2', 'sub'), run('Compute2:TestMul2', 'sub', true),
    ];
    const {root, unmount} = mount(BulkExportDialog, {locationName: 'Top', runs, spaceId: 'top', filtered: true,
      unlisted: [{spaceId: 'broken', name: 'Broken'}]});
    const warnings = () => Array.from(root.querySelectorAll('.c2-export-warning')).map((el) => el.textContent);
    try {
      const cells = () => Array.from(root.querySelectorAll('.c2-export-models > span')).map((el) => el.textContent);
      await awaitCheck(() => !cells().includes('Loading…'), 'the exports should load', 5000);
      expect(root.querySelector<HTMLSelectElement>('.c2-export-models select')!.value, DEFAULT_EXPORT);
      // a row is the model, its runs and its export
      const rows = cells().reduce<(string | null)[][]>((acc, cell, idx) =>
        idx % 3 === 0 ? [...acc, [cell]] : [...acc.slice(0, -1), [...acc.at(-1)!, cell]], []);
      expectDeepEqual(rows.map(([, count, choice]) => [count, choice]), [
        ['3 runs, 1 not completed', 'Default Excelrecnonehook'], ['1 run', 'Default Excel'], ['1 run', 'Not completed'],
      ]);
      expect(root.querySelector('.c2-runs-footer-note')!.textContent, '3 runs to export');
      expectDeepEqual(warnings(), ['Couldn\'t list Broken, so its runs aren\'t exported.',
        'Search or filters are on: only the runs they show are exported.']);

      const toggle = root.querySelector<HTMLInputElement>('.c2-export-toggle input')!;
      toggle.checked = false;
      toggle.dispatchEvent(new Event('change'));
      await Vue.nextTick();
      expect(root.querySelector('.c2-runs-footer-note')!.textContent, '1 run to export');
      expect(warnings().length, 1, 'a subspace that couldn\'t be listed doesn\'t matter without subspaces');
    } finally {
      unmount();
    }
  });

  test('Export is grey without runs, and takes what the search keeps', () => withRunStorage(async (f) => {
    const name = uniqueName('Bulk browser');
    await f.save(callOf(EXPORT_MODEL, {a: 1}, false), {metadata: {name}});
    await f.save(callOf(EXPORT_MODEL, {a: 2}, true), {metadata: {name: `${name} draft`}});
    const {root, unmount} = mount(RunBrowser, {model: DG.Func.byName(EXPORT_MODEL)});
    const exportButton = () => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-buttons button'))
      .find((el) => el.textContent === 'Export')!;
    const search = (value: string) => {
      const input = root.querySelector<HTMLInputElement>('.c2-runs-search')!;
      input.value = value;
      input.dispatchEvent(new Event('input'));
    };
    const row = (text: string) => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-list .c2-runs-row'))
      .find((el) => el.querySelector('.c2-runs-row-name')?.textContent === text);
    try {
      await awaitCheck(() => !exportButton().classList.contains('d4-disabled'), 'Export should take My runs', 15000);
      await awaitCheck(() => !!row(`${name} draft`), 'the draft should be listed', 15000);
      expect(!!row(`${name} draft`)!.querySelector('.c2-runs-not-completed'), true, 'the draft is not completed');
      expect(!!row(name)!.querySelector('.c2-runs-not-completed'), false);
      search(uniqueName('Nothing like it'));
      await awaitCheck(() => exportButton().classList.contains('d4-disabled'), 'Export should be grey', 5000);
      expect(exportButton().title, 'No runs to export');

      search(name);
      await awaitCheck(() => !exportButton().classList.contains('d4-disabled'), 'Export should take the run', 5000);
      exportButton().click();
      const dialog = () => document.querySelector<HTMLElement>('.c2-export-dialog');
      await awaitCheck(() => !!dialog()?.querySelector('.c2-export-warning'), 'the dialog should warn', 5000);
      expect(dialog()!.querySelector('.c2-runs-footer-note')!.textContent, '1 run to export');
      expect(dialog()!.querySelector('.c2-export-models')!.textContent!.includes('2 runs, 1 not completed'), true);
      expect(!!dialog()!.querySelector('.c2-export-toggle'), false);
      Array.from(dialog()!.querySelectorAll<HTMLElement>('button')).find((el) => el.textContent === 'Cancel')!.click();
      await awaitCheck(() => !dialog(), 'the dialog should close', 5000);
    } finally {
      unmount();
    }
  }));

  test('Saved runs export as their menus do', () => withRunStorage(async (f) => {
    const call = DG.Func.byName(EXPORT_MODEL).prepare({a: 5});
    await call.call();
    const saved = await f.save(call, {metadata: {name: uniqueName('Bulk call')}});
    expect(text(await exportSavedCall(saved.id, 'rec')), `${EXPORT_MODEL}|5`);
    const excel = await exportSavedCall(saved.id, DEFAULT_EXPORT);
    expect(excel.fileName.startsWith(`${EXPORT_MODEL} - `), true, `unexpected file ${excel.fileName}`);
    expect(excel.fileName.endsWith('.xlsx') && excel.data.length > 0, true);

    const workflow = await saveWorkflow('Compute2:MockTreeWizardPipeline');
    try {
      expect(text((await exportSavedWorkflow(workflow.id, 'recorder'))!).startsWith('Exported treeWizardRoot'), true);
      const zip = (await exportSavedWorkflow(workflow.id, DEFAULT_EXPORT))!;
      expect(zip.fileName.endsWith('.zip'), true);
      expect(Object.keys(unzipSync(zip.data)).includes('000_summary.xlsx'), true);
    } finally {
      await deleteSavedWorkflow(workflow.id);
    }

    const compact = await saveWorkflow('Compute2:MockSingleStepPipeline');
    try {
      const step = (await exportSavedWorkflow(compact.id, DEFAULT_EXPORT))!;
      expect(step.fileName.endsWith('.xlsx'), true, `a compact view exports its step, got ${step.fileName}`);
    } finally {
      await deleteSavedWorkflow(compact.id);
    }
  }));

  test('Runs of a custom view export with the exports of its function', () => withRunStorage(async (f) => {
    const customView = DG.Func.byName(CUSTOM_VIEW);
    expectDeepEqual((await modelExports(customView, modelKind(customView, false))).map(({id}) => id),
      [DEFAULT_EXPORT, 'rec']);
    const call = DG.Func.byName('Compute2:TestAdd2').prepare({a: 4, b: 1});
    await call.call();
    call.options[CUSTOM_VIEW_OPTION] = CUSTOM_VIEW;
    const saved = await f.save(call, {metadata: {name: uniqueName('Bulk custom view')}});
    const {files, failed} = await exportRuns([{item: {run: saved, model: customView, folders: []}, kind: 'call',
      exportId: 'rec'}], {onProgress: () => {}, isCancelled: () => false});
    expectDeepEqual(failed, []);
    expect(text(files[0].file), 'Compute2:TestAdd2|4');

    const hooksOf: string[] = [];
    await exportSavedCall(saved.id, DEFAULT_EXPORT, customView, async (func) => {
      hooksOf.push(func.name);
      return undefined;
    });
    expectDeepEqual(hooksOf, ['TestAdd2'], {prefix: 'the report takes the hook of the run\'s function'});
  }));

  test('A step\'s history exports with the step\'s hook in place of its function\'s', () => withRunStorage(async (f) => {
    const model = DG.Func.byName('Compute2:ViewersHookModel');
    const saved = await f.save(await model.prepare({n: 2}).call(), {metadata: {name: uniqueName('Bulk step hook')}});
    const applied: string[] = [];
    const {files, failed} = await exportRuns([{item: {run: saved, model, folders: []}, kind: 'call',
      exportId: DEFAULT_EXPORT}], {onProgress: () => {}, isCancelled: () => false,
      viewersHook: (ioName: string) => applied.push(ioName)});
    expectDeepEqual(failed, []);
    expect(files.length, 1);
    expectDeepEqual(applied, ['df']);
  }));

  test('A workflow run saves whether it is completed, and one with outdated steps is not exported', async () => {
    const run = await saveWorkflow('Compute2:MockTreeWizardPipeline');
    try {
      expect(readRunCompleted(run), true);
    } finally {
      await deleteSavedWorkflow(run.id);
    }

    const config = (await loadWorkflowConfig('Compute2:MockPipeline2').toPromise())!;
    const tree = StateTree.fromPipelineConfig({config});
    await tree.init().toPromise();
    const unrun: DG.FuncCall = (await tree.save().toPromise())!;
    tree.close();
    try {
      expect(readRunCompleted(unrun), false);
      expect(await exportSavedWorkflow(unrun.id, DEFAULT_EXPORT), null);
    } finally {
      await deleteSavedWorkflow(unrun.id);
    }
  });
});
