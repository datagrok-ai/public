import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {first} from 'rxjs/operators';
import {unzipSync} from 'fflate';
import dayjs from 'dayjs';
import type ExcelJS from 'exceljs';
import {category, test, before, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import type {
  ExportCbInput, ExportSummaryItem, ExportUtils,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import type {
  FuncCallStateInfo,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {getExportSummary, reportSummary, reportTree, SUMMARY_FILE_NAME} from '../components/TreeWizard/tree-export';
import {DEFAULT_EXPORT, TreeExportContext, exportCall, exportWorkflow, treeExportUtils} from '../run-exports';
import {friendlyIoName} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import {replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';

const callState = (isOutputOutdated: boolean, runError?: string): FuncCallStateInfo =>
  ({isRunning: false, isRunnable: true, isOutputOutdated, runError, pendingDependencies: []});

const step = (uuid: string, funcCall?: DG.FuncCall) =>
  ({type: 'funccall', uuid, configId: uuid, friendlyName: uuid, isReadonly: false, funcCall});

const workflow = (uuid: string, steps: any[]) =>
  ({type: 'static', uuid, configId: uuid, friendlyName: uuid, isReadonly: false, steps}) as unknown as PipelineState;

async function loadWorkbook(data: Uint8Array): Promise<ExcelJS.Workbook> {
  await DG.Utils.loadJsCss(['/js/common/exceljs.min.js']);
  //@ts-ignore
  const wb = new window.ExcelJS.Workbook() as ExcelJS.Workbook;
  await wb.xlsx.load(data.slice().buffer);
  return wb;
}

const sheetRows = (wb: ExcelJS.Workbook) => {
  const rows: string[][] = [];
  wb.getWorksheet('Summary')!.eachRow((row) => rows.push((row.values as any[]).slice(1).map((v) => String(v ?? ''))));
  return rows;
};

category('Export: workflow summary', () => {
  let fc: DG.FuncCall;
  let treeState: PipelineState;
  let states: any;
  let io: (name: string) => string;

  before(async () => {
    const call = async () => await DG.Func.byName('Compute2:TestAdd2').prepare({a: 1, b: 2}).call();
    fc = await call();
    io = (name) => friendlyIoName(fc, name);
    treeState = workflow('root', [
      step('ok', await call()),
      step('outdated', await call()),
      step('failed', await call()),
      workflow('nested', [
        step('errors', await call()),
        step('inconsistent', await call()),
        step('notLoaded'),
      ]),
    ]);
    states = {
      callInfoStates: {
        ok: callState(false),
        outdated: callState(true),
        failed: callState(true, 'Error: boom'),
        errors: callState(true),
        inconsistent: callState(false),
      },
      validationStates: {
        outdated: {a: {warnings: ['Too high'], notifications: [{description: 'Rounded'}]}},
        errors: {b: {errors: ['Required']}},
      },
      consistencyStates: {
        inconsistent: {a: {restriction: 'restricted', inconsistent: true, assignedValue: 5}},
      },
      pipelineValidations: {nested: {errors: ['Need two items']}},
      descriptions: {ok: {title: 'first'}},
    };
  });

  test('Summary lists workflows and steps in tree order', async () => {
    const items = getExportSummary(treeState, states);
    expectDeepEqual(items.map((item) => [item.kind, item.name, item.path.join('/')]), [
      ['workflow', 'root', ''],
      ['step', 'ok', ''],
      ['step', 'outdated', ''],
      ['step', 'failed', ''],
      ['workflow', 'nested', '004_nested'],
      ['step', 'errors', '004_nested'],
      ['step', 'inconsistent', '004_nested'],
      ['step', 'notLoaded', '004_nested'],
    ]);
  });

  test('Step items carry the tree status, run error and messages', async () => {
    const byName = Object.fromEntries(getExportSummary(treeState, states).map((item) => [item.name, item]));
    expectDeepEqual(byName.ok.status, 'succeeded');
    expectDeepEqual(byName.ok.title, 'first');
    expectDeepEqual(byName.outdated.status, 'next warn');
    expectDeepEqual(byName.outdated.warnings, [`${io('a')}: Too high`]);
    expectDeepEqual(byName.outdated.notifications, [`${io('a')}: Rounded`]);
    expectDeepEqual(byName.failed.status, 'failed');
    expectDeepEqual(byName.failed.runError, 'Error: boom');
    expectDeepEqual(byName.errors.status, 'next error');
    expectDeepEqual(byName.errors.errors, [`${io('b')}: Required`]);
    expectDeepEqual(byName.inconsistent.status, 'succeeded inconsistent');
    expectDeepEqual(byName.inconsistent.inconsistentInputs, [io('a')]);
  });

  test('A title given as a list is left out, as in the tree', async () => {
    const listTitle = {...states, descriptions: {...states.descriptions, ok: {title: ['a', 'b']}}};
    const item = getExportSummary(treeState, listTitle).find((item) => item.name === 'ok')!;
    expect(item.title === undefined, true, 'no title');
    expect(item.fileName!.includes('a,b'), false, `file name ${item.fileName}`);
  });

  test('A step without a call is not loaded', async () => {
    const item = getExportSummary(treeState, states).find((item) => item.name === 'notLoaded')!;
    expect(item.fileName === undefined, true, 'no file name');
    expect(item.status === undefined, true, 'no status');
  });

  test('A step with a call state but no call is not loaded', async () => {
    const withCallState = {...states, callInfoStates: {...states.callInfoStates, notLoaded: callState(true)}};
    const items = getExportSummary(treeState, withCallState);
    expect(items.find((item) => item.name === 'notLoaded')!.status === undefined, true, 'no status');
    expectDeepEqual(items.find((item) => item.name === 'nested')!.rollup?.notLoaded, 1);
  });

  test('Workflow items carry their validations and a rollup of their steps', async () => {
    const byName = Object.fromEntries(getExportSummary(treeState, states).map((item) => [item.name, item]));
    expectDeepEqual(byName.nested.errors, ['Need two items']);
    expectDeepEqual(byName.nested.rollup,
      {steps: 3, notLoaded: 1, failed: 0, outdated: 1, withErrors: 1, withWarnings: 0, inconsistent: 1});
    expectDeepEqual(byName.root.rollup,
      {steps: 6, notLoaded: 1, failed: 1, outdated: 2, withErrors: 1, withWarnings: 1, inconsistent: 1});
  });

  test('Tree export adds the summary file next to the step files', async () => {
    const cbInputs: ExportCbInput[] = [];
    const [, zip, summary] = await reportTree({
      treeState, ...states, cb: async (input) => {cbInputs.push(input);},
    });
    const stepFiles = summary.filter((item) => item.fileName).map((item) => [...item.path, item.fileName].join('/'));
    expectDeepEqual(Object.keys(zip).sort(), [SUMMARY_FILE_NAME, ...stepFiles].sort());
    expectDeepEqual(summary, getExportSummary(treeState, states));
    expectDeepEqual(cbInputs.map((input) => [input.fileName, input.status]),
      summary.filter((item) => item.fileName).map((item) => [item.fileName, item.status]));
  });

  test('Summary workbook has one row per item', async () => {
    const [, zip] = await reportTree({treeState, ...states});
    const rows = sheetRows(await loadWorkbook((zip[SUMMARY_FILE_NAME] as [Uint8Array, any])[0]));
    expectDeepEqual(rows.length, 9);
    expectDeepEqual(rows[0].slice(0, 6), ['Path', 'Type', 'Name', 'File', 'Status', 'Run error']);
    const byName = Object.fromEntries(rows.slice(1).map((row) => [row[2], row]));
    expectDeepEqual(byName['ok - first'][4], 'This step is succeeded');
    expectDeepEqual(byName.failed[5], 'Error: boom');
    expectDeepEqual(byName.notLoaded[4], 'Not loaded');
    expectDeepEqual(byName.nested.slice(0, 2), ['004_nested', 'Workflow']);
    expectDeepEqual(byName.nested[6], 'Need two items');
    expectDeepEqual(byName.nested.slice(10), ['3', '1', '0', '1', '1', '0', '1']);
  });

  test('Standalone summary workbook matches the exported one', async () => {
    const items: ExportSummaryItem[] = getExportSummary(treeState, states);
    const [blob] = await reportSummary(items);
    const [, zip] = await reportTree({treeState, ...states});
    const standalone = sheetRows(await loadWorkbook(new Uint8Array(await blob.arrayBuffer())));
    const exported = sheetRows(await loadWorkbook((zip[SUMMARY_FILE_NAME] as [Uint8Array, any])[0]));
    expectDeepEqual(standalone, exported);
  });
});

category('Export: custom export utils', () => {
  let fc: DG.FuncCall;
  const zipped = (file: DG.FileInfo) => unzipSync(file.data);
  const callName = (call: DG.FuncCall) => `${replaceForWindowsPath(`${call.func.nqName} - ${call.started}`)}.xlsx`;
  // a workflow whose custom export hands back what its handler builds with the utils
  const withExport = (steps: any[],
    handler: (state: PipelineState, utils: ExportUtils, startDownload: boolean) => Promise<DG.FileInfo | Blob>) =>
    ({...(workflow('root', steps) as any), customExports: [{id: 'custom', handler}]}) as PipelineState;

  before(async () => {
    fc = await DG.Func.byName('Compute2:ViewersHookModel').prepare({n: 2}).call();
  });

  test('reportFuncCallExcel applies the step\'s viewers hook and meta, as the step file does', async () => {
    const applied: any[] = [];
    const state = withExport([{...step('hooked', fc), viewersHook: (ioName: string, type: string, _v: any, meta: any) =>
      applied.push([ioName, type, meta])}], async (_state, utils) => {
      const [file] = await utils.reportFuncCallExcel(fc, 'hooked');
      return file;
    });
    await exportWorkflow(state, 'custom', {
      states: {callInfoStates: {hooked: callState(false)}},
      metaStates: {hooked: {df: new BehaviorSubject<any>({color: 'red'})}},
    });
    expectDeepEqual(applied, [['df', 'Grid', {color: 'red'}]]);
  });

  test('reportStateExcel writes a callback\'s workbook changes, unless it replaced the file', async () => {
    const state = withExport([step('first', fc), step('second', fc)], async (state, utils) => {
      const [file] = await utils.reportStateExcel(state, async (input: ExportCbInput) => {
        if (input.fileName.startsWith('001_'))
          input.wb.addWorksheet('Extra');
        else
          input.archive[[...input.path, input.fileName].join('/')] = [new TextEncoder().encode('replaced'), {}];
      });
      return file;
    });
    const zip = zipped(await exportWorkflow(state, 'custom', {
      states: {callInfoStates: {first: callState(false), second: callState(false)}},
    }));
    const file = (prefix: string) => zip[Object.keys(zip).find((key) => key.startsWith(prefix))!];
    expect(!!(await loadWorkbook(file('001_'))).getWorksheet('Extra'), true, 'the added sheet should be in the file');
    expect(new TextDecoder().decode(file('002_')), 'replaced');
  });

  test('A workflow export is told not to download and fails for a Blob', async () => {
    let asked: boolean | undefined;
    const state = withExport([], async (_state, _utils, startDownload) => {
      asked = startDownload;
      return new Blob(['zip']);
    });
    await expectExceptionAsync(async () => {await exportWorkflow(state, 'custom', {states: {}});});
    expect(asked, false);
  });

  test('A call file is named after its run date, or Not completed', async () => {
    expect((await exportCall(fc, DEFAULT_EXPORT, {completed: true})).fileName,
      callName(fc));
    expect((await exportCall(fc, DEFAULT_EXPORT, {completed: false})).fileName,
      'Compute2_ViewersHookModel - Not completed.xlsx');
  });

  test('A workflow file is named after its run date, Not completed, or edited', async () => {
    const started = dayjs();
    const name = async (ctx: Partial<TreeExportContext>) => (await exportWorkflow(workflow('root', [step('s', fc)]),
      DEFAULT_EXPORT, {states: {callInfoStates: {s: callState(false)}}, meta: {title: 't', started}, ...ctx})).fileName;
    expect(await name({completed: true}), `root - t - ${replaceForWindowsPath(`${started}`)}.zip`);
    expect(await name({completed: false}), 'root - t - Not completed.zip');
    expect(await name({completed: true, hasNotSavedEdits: true}), 'root - t - edited.zip');
    expect(await name({completed: false, meta: {title: 't'}}), 'root - t - Not completed.zip');
  });

  test('The default report utils return files named after the run', async () => {
    const names: string[] = [];
    const state = withExport([step('s', fc)], async (state, utils) => {
      const [stepFile] = await utils.reportFuncCallExcel(fc, 's');
      const [zip] = await utils.reportStateExcel(state);
      const [summary] = await utils.reportSummaryExcel(state);
      names.push(stepFile.fileName, zip.fileName, summary.fileName);
      return zip;
    });
    await exportWorkflow(state, 'custom',
      {states: {callInfoStates: {s: callState(false)}}, meta: {title: 't', started: dayjs()}, completed: false});
    expectDeepEqual(names, [callName(fc), 'root - t - Not completed.zip',
      'root - t - Not completed - summary.xlsx']);
  });

  test('A step file is named after its run date, failed, or Not completed', async () => {
    const keys = async (info: FuncCallStateInfo) => Object.keys((await reportTree({
      treeState: workflow('root', [step('s', fc)]), callInfoStates: {s: info}}))[1]);
    expect((await keys(callState(false))).includes(`001_${replaceForWindowsPath(`s - ${fc.started}`)}.xlsx`), true);
    expect((await keys(callState(true))).includes('001_s - Not completed.xlsx'), true);
    expect((await keys(callState(true, 'boom'))).includes('001_s - failed.xlsx'), true);
  });

  test('getFuncCallCustomExports lists the step exports in menu order', async () => {
    const call = DG.Func.byName('Compute2:TestPreferredExportModel').prepare({a: 1});
    expectDeepEqual(treeExportUtils(workflow('root', []), {states: {}}).getFuncCallCustomExports(call),
      ['none', 'hook', 'rec']);
  });

  test('Default Excel takes a table with no columns', async () => {
    const empty = await DG.Func.byName('Compute2:ViewersHookModel').prepare({n: 2}).call();
    const wb = await loadWorkbook((await exportCall(empty, DEFAULT_EXPORT)).data);
    expect(wb.worksheets.length > 0, true);
  });

  test('ReportFuncCallExcel takes a step\'s hook and meta in place of the function\'s hook', async () => {
    const applied: any[] = [];
    await DG.Func.byName('Compute2:ReportFuncCallExcel').apply({funcCall: fc,
      viewersHook: (ioName: string, type: string, _v: any, meta: any) => applied.push([ioName, type, meta]),
      metaState: {df: new BehaviorSubject<any>({color: 'red'})}});
    expectDeepEqual(applied, [['df', 'Grid', {color: 'red'}]]);
  });

  test('A step\'s custom export gets the step\'s hook and meta', async () => {
    const call = DG.Func.byName('Compute2:TestCustomExportModel').prepare({a: 1});
    const hook = () => {};
    const meta = {a: new BehaviorSubject<any>('m')};
    const text = (file: DG.FileInfo) => new TextDecoder().decode(file.data);
    expect(text(await exportCall(call, 'hook', {viewersHook: hook, metaState: meta})), 'function|m');
    const state = withExport([{...step('s', call), viewersHook: hook}],
      async (_state, utils) => utils.runFuncCallCustomExport(call, 's', 'hook'));
    expect(text(await exportWorkflow(state, 'custom', {
      states: {callInfoStates: {s: callState(false)}}, metaStates: {s: meta},
    })), 'function|m');
  });

  test('ReportFuncCallExcel builds a run\'s default Excel with its viewers hook', async () => {
    const applied = grok.events.onCustomEvent('compute2-test-viewers-hook-applied').pipe(first()).toPromise();
    const [file, wb] = await DG.Func.byName('Compute2:ReportFuncCallExcel').apply({funcCall: fc});
    expect(file instanceof DG.FileInfo && file.data.length > 0, true, 'a workbook file');
    expect(file.fileName, callName(fc));
    const [outdated] = await DG.Func.byName('Compute2:ReportFuncCallExcel')
      .apply({funcCall: fc, isOutputOutdated: true});
    expect(outdated.fileName, 'Compute2_ViewersHookModel - Not completed.xlsx');
    expect(wb.worksheets.length > 0, true, 'the workbook itself, to add to');
    expectDeepEqual(await applied, {ioName: 'df', type: 'Grid'});
  });
});
