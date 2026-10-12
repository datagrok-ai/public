import * as DG from 'datagrok-api/dg';
import type ExcelJS from 'exceljs';
import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import type {
  ExportCbInput, ExportSummaryItem,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import type {
  FuncCallStateInfo,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {friendlyIoName, getExportSummary, reportSummary, reportTree, SUMMARY_FILE_NAME} from '../utils';

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

  test('A step without a call is not loaded', async () => {
    const item = getExportSummary(treeState, states).find((item) => item.name === 'notLoaded')!;
    expect(item.fileName === undefined, true, 'no file name');
    expect(item.status === undefined, true, 'no status');
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
    const [, zip, , summary] = await reportTree({
      startDownload: false, treeState, ...states, cb: async (input) => {cbInputs.push(input);},
    });
    const stepFiles = summary.filter((item) => item.fileName).map((item) => [...item.path, item.fileName].join('/'));
    expectDeepEqual(Object.keys(zip).sort(), [SUMMARY_FILE_NAME, ...stepFiles].sort());
    expectDeepEqual(summary, getExportSummary(treeState, states));
    expectDeepEqual(cbInputs.map((input) => [input.fileName, input.status]),
      summary.filter((item) => item.fileName).map((item) => [item.fileName, item.status]));
  });

  test('Summary workbook has one row per item', async () => {
    const [, zip] = await reportTree({startDownload: false, treeState, ...states});
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
    const [, zip] = await reportTree({startDownload: false, treeState, ...states});
    const standalone = sheetRows(await loadWorkbook(new Uint8Array(await blob.arrayBuffer())));
    const exported = sheetRows(await loadWorkbook((zip[SUMMARY_FILE_NAME] as [Uint8Array, any])[0]));
    expectDeepEqual(standalone, exported);
  });
});
