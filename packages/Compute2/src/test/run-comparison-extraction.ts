import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import type {
  PipelineSerializedState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {CONFIG_PATH} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {
  parseComparisonDefaults, entryFromDataFrame, collectWorkflowSteps, entryFromFuncCall,
} from '../components/RunComparison/entry-extraction';

category('RunComparison: annotation defaults', () => {
  test('parses the comparison JSON annotation', async () => {
    const defaults = parseComparisonDefaults({
      comparison: '{"index": "time", "split": "species", "mode": "timeseries", "units": "s"}',
    });
    expect(defaults.index, 'time');
    expect(defaults.split, 'species');
    expect(defaults.mode, 'timeseries');
    expect(defaults.units, 's');
  });

  test('ignores unknown mode and units values', async () => {
    const defaults = parseComparisonDefaults({comparison: '{"mode": "spiral", "units": "weeks"}'});
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });

  test('malformed JSON produces no defaults', async () => {
    const defaults = parseComparisonDefaults({comparison: '{index: time'});
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
  });

  test('null annotation produces no defaults', async () => {
    // JSON.parse('null') succeeds with null — must not crash on property access
    const defaults = parseComparisonDefaults({comparison: 'null'});
    expect(defaults.index == null, true);
    expect(defaults.mode == null, true);
  });

  test('missing annotations produce no defaults', async () => {
    const defaults = parseComparisonDefaults(undefined);
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });
});

category('RunComparison: raw table entries', () => {
  const makeDf = () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'time', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20, 30]),
    ]);
    df.name = 'data';
    return df;
  };

  test('distinct tables with equal names and row counts get distinct ids', async () => {
    expect(entryFromDataFrame(makeDf()).id === entryFromDataFrame(makeDf()).id, false);
  });

  test('the same table maps to the same id', async () => {
    const df = makeDf();
    expect(entryFromDataFrame(df).id, entryFromDataFrame(df).id);
  });
});

category('RunComparison: workflow steps', () => {
  const step = (configId: string, friendlyName: string, funcCallId: string) =>
    ({type: 'funccall', configId, friendlyName, funcCallId});
  // only the fields the step walk reads
  const state = {
    type: 'static', configId: 'wf', friendlyName: 'Workflow', steps: [
      step('config', 'Configuration', 'c'),
      {type: 'dynamic', configId: 'days', friendlyName: 'Daily cultivation', steps: [
        step('day', 'Day', 'd1'), step('day', 'Day', 'd2'), step('other', 'Other', 'o'),
      ]},
    ],
  } as unknown as PipelineSerializedState;

  test('repeated steps get positional paths and labels', async () => {
    const steps = collectWorkflowSteps(state);
    expectArray(steps.map((s) => s.path), ['wf/config', 'wf/days/day[1]', 'wf/days/day[2]', 'wf/days/other']);
    expectArray(steps.map((s) => s.friendlyPath),
      ['Configuration', 'Daily cultivation · Day 1', 'Daily cultivation · Day 2', 'Daily cultivation · Other']);
    expectArray(steps.map((s) => s.funcCallId), ['c', 'd1', 'd2', 'o']);
  });
});

category('RunComparison: workflow step failures', () => {
  test('steps that cannot be loaded are reported, not thrown', async () => {
    const call = DG.Script.create(
      '//name: rcStepFailureProbe\n//language: javascript\n//output: int x\nx = 1;').prepare();
    const missing = (configId: string, friendlyName: string) =>
      ({type: 'funccall', configId, friendlyName, funcCallId: `${configId}-run-that-does-not-exist`});
    call.options[CONFIG_PATH] = serialize({
      type: 'static', configId: 'wf', friendlyName: 'Workflow',
      steps: [missing('config', 'Configuration'), missing('result', 'Result')],
    });
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'workflow');
    expectArray(entry.failedSteps ?? [], ['Configuration', 'Result']);
    expect(entry.nodes.tables.length + entry.nodes.scalars.length, 0);
  });
});

category('RunComparison: workflow step paths', () => {
  test('nested pipelines join ids and names, unnamed steps fall back to their id', async () => {
    const state = {
      type: 'static', configId: 'wf', friendlyName: 'Workflow', steps: [
        {type: 'static', configId: 'stage', friendlyName: 'Stage', steps: [
          {type: 'static', configId: 'inner', steps: [
            {type: 'funccall', configId: 'calc', friendlyName: 'Calc', funcCallId: 'x'},
          ]},
          {type: 'funccall', configId: 'pending', friendlyName: 'Pending'},
        ]},
      ],
    } as unknown as PipelineSerializedState;
    const steps = collectWorkflowSteps(state);
    expectArray(steps.map((s) => s.path), ['wf/stage/inner/calc']);
    expectArray(steps.map((s) => s.friendlyPath), ['Stage · inner · Calc']);
  });

  test('a root that is a single step uses its own id and name', async () => {
    const solo = {type: 'funccall', configId: 'solo', funcCallId: 's'} as unknown as PipelineSerializedState;
    const [step] = collectWorkflowSteps(solo);
    expect(step.path, 'solo');
    expect(step.friendlyPath, 'solo');
    expect(step.funcCallId, 's');
  });
});

category('RunComparison: function run extraction', () => {
  const probeCall = () => DG.Script.create([
    '//name: rcExtractionProbe',
    '//language: javascript',
    '//input: double a {caption: Data}',
    '//input: double x',
    '//input: string label',
    '//output: double b {caption: Data}',
    '//output: dataframe result',
    '//output: dataframe empty',
  ].join('\n')).prepare({a: 1, x: 2, label: 'l'});

  test('scalars and tables come from inputs and outputs, shared captions get the property name', async () => {
    const call = probeCall();
    call.inputs['x'] = NaN;
    call.outputs['b'] = 3;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1])]);
    call.outputs['result'] = df;
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'function');
    const scalars = new Map(entry.nodes.scalars.map((s) => [s.path, s]));
    expectArray([...scalars.keys()], ['a', 'x', 'b']);
    expect(scalars.get('a')!.name, 'Data (a)');
    expect(scalars.get('b')!.name, 'Data (b)');
    expect(scalars.get('a')!.value, 1);
    expect(scalars.get('x')!.value === null, true);
    expect(scalars.get('b')!.value, 3);
    expectArray(entry.nodes.tables.map((t) => t.path), ['result']);
    expect(entry.dataFrames.get('result')!.dart === df.dart, true);
  });

  test('a malformed workflow config degrades to an entry without steps', async () => {
    const call = probeCall();
    call.options[CONFIG_PATH] = '{not json';
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'workflow');
    expect(entry.nodes.scalars.length + entry.nodes.tables.length, 0);
    expect(entry.failedSteps === undefined, true);
  });
});
