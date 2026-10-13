import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import type {
  PipelineSerializedState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {
  parseComparisonDefaults, entryFromDataFrame, collectWorkflowSteps,
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
