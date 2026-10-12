import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {createTestScheduler} from '../../../test-utils';

const config: PipelineConfiguration = {
  id: 'pipeline',
  type: 'static',
  steps: [{id: 'step', nqName: 'LibTests:TestToleranceInputs'}],
};

category('ComputeUtils: Driver consistency tolerance', async () => {
  let testScheduler: TestScheduler;
  let pconf: PipelineConfigurationProcessed;

  before(async () => {
    testScheduler = createTestScheduler();
    pconf = await getProcessedConfig(config);
  });

  // a link assigns `assigned`, then the user edits the input to `edited`
  function isInconsistent(input: string, assigned: any, edited: any) {
    let inconsistent: boolean | undefined;
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step = () => tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      cold('-a').subscribe(() => step().getStateStore().setState(input, assigned));
      cold('--a').subscribe(() => step().getStateStore().editState(input, edited));
      cold('----a').subscribe(() => inconsistent = step().consistencyInfo$.value[input]?.inconsistent);
    });
    return inconsistent;
  }

  const jsonError = (text: string) => {
    try {
      JSON.parse(text);
    } catch (e) {
      return (e as Error).message;
    }
  };

  const df = (value: number) => DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', new Float32Array([value]))]);

  test('Without an annotation small changes stay within the default tolerance', async () => {
    expectDeepEqual(isInconsistent('plain', 1e-6, 5e-5), false);
    expectDeepEqual(isInconsistent('plain', 1e-6, 1e-3), true);
  });

  test('An absolute tolerance annotation replaces the default', async () => {
    expectDeepEqual(isInconsistent('abs', 1e-6, 5e-5), true);
    expectDeepEqual(isInconsistent('abs', 1e-6, 1e-6 + 1e-10), false);
  });

  test('A relative tolerance annotation scales with the value', async () => {
    expectDeepEqual(isInconsistent('plain', 1e6, 1e6 + 0.5), true, {prefix: 'Default flags float noise'});
    expectDeepEqual(isInconsistent('rel', 1e6, 1e6 + 0.5), false);
    expectDeepEqual(isInconsistent('rel', 1e6, 1e6 + 5), true);
  });

  test('A tolerance annotation applies to dataframe cells', async () => {
    expectDeepEqual(isInconsistent('df', df(1e-6), df(5e-5)), true);
    expectDeepEqual(isInconsistent('df', df(1e-6), df(1e-6)), false);
  });

  test('A relative tolerance wins over an absolute one', async () => {
    expectDeepEqual(isInconsistent('both', 1e-3, 0.5), true);
    expectDeepEqual(isInconsistent('both', 1e6, 1e6 + 0.5), false);
  });

  const cols = ({conc = 1, other = 1, big = 1e6, n = 1}) => DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'conc', [conc]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'other', [other]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'big', [big]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', [n]),
  ]);

  test('Column tolerances match names ignoring case', async () => {
    expectDeepEqual(isInconsistent('dfCols', cols({conc: 1e-6}), cols({conc: 5e-5})), true);
    expectDeepEqual(isInconsistent('dfCols', cols({conc: 1e-6}), cols({conc: 1e-6 + 1e-10})), false);
  });

  test('The wildcard covers columns not listed', async () => {
    expectDeepEqual(isInconsistent('dfCols', cols({other: 1}), cols({other: 2})), false);
    expectDeepEqual(isInconsistent('dfCols', cols({other: 1}), cols({other: 4})), true);
  });

  test('A relative column tolerance wins over the absolute wildcard', async () => {
    expectDeepEqual(isInconsistent('dfCols', cols({big: 1e6}), cols({big: 1e6 + 0.5})), false);
    expectDeepEqual(isInconsistent('dfCols', cols({big: 1}), cols({big: 2})), true);
  });

  test('Int columns compare exactly', async () => {
    expectDeepEqual(isInconsistent('dfCols', cols({n: 1}), cols({n: 2})), true);
    expectDeepEqual(isInconsistent('dfCols', cols({n: 1}), cols({n: 1})), false);
  });

  test('Malformed tolerance annotations are warned about and ignored', async () => {
    const logger = new DriverLogger();
    const bad = await getProcessedConfig({
      id: 'pipeline', type: 'static', steps: [{id: 'step', nqName: 'LibTests:TestToleranceBadInputs'}],
    }, logger) as any;
    expectDeepEqual(bad.steps[0].io.filter((io: any) => io.tolerance).map((io: any) => io.id), []);
    expectDeepEqual(logger.errors.map((item) => [item.severity, item.message]), [
      ['warning', 'LibTests:TestToleranceBadInputs, input notNumber: small is not a number'],
      ['warning', 'LibTests:TestToleranceBadInputs, input columns: column tolerances need a dataframe input'],
      ['warning',
        'LibTests:TestToleranceBadInputs, input count: ' +
        'consistency tolerance applies to double and dataframe inputs only'],
      ['warning', `LibTests:TestToleranceBadInputs, input badJson: ${jsonError('{"x": small}')}`],
      ['warning', 'LibTests:TestToleranceBadInputs, input badColumn: column x: 1e-9 is not a number'],
    ]);
  });
});
