import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {makeValidationResult, PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import { PipelineStateStatic, StepFunCallSerializedState } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

category('ComputeUtils: Driver links data', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      expectDeepEqual(actual, expected);
      // console.log(actual, expected);
    });
  });

  test('Propagate restriction info to serialized state ', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
        handler({controller}) {
          const in1 = controller.getFirst('in1');
          controller.setAll('out1', in1 + 1, 'restricted');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      link.wire(tree);
      cold('-a').subscribe(() => {
        link.enable();
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('--a').subscribe(() => {
        const state = tree.toSerializedState();
        const expected = {
          "a": {
            "type": "restricted",
            "assignedValue": 2
          }
        }
        expectDeepEqual(((state as PipelineStateStatic).steps[1] as StepFunCallSerializedState).inputRestrictions, expected);
      });
    });
  });

  test('Propagate restriction df info to serialized state ', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['a', 'b', 'c'])]);

    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestDF1',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/df',
        handler({controller}) {
          controller.setAll('out1', df, 'restricted');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromConfig({config: pconf, mockMode: true});
      tree.initAll().subscribe();
      const ls = new LinksState();
      const [link] = ls.createLinks(tree);
      const inNode = tree.getNode([{idx: 0}]);
      link.wire(tree);
      cold('-a').subscribe(() => {
        link.enable();
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('--a').subscribe(() => {
        const state = tree.toSerializedState();
        const expected = {
          "df": {
            "type": "restricted",
            "assignedValue": df
          }
        }
        expectDeepEqual(((state as PipelineStateStatic).steps[1] as StepFunCallSerializedState).inputRestrictions, expected);
      });
    });
  });

});
