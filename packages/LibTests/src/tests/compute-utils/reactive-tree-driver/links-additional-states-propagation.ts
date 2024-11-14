import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {snapshotCompare} from '../../../test-utils';
import {LoadedPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {TestScheduler} from 'rxjs/testing';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {switchMap} from 'rxjs/operators';

category('ComputeUtils: Driver links additional states propagation', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Propagate validatations to the state', async () => {
    const config2: PipelineConfiguration = {
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
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['test warn']}));
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      const validations$ = tree.getValidations()[inNode.getItem().uuid];
      expectObservable(validations$).toBe('a 250ms b', {
        a: {},
        b: {
          'a': {
            'warnings': [
              {
                'description': 'test warn',
              },
            ],
            'errors': [],
            'notifications': [],
          },
        },
      });
    });
  });

  test('Propagate meta info to the state', async () => {
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
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'val'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);

      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      const meta$ = tree.getMeta()[outNode.getItem().uuid];
      expectObservable(meta$.pipe(
        switchMap((x) => x.a),
      )).toBe('ab', {
        a: undefined,
        b: {
          'key': 'val',
        },
      });
    });
  });


  test('Propagate consistency info to the state', async () => {
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
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      cold('--a').subscribe(() => {
        outNode.getItem().getStateStore().editState('a', 3);
      });
      const consistency$ = tree.getConsistency()[outNode.getItem().uuid];
      expectObservable(consistency$).toBe('a bc', {
        a: {},
        b: {
          'a': {
            'restriction': 'restricted',
            'inconsistent': false,
            'assignedValue': 2,
          },
        },
        c: {
          'a': {
            'restriction': 'restricted',
            'inconsistent': true,
            'assignedValue': 2,
          },
        },
      });
    });

  });
});
