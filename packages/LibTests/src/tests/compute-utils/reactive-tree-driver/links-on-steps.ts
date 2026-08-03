import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineStepConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver links on steps', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Default data link on a step propagates within the step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfLink',
            from: 'in1:a',
            to: 'out1:b',
          }],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const stepNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        stepNode.getItem().getStateStore().setState('a', 5);
      });
      expectObservable(stepNode.getItem().getStateStore().getStateChanges('b'), '^ 1000ms !').toBe('a b', {a: undefined, b: 5});
    });
  });

  test('Custom handler on a step-level link', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfLink',
            from: 'in1:a',
            to: 'out1:b',
            handler({controller}) {
              const v = controller.getFirst('in1');
              controller.setAll('out1', v * 10);
            },
          }],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const stepNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        stepNode.getItem().getStateStore().setState('a', 3);
      });
      expectObservable(stepNode.getItem().getStateStore().getStateChanges('b'), '^ 1000ms !').toBe('a b', {a: undefined, b: 30});
    });
  });

  test('Validator link on a step fires on referenced input change', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfValidator',
            from: 'in1:a',
            to: 'out1:a',
            type: 'validator',
            handler({controller}) {
              controller.setValidation('out1', {warnings: [{description: 'step warn'}]});
              return;
            },
          }],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const stepNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        link.setActive();
        stepNode.getItem().getStateStore().setState('a', 1);
      });
      expectObservable((stepNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$, '^ 1000ms !').toBe('a 250ms b', {
        a: {},
        b: {
          [link.uuid]: {
            'a': {
              'warnings': [
                {'description': 'step warn'},
              ],
            },
          },
        },
      });
    });
  });

  test('Step-level link does not leak across siblings', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfLink',
            from: 'in1:a',
            to: 'out1:b',
          }],
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step1 = tree.nodeTree.getNode([{idx: 0}]);
      const step2 = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        step2.getItem().getStateStore().setState('a', 7);
      });
      // step1.b must stay untouched: only step1's own 'a' would trigger the link.
      expectObservable(step1.getItem().getStateStore().getStateChanges('b'), '^ 1000ms !').toBe('a', {a: undefined});
      // step2 has no link, so its 'b' must stay undefined too.
      expectObservable(step2.getItem().getStateStore().getStateChanges('b'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('Same step can carry both links and actions', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfLink',
            from: 'in1:a',
            to: 'out1:b',
          }],
          actions: [{
            id: 'action1',
            from: 'in1:a',
            to: 'out1:a',
            position: 'none',
            handler({controller}) {
              controller.setAll('out1', 99);
              return;
            },
          }],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
    await tree.init().toPromise();
    const stepNode = tree.nodeTree.getNode([{idx: 0}]);
    expectDeepEqual([...tree.linksState.links.values()].length, 1);
    expectDeepEqual([...tree.linksState.actions.values()].length, 1);
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    // Action set a to 99; link should have run as a side effect of that change,
    // setting b to 99 (default copy handler).
    expectDeepEqual(stepNode.getItem().getStateStore().getState('a'), 99);
    expectDeepEqual(stepNode.getItem().getStateStore().getState('b'), 99);
  });

  test('processStepConfig emits processed links', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          links: [{
            id: 'selfLink',
            from: 'in1:a',
            to: 'out1:b',
          }],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    const stepConf = (pconf as any).steps[0] as PipelineStepConfigurationProcessed;
    expectDeepEqual(stepConf.links?.length, 1);
    expectDeepEqual(stepConf.links?.[0].id, 'selfLink');
  });
});
