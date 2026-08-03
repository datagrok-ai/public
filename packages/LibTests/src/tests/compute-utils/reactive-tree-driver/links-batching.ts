import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {skip, take} from 'rxjs/operators';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver links batching', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Batched default validators via runLinks', async () => {
    const pconf = await getProcessedConfig({
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
    });

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true, batchLinks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const n1 = tree.nodeTree.getNode([{idx: 0}]);
      const n2 = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        n1.getItem().getStateStore().setState('a', 1);
        n1.getItem().getStateStore().setState('b', 2);
        n2.getItem().getStateStore().setState('a', 3);
        n2.getItem().getStateStore().setState('b', 4);
      });
      // With batching, all validators should complete — verify final state is correct
      expectObservable((n1.getItem() as FuncCallNode).funcCallState$.pipe(skip(1), take(1)), '-^ 1000ms !').toBe('-(a|)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
      expectObservable((n2.getItem() as FuncCallNode).funcCallState$.pipe(skip(1), take(1)), '-^ 1000ms !').toBe('-(a|)', {
        a: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
    });
  });

  test('Batched validators via reactive path', async () => {
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
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'test warn'}]}));
        },
      }, {
        id: 'link2',
        from: 'in1:step2/a',
        to: 'out1:step2/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'test warn 2'}]}));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const n1 = tree.nodeTree.getNode([{idx: 0}]);
      const n2 = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        n1.getItem().getStateStore().setState('a', 1);
        n2.getItem().getStateStore().setState('a', 2);
      });
      // Both validators should fire via the batch scheduler after a single debounceTime(0)
      expectObservable((n1.getItem().getStateStore() as FuncCallInstancesBridge).validations$.pipe(skip(1), take(1)), '-^ 1000ms !').toBe('-(a|)', {
        a: {
          [tree.linksState.links.values().next().value!.uuid]: {
            'a': {
              'warnings': [
                {
                  'description': 'test warn',
                },
              ],
            },
          },
        },
      });
    });
  });

  test('Sequential opt-out prevents batching', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        sequential: true,
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'sequential warn'}]}));
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      // The sequential validator should NOT be batchable
      expectDeepEqual(links[0].isBatchable, false);
    });
  });

  test('Batched meta links', async () => {
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
        from: 'in1:step1/a',
        to: 'out1:step2/a',
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'val1'});
        },
      }, {
        id: 'link2',
        from: 'in1:step1/b',
        to: 'out1:step2/b',
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'val2'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      tree.init().subscribe();
      const n1 = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      // Verify both meta links are batchable
      const links = [...tree.linksState.links.values()];
      expectDeepEqual(links[0].isBatchable, true);
      expectDeepEqual(links[1].isBatchable, true);
      cold('-a').subscribe(() => {
        n1.getItem().getStateStore().setState('a', 1);
        n1.getItem().getStateStore().setState('b', 2);
      });
      // Both meta links should fire and set their meta values
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).meta.a
      ).toBe('ab', {
        a: undefined,
        b: {
          'key': 'val1',
        },
      });
    });
  });

  test('isBatchable false when batchLinks disabled', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'test'}]}));
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: false});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      // With batchLinks disabled, validators get default debounce (250ms) and are not batchable
      expectDeepEqual(links[0].isBatchable, false);
      expectDeepEqual(tree.linksState.batchLinks, false);
    });
  });

  test('Mixed batch and sequential validators', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'batched warn'}]}));
        },
      }, {
        id: 'link2',
        from: 'in1:step1/b',
        to: 'out1:step1/b',
        type: 'validator',
        sequential: true,
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'sequential warn'}]}));
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      // First link is batchable, second is sequential
      expectDeepEqual(links[0].isBatchable, true);
      expectDeepEqual(links[1].isBatchable, false);
    });
  });

  test('Batched nodemeta links', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'selector1',
        type: 'selector',
        from: 'in:step1/a',
        to: ['out1:title', 'out2:description'],
        handler({controller}) {
          const val = controller.getFirst('in');
          controller.setDescriptionItem('out1', `Title ${val}`);
          controller.setDescriptionItem('out2', `Desc ${val}`);
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, batchLinks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const links = [...tree.linksState.links.values()];
      // selector/nodemeta link should be batchable
      expectDeepEqual(links[0].isBatchable, true);
    });
  });
});
