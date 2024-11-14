import * as DG from 'datagrok-api/dg';
import {category, test, before, delay} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {callHandler, makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';

const config1: PipelineConfiguration = {
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
    handler({controller}) {
      controller.setAll('out1', 2, 'restricted');
      return;
    },
  }],
};


category('ComputeUtils: Driver instance additional states', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Propagate validations info to view state', async () => {
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
          controller.setValidation('out1', makeValidationResult({warnings: ['test warn']}));
          return;
        },
      }, {
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step1/b',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['another test warn']}));
          return;
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const ls = new LinksState();
      const [link1, link2] = ls.createStateLinks(tree.nodeTree);
      link1.wire(tree.nodeTree);
      link2.wire(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const validations = tree.getValidations();
      cold('-a').subscribe(() => {
        link1.trigger();
        link2.trigger();
      });
      const a = {};
      const b = {
        'a': {
          'errors': [],
          'warnings': [
            {
              'description': 'test warn',
            },
          ],
          'notifications': [],
        },
      };
      const c = {
        'a': {
          'errors': [],
          'warnings': [
            {
              'description': 'test warn',
            },
          ],
          'notifications': [],
        },
        'b': {
          'errors': [],
          'warnings': [
            {
              'description': 'another test warn',
            },
          ],
          'notifications': [],
        },
      };
      expectObservable(validations[inNode.getItem().uuid], '^ 1000ms !').toBe('a(bc)', {a, b, c});
    });
  });


  test('Propagate consistency info to view state', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const ls = new LinksState();
      const [link1] = ls.createStateLinks(tree.nodeTree);
      link1.wire(tree.nodeTree);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      const consistency = tree.getConsistency();
      cold('-a').subscribe(() => {
        link1.trigger();
      });
      cold('--a').subscribe(() => {
        outNode.getItem().getStateStore().editState('a', 3);
      });
      const a = {};
      const b = {};
      const c = {
        'a': {
          'restriction': 'restricted',
          'inconsistent': true,
          'assignedValue': 2,
        },
      };
      expectObservable(consistency[outNode.getItem().uuid], '^ 1000ms !').toBe('abc', {a, b, c});
    });
  });

  test('Propagate funccalls state to view state', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      const states = tree.getFuncCallStates();
      cold('-a').subscribe(() => {
        const fcnode = node.getItem() as FuncCallNode;
        fcnode.getStateStore().setState('a', 1);
        fcnode.getStateStore().setState('b', 2);
        fcnode.getStateStore().run({'res': 3}, 5).subscribe();
      });
      const a = {
        'isRunning': false,
        'isRunnable': true,
        'isOutputOutdated': true,
        'pendingDependencies': [],
        'runError': undefined,
      };
      const b = {
        'isRunning': true,
        'isRunnable': true,
        'isOutputOutdated': true,
        'pendingDependencies': [],
        'runError': undefined,
      };
      const c = {
        'isRunning': true,
        'isRunnable': false,
        'isOutputOutdated': true,
        'pendingDependencies': [],
        'runError': undefined,
      };
      const d = {
        'isRunning': true,
        'isRunnable': false,
        'isOutputOutdated': false,
        'pendingDependencies': [],
        'runError': undefined,
      };
      const e = {
        'isRunning': false,
        'isRunnable': false,
        'isOutputOutdated': false,
        'pendingDependencies': [],
        'runError': undefined,
      };
      const f = {
        'isRunning': false,
        'isRunnable': true,
        'isOutputOutdated': false,
        'pendingDependencies': [],
        'runError': undefined,
      };
      expectObservable(states[node.getItem().uuid], '^ 1000ms !').toBe('a(bc)-(def)', {a, b, c, d, e, f});
    });
  });

  test('Propage nodes descriptions', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
      ],
      links: [{
        id: 'selector',
        type: 'selector',
        from: 'in:step1/a',
        to: ['out1:name', 'out2:description', 'out3:tags'],
        handler({controller}) {
          const val = controller.getFirst('in');
          controller.setDescriptionItem('out1', `Title ${val}`);
          controller.setDescriptionItem('out2', `Description ${val}`);
          controller.setDescriptionItem('out3', [`tag ${val}`]);
        }
      }]
        }
    const pconf = await getProcessedConfig(config5);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      const pipeline = tree.nodeTree.root;
      cold('--a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      cold('----a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 2);
      });
      expectObservable(tree.getNodesDescriptions()[pipeline.getItem().uuid]).toBe('a-b-c', {
        a: {
          "name": undefined,
          "description": undefined,
          "tags": undefined,
        },
        b: {
          "name": "Title 1",
          "description": "Description 1",
          "tags": [
            "tag 1"
          ]
        },
        c: {
          "name": "Title 2",
          "description": "Description 2",
          "tags": [
            "tag 2"
          ]
        }
      });
    });
  });

  test('Propagate consistency info to RO view state', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, isReadonly: true});
      tree.init().subscribe();
      const ls = new LinksState();
      const [link1] = ls.createStateLinks(tree.nodeTree);
      link1.wire(tree.nodeTree);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      const consistency = tree.getConsistency();
      cold('-a').subscribe(() => {
        link1.trigger();
      });
      cold('--a').subscribe(() => {
        outNode.getItem().getStateStore().editState('a', 3);
      });
      const a = {};
      const b = {
        'a': {
          'restriction': 'restricted',
          'inconsistent': true,
          'assignedValue': 2,
        },
      };
      expectObservable(consistency[outNode.getItem().uuid], '^ 1000ms !').toBe('abb', {a, b});
    });
  });

  test('Restore saved consistency state', async () => {
    const conf = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(conf);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const outNode = tree.nodeTree.getNode([{idx: 1}]);
    outNode.getItem().getStateStore().setState('a', 10, 'restricted');
    const metaCallSaved = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCallSaved!.id, config: pconf, isReadonly: false}).toPromise();
    await loadedTree.init().toPromise();
    const outNodeLoaded = loadedTree.nodeTree.getNode([{idx: 1}]);
    outNodeLoaded.getItem().getStateStore().editState('a', 3);
    const consistency = loadedTree.getConsistency();
    expectDeepEqual(consistency[outNodeLoaded.getItem().uuid]?.value, {
      'a': {
        'restriction': 'restricted',
        'inconsistent': true,
        'assignedValue': 10,
      },
    });
  });

  test('Restore saved outdated state', async () => {
    const conf = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(conf);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const node = tree.nodeTree.getNode([{idx: 0}]);
    const fcnode = node.getItem() as FuncCallNode;
    fcnode.getStateStore().setState('a', 1);
    fcnode.getStateStore().setState('b', 2);
    await fcnode.getStateStore().run().toPromise();
    const metaCallSaved = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCallSaved!.id, config: pconf, isReadonly: false}).toPromise();
    await loadedTree.init().toPromise();
    const nodeLoaded = loadedTree.nodeTree.getNode([{idx: 0}]);
    const states = loadedTree.getFuncCallStates();
    expectDeepEqual(states[nodeLoaded.getItem().uuid]?.value, {
      'isRunning': false,
      'isRunnable': true,
      'isOutputOutdated': false,
    });
  });

  test('Handle script errors', async () => {
    const config2 = await callHandler<PipelineConfiguration>('LibTests:MockProvider4', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config2);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const node = tree.nodeTree.getNode([{idx: 0}]);
    const fcnode = node.getItem() as FuncCallNode;
    fcnode.getStateStore().setState('a', -1);
    fcnode.getStateStore().setState('b', -2);
    await fcnode.getStateStore().run().toPromise();
    expectDeepEqual(fcnode.funcCallState$.value, {
      'isRunning': false,
      'isRunnable': true,
      'isOutputOutdated': true,
      'runError': 'Error: Test error',
      'pendingDependencies': [],
    });
  });

  test('Remove script errors on success', async () => {
    const config2 = await callHandler<PipelineConfiguration>('LibTests:MockProvider4', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config2);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const node = tree.nodeTree.getNode([{idx: 0}]);
    const fcnode = node.getItem() as FuncCallNode;
    fcnode.getStateStore().setState('a', -1);
    fcnode.getStateStore().setState('b', -2);
    await fcnode.getStateStore().run().toPromise();
    fcnode.getStateStore().setState('a', 1);
    fcnode.getStateStore().setState('b', 2);
    await fcnode.getStateStore().run().toPromise();
    expectDeepEqual(fcnode.funcCallState$.value, {
      'isRunning': false,
      'isRunnable': true,
      'isOutputOutdated': false,
      'pendingDependencies': [],
    });
  });

  test('Restore saved error state', async () => {
    const config2 = await callHandler<PipelineConfiguration>('LibTests:MockProvider4', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config2);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const node = tree.nodeTree.getNode([{idx: 0}]);
    const fcnode = node.getItem() as FuncCallNode;
    fcnode.getStateStore().setState('a', -1);
    fcnode.getStateStore().setState('b', -2);
    await fcnode.getStateStore().run().toPromise();
    const metaCallSaved = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCallSaved!.id, config: pconf, isReadonly: false}).toPromise();
    await loadedTree.init().toPromise();
    const nodeLoaded = loadedTree.nodeTree.getNode([{idx: 0}]);
    expectDeepEqual((nodeLoaded.getItem() as FuncCallNode).funcCallState$.value, {
      'isRunning': false,
      'isRunnable': true,
      'isOutputOutdated': true,
      'runError': 'Error: Test error',
      'pendingDependencies': [],
    });

  });
});
