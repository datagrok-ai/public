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


category('ComputeUtils: Driver links reactivity: validators', async () => {
  let testScheduler: TestScheduler;

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
        controller.setValidation('out1', ({warnings: [{description: 'test warn'}]}));
        return;
      },
    }],
  };

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Run validators with debounce', async () => {
    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        link.setActive();
        inNode.getItem().getStateStore().setState('a', 1);
      });
      expectObservable(link.isRunning$, '^ 1000ms !').toBe('a 250ms (bc)', {a: false, b: true, c: false});
      expectObservable((inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$).toBe('a 250ms b', {
        a: {},
        b: {
          [link.uuid]: {
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

  test('Run validators on trigger', async () => {
    const pconf = await getProcessedConfig(config2);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        link.trigger();
      });
      expectObservable(link.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      expectObservable((inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$).toBe('a b', {
        a: {},
        b: {
          [link.uuid]: {
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

  test('Links propagate multiple validators results', async () => {
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
          controller.setValidation('out1', ({warnings: [{description: 'some warn'}]}));
          return;
        },
      }, {
        id: 'link2',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'another warn'}]}));
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link1, link2] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link1.wire(tree.nodeTree);
      link2.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        link1.trigger();
        link2.trigger();
      });
      expectObservable(link1.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      expectObservable(link2.isRunning$, '^ 1000ms !').toBe('a (bc)', {a: false, b: true, c: false});
      expectObservable((inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$).toBe('a(bc)', {
        a: {},
        b: {
          [link1.uuid]: {
            a: {'warnings': [
              {
                'description': 'some warn',
              }],
            }},
        },
        c: {
          [link1.uuid]: {
            a: {'warnings': [
              {
                'description': 'some warn',
              }],
            }},
          [link2.uuid]: {
            a: {'warnings': [
              {
                'description': 'another warn',
              }],
            },
          },
        },
      });
    });
  });

  test('Run default validators', async () => {
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
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
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
      cold('--a').subscribe(() => {
        n2.getItem().getStateStore().setState('b', undefined);
      });
      expectObservable((n1.getItem() as FuncCallNode).validationInfo$, '-^ 1000ms !').toBe('-(abc)', {
        a: {
          'a': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
          'b': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
        },
        b: {
          'b': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
        },
        c: {},
      });
      expectObservable((n2.getItem() as FuncCallNode).validationInfo$.pipe(take(3)), '-^ 1000ms !').toBe('-(abc|)', {
        a: {
          'a': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
          'b': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
        },
        b: {
          'b': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
        },
        c: {},
      });
      expectObservable((n2.getItem() as FuncCallNode).validationInfo$.pipe(skip(3)), '-^ 1000ms !').toBe('--a', {
        a: {
          'b': {
            'errors': [
              {
                'description': 'Missing value',
              },
            ],
            'warnings': [],
            'notifications': [],
          },
        },
      });
      expectObservable((n1.getItem() as FuncCallNode).funcCallState$, '-^ 1000ms !').toBe('-(ab)', {
        a: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
      expectObservable((n2.getItem() as FuncCallNode).funcCallState$.pipe(take(2)), '-^ 1000ms !').toBe('-(ab|)', {
        a: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
        b: {
          'isRunning': false,
          'isRunnable': true,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
      expectObservable((n2.getItem() as FuncCallNode).funcCallState$.pipe(skip(2)), '-^ 1000ms !').toBe('--a', {
        a: {
          'isRunning': false,
          'isRunnable': false,
          'isOutputOutdated': true,
          'pendingDependencies': [],
        },
      });
    });
  });
});
