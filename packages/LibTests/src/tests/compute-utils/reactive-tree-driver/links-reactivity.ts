import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of, Subject} from 'rxjs';
import {delay, filter, mapTo, skip, switchMap, take} from 'rxjs/operators';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {makeValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {snapshotCompare} from '../../../test-utils';


category('ComputeUtils: Driver links reactivity', async () => {
  let testScheduler: TestScheduler;

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
      from: 'in1:step1/b',
      to: 'out1:step2/a',
    }],
  };

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

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Run default handler on trigger', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
        link.trigger();
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
    });
  });

  test('Run enabled default handler', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
    });
  });

  test('Dont run disabled handlers', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('Multiple changes buffering', async () => {
    const pconf = await getProcessedConfig(config1);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 1);
        inNode.getItem().getStateStore().setState('b', 2);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
    });
  });

  test('Run custom handlers', async () => {
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
          controller.setAll('out1', in1 + 1);
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
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
    });
  });

  test('Run custom async handlers', async () => {
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
          controller.setAll('out1', in1 + 1);
          return of(null).pipe(delay(100), mapTo(void(0)));
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
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 100ms b', {a: undefined, b: 2});
    });
  });

  test('Properly cancel async handlers', async () => {
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
          controller.setAll('out1', in1 + 1);
          return of(null).pipe(delay(100), mapTo(void(0)));
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
      cold('50ms a').subscribe(() => {
        inNode.getItem().getStateStore().setState('b', 2);
      });
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a 149ms b', {a: undefined, b: 3});
    });
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

  test('Propagate restriction info in a default handler', async () => {
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
        defaultRestrictions: {
          out1: 'restricted',
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
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 1});
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$).toBe('a b', {
        a: {},
        b: {
          'a': {
            'type': 'restricted',
            'assignedValue': 1,
          },
        },
      });
    });
  });

  test('Propagate restriction info in a custom handler', async () => {
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
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a b', {a: undefined, b: 2});
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$).toBe('a b', {
        a: {},
        b: {
          'a': {
            'type': 'restricted',
            'assignedValue': 2,
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
          controller.setValidation('out1', makeValidationResult({warnings: ['some warn']}));
          return;
        },
      }, {
        id: 'link2',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', makeValidationResult({warnings: ['another warn']}));
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


  test('Propagate meta info', async () => {
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
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).meta$.pipe(
        switchMap((x) => x.a),
      )).toBe('ab', {
        a: undefined,
        b: {
          'key': 'val',
        },
      });
    });
  });

  test('Get and run pipeline validation actions', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
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
        actions: 'actions',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
      actions: [{
        id: 'action1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        position: 'none',
        handler({controller}) {
          controller.setAll('out1', 10);
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid);
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 10});
    });
  });

  test('Get and run pipeline validation actions with params', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
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
        actions: 'actions',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
      actions: [{
        id: 'action1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        position: 'none',
        handler({controller}) {
          const val = controller.getAdditionalParam('val');
          controller.setAll('out1', val);
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid, { val: 100 });
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 100});
    });
  });

  test('Get and run funcall validation actions', async () => {
    const s = new Subject<string>();
    const config3: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          actions: [{
            id: 'action1',
            from: 'in1:a',
            to: 'out1:a',
            position: 'none',
            handler({controller}) {
              controller.setAll('out1', 10);
              return;
            },
          }],
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
        actions: 'actions:step1',
        type: 'validator',
        handler({controller}) {
          const v = controller.getFirst('in1');
          if (v === 1) {
            const action = controller.getValidationAction('actions', 'action1');
            s.next(action);
          }
          return;
        },
      }],
    };

    const pconf = await getProcessedConfig(config3);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().setState('a', 1);
      });
      s.pipe(
        delay(1),
      ).subscribe((uuid) => {
        tree.runAction(uuid);
      });
      expectObservable(inNode.getItem().getStateStore().getStateChanges('a')).toBe('ab 250ms c', {a: undefined, b: 1, c: 10});
    });
  });

  test('Run pipeline mutation actions', async () => {
    const config4: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline',
          type: 'parallel',
          stepTypes: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestMul2',
            },
          ],
          initialSteps: [
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'stepr',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1:nestedPipeline',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', {
            id: 'nestedPipeline',
            steps: [
              {
                id: 'step2',
                initialValues: {
                  a: 5,
                },
              },
              {
                id: 'step1',
                initialValues: {
                  a: 10,
                },
              },
            ],
          });
        },
      }],
      links: [{
        id: 'link1',
        from: 'in1:nestedPipeline/all(step1|step2)/a',
        to: 'out1:stepr/a',
        handler({controller}) {
          const v = controller.getAll<number>('in1')!;
          const r = v.reduce((acc, val) => acc + val, 0);
          controller.setAll('out1', r);
        },
      }],
    };

    const pconf = await getProcessedConfig(config4);
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const out = tree.nodeTree.getNode([{idx: 1}]);
      const action = [...tree.linksState.actions.values()][0];
      cold('-a').subscribe(() => {
        tree.runAction(action.uuid).subscribe();
      });
      expectObservable(out.getItem().getStateStore().getStateChanges('a')).toBe('ab', {a: undefined, b: 15});
    });
  });

  test('Run pipeline mutation actions on root pipeline', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'parallel',
      stepTypes: [
        {
          id: 'nestedPipeline1',
          type: 'parallel',
          stepTypes: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestMul2',
            },
          ],
          initialSteps: [
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'step',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', {
            id: 'pipeline1',
            steps: [{
              id: 'nestedPipeline1',
              steps: [
                {
                  id: 'step2',
                  initialValues: {
                    a: 5,
                  },
                },
                {
                  id: 'step1',
                  initialValues: {
                    a: 10,
                  },
                },
              ],
            }]
          });
        },
      }],
    };

    const pconf = await getProcessedConfig(config5);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    await snapshotCompare(state, 'Run pipeline mutation actions on root pipeline');
  });

  test('Run pipeline mutation actions on pipeline', async () => {
    const config5: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'nestedPipeline1',
          type: 'parallel',
          stepTypes: [
            {
              id: 'step1',
              nqName: 'LibTests:TestAdd2',
            },
            {
              id: 'step2',
              nqName: 'LibTests:TestMul2',
            },
          ],
          initialSteps: [
            {
              id: 'step1',
            },
          ],
        },
        {
          id: 'step',
          nqName: 'LibTests:TestSub2',
        },
      ],
      actions: [{
        id: 'action1',
        from: [],
        position: 'none',
        to: 'out1:nestedPipeline1',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1',
            {
              id: 'nestedPipeline1',
              steps: [
                {
                  id: 'step2',
                  initialValues: {
                    a: 5,
                  },
                },
                {
                  id: 'step1',
                  initialValues: {
                    a: 10,
                  },
                },
              ],
            }
          );
        },
      }],
    };

    const pconf = await getProcessedConfig(config5);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    await snapshotCompare(state, 'Run pipeline mutation actions on pipeline');
  });

  test('Select state descriptions', async () => {
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
        to: ['out1:title', 'out2:description', 'out3:tags'],
        handler({controller}) {
          const val = controller.getFirst('in');
          controller.setDescriptionItem('out1', `Title ${val}`);
          controller.setDescriptionItem('out2', `Description ${val}`);
          controller.setDescriptionItem('out3', [`tag ${val}`]);
        },
      }],
    };
    const pconf = await getProcessedConfig(config5);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]);
      const pipeline = tree.nodeTree.root;
      cold('-a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 1);
      });
      cold('--a').subscribe(() => {
        node.getItem().getStateStore().setState('a', 2);
      });
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('title')).toBe('abc',
        {a: undefined, b: 'Title 1', c: 'Title 2'});
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('description')).toBe('abc',
        {a: undefined, b: 'Description 1', c: 'Description 2'});
      expectObservable(pipeline.getItem().nodeDescription.getStateChanges('tags')).toBe('abc',
        {a: undefined, b: ['tag 1'], c: ['tag 2']});
    });
  });

  test('FuncCall pipeline actions', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          initialValues: {
            a: 1,
            b: 2,
          },
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          initialValues: {
            a: 3,
            b: 4,
          },
        },
      ],
      actions: [{
        id: 'action1',
        type: 'funccall',
        from: ['a:step1/a', 'b:step1/b', 'fc(call):step2'],
        to: 'out(call):step2',
        position: 'none',
        async handler({controller}) {
          const a = controller.getFirst('a');
          const b = controller.getFirst('b');
          const fc = controller.getFirst('fc');
          expectDeepEqual(a, 1);
          expectDeepEqual(b, 2);
          expectDeepEqual(fc instanceof DG.FuncCall, true);
          const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
          const nfc = func.prepare({a: 33, b: 44});
          controller.setFuncCall('out', nfc);
        },
      }],
    };

    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    await tree.treeMutationsLocked$.pipe(filter((x) => !x), take(1)).toPromise();
    const node = tree.nodeTree.getNode([{idx: 1}]);
    const a = node.getItem().getStateStore().getState('a');
    const b = node.getItem().getStateStore().getState('b');
    expectDeepEqual(a, 33);
    expectDeepEqual(b, 44);
  });

  test('FuncCall step actions', async () => {
    const config1: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          initialValues: {
            a: 1,
            b: 2,
          },
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
          initialValues: {
            a: 3,
            b: 4,
          },
          actions: [{
            id: 'action1',
            type: 'funccall',
            from: 'fc(call)',
            to: 'out(call)',
            position: 'none',
            async handler({controller}) {
              const fc = controller.getFirst('fc');
              expectDeepEqual(fc instanceof DG.FuncCall, true);
              const func: DG.Func = await grok.functions.eval('LibTests:simpleInputs');
              const nfc = func.prepare({a: 33, b: 44});
              controller.setFuncCall('out', nfc);
            },
          }],
        },
      ],
    };

    const pconf = await getProcessedConfig(config1);
    const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: false});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    await tree.treeMutationsLocked$.pipe(filter((x) => !x), take(1)).toPromise();
    const node = tree.nodeTree.getNode([{idx: 1}]);
    const a = node.getItem().getStateStore().getState('a');
    const b = node.getItem().getStateStore().getState('b');
    expectDeepEqual(a, 33);
    expectDeepEqual(b, 44);
  });

});
