import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {FuncCallIODescription} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {LinksState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {FuncCallMockAdapter} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallAdapters';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver readonly mock adapter', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('setState is ignored in readonly mode', async () => {
    const adapter = new FuncCallMockAdapter([{id: 'a'}, {id: 'b'}], true);
    adapter.setState('a', 5);
    adapter.setState('b', 10);
    expectDeepEqual(adapter.getState('a'), undefined);
    expectDeepEqual(adapter.getState('b'), undefined);
  });

  test('editState bypasses readonly guard', async () => {
    const adapter = new FuncCallMockAdapter([{id: 'a'}, {id: 'b'}], true);
    adapter.editState('a', 5);
    adapter.editState('b', 10);
    expectDeepEqual(adapter.getState('a'), 5);
    expectDeepEqual(adapter.getState('b'), 10);
  });

  test('getStateChanges does not emit on setState', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter([{id: 'a'}], true);
      cold('-a').subscribe(() => {
        adapter.setState('a', 5);
      });
      expectObservable(adapter.getStateChanges('a'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('getStateChanges emits on editState', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter([{id: 'a'}], true);
      cold('-a').subscribe(() => {
        adapter.editState('a', 5);
      });
      expectObservable(adapter.getStateChanges('a'), '^ 1000ms !').toBe('ab', {a: undefined, b: 5});
    });
  });
});


category('ComputeUtils: Driver readonly bridge', async () => {
  let testScheduler: TestScheduler;

  const io: FuncCallIODescription[] = [
    {direction: 'input', id: 'arg1', type: 'string', nullable: true},
    {direction: 'input', id: 'arg2', type: 'string', nullable: true},
    {direction: 'output', id: 'res', type: 'string', nullable: true},
  ];

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('setState propagates restrictions but not adapter value', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, true);
      const bridge = new FuncCallInstancesBridge(io, [], true);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setState('arg1', 'hello', 'restricted');
      });
      // restriction is recorded
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('abc', {
        a: {},
        b: {},
        c: {
          arg1: {
            assignedValue: 'hello',
            type: 'restricted',
          },
        },
      });
      // restriction update is emitted (RO path fires inputRestrictionsUpdates$)
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('--a', {
        a: [
          'arg1',
          {
            assignedValue: 'hello',
            type: 'restricted',
          },
        ],
      });
      // but adapter value stays undefined — setState was blocked
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('setToConsistent is no-op in readonly mode', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, true);
      const bridge = new FuncCallInstancesBridge(io, [], true);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {
          arg1: {
            assignedValue: 10,
            type: 'restricted',
          },
        }, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setToConsistent('arg1');
      });
      // restrictions preserved
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {
        a: {},
        b: {
          arg1: {
            assignedValue: 10,
            type: 'restricted',
          },
        },
      });
      // adapter value stays undefined — setToConsistent skipped due to isReadonly
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('editState works through bridge in readonly mode', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, true);
      const bridge = new FuncCallInstancesBridge(io, [], true);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.editState('arg1', 'edited');
      });
      // editState bypasses readonly — value changes
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('abc', {a: undefined, b: undefined, c: 'edited'});
    });
  });

  test('Validators work in readonly bridge', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, true);
      const bridge = new FuncCallInstancesBridge(io, [], true);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setValidation('arg1', 'validator1', ({errors: ['required']}));
      });
      // validations work in RO — they write to validations$, not to the adapter
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a-b', {
        a: {},
        b: {
          validator1: {
            arg1: {
              errors: ['required'],
            },
          },
        },
      });
      // isRunable reacts to validation errors
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('abc', {a: false, b: true, c: false});
    });
  });
});


category('ComputeUtils: Driver readonly tree', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Data link does not update output state in readonly tree', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, isReadonly: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().editState('b', 1);
        link.trigger();
      });
      // output state stays undefined — setState on RO adapter is blocked
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('Restriction link propagates restrictions in readonly tree', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
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
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, isReadonly: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().editState('b', 1);
      });
      // restrictions propagate even in RO
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).inputRestrictions$).toBe('a b', {
        a: {},
        b: {
          'a': {
            'type': 'restricted',
            'assignedValue': 1,
          },
        },
      });
      // but state itself stays undefined
      expectObservable(outNode.getItem().getStateStore().getStateChanges('a'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('Validator link runs in readonly tree', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step1/a',
        type: 'validator',
        handler({controller}) {
          controller.setValidation('out1', ({warnings: [{description: 'readonly warn'}]}));
          return;
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, isReadonly: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      const ls = new LinksState();
      const [link] = ls.createStateLinks(tree.nodeTree);
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      link.wire(tree.nodeTree);
      cold('-a').subscribe(() => {
        link.trigger();
      });
      // validators write to validations$, not to adapter — works in RO
      expectObservable((inNode.getItem().getStateStore() as FuncCallInstancesBridge).validations$).toBe('a b', {
        a: {},
        b: {
          [link.uuid]: {
            'a': {
              'warnings': [
                {
                  'description': 'readonly warn',
                },
              ],
            },
          },
        },
      });
    });
  });

  test('Meta link runs in readonly tree', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/b',
        to: 'out1:step2/a',
        type: 'meta',
        handler({controller}) {
          controller.setViewMeta('out1', {key: 'ro-meta'});
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, isReadonly: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outNode = tree.nodeTree.getNode([{idx: 1}]);
      cold('-a').subscribe(() => {
        inNode.getItem().getStateStore().editState('b', 1);
      });
      // meta writes to metaStates, not to adapter — works in RO
      expectObservable((outNode.getItem().getStateStore() as FuncCallInstancesBridge).meta.a
      ).toBe('ab', {
        a: undefined,
        b: {
          'key': 'ro-meta',
        },
      });
    });
  });

  test('Consistency shows inconsistent after editState in readonly tree', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/a',
        to: 'out1:step2/a',
        handler({controller}) {
          controller.setAll('out1', 2, 'restricted');
        },
      }],
    };
    const pconf = await getProcessedConfig(config);

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
});
