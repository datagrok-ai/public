import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {FuncCallMockAdapter} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallAdapters';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {TestScheduler} from 'rxjs/testing';
import {FuncCallIODescription} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';

category('ComputeUtils: Driver instance bridge', async () => {
  let testScheduler: TestScheduler;

  const io: FuncCallIODescription[] = [
    {direction: 'input', id: 'arg1', type: 'string', nullable: true},
    {direction: 'input', id: 'arg2', type: 'string', nullable: true},
    {direction: 'output', id: 'res', type: 'string', nullable: true},
  ];

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Init', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Init with values', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      bridge.setPreInitialData({
        initialRestrictions: {},
        initialValues: {
          arg1: 1,
          arg2: 2,
        },
      });
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: true});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: 1});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: 2});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Init with values and restrictions', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      bridge.setPreInitialData({
        initialRestrictions: {
          arg1: {
            assignedValue: 1,
            type: 'restricted',
          },
          arg2: {
            assignedValue: 3,
            type: 'info',
          },
        },
        initialValues: {
          arg1: 1,
          arg2: 2,
        },
      });
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: true});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {
        arg1: {
          assignedValue: 1,
          type: 'restricted',
        },
        arg2: {
          assignedValue: 3,
          type: 'info',
        },
      }});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: 1});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: 2});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });


  test('Init with previous state', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {
          arg1: {
            assignedValue: 10,
            type: 'restricted',
          },
        }, isOutputOutdated: true, initValues: true});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {
        a: {},
        b: {
          arg1: {
            assignedValue: 10,
            type: 'restricted',
          },
        },
      });
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: 10});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Init RO with previous state', async () => {
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
        }, isOutputOutdated: true, initValues: true});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {
        a: {},
        b: {
          arg1: {
            assignedValue: 10,
            type: 'restricted',
          },
        },
      });
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('SetState changes', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setState('arg1', 1);
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('abc', {a: {}, b: {}, c: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('abc', {a: undefined, b: undefined, c: 1});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('SetState additional states changes', async () => {
    testScheduler.run((helpers) => {
      const { expectObservable, cold } = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [{id: 'someState'}], false);
      cold('-a').subscribe(() => {
        bridge.init({ adapter, restrictions: {}, isOutputOutdated: true, initValues: false });
      });
      cold('--a').subscribe(() => {
        bridge.setState('arg1', 1);
      });
      cold('----a').subscribe(() => {
        bridge.setState('someState', {x: 1});
      });

      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('abc', { a: {}, b: {}, c: {} });
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', { a: true });
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', { a: {} });
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', { a: false });
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', { a: false, b: true });
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('abc', { a: undefined, b: undefined, c: 1 });
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', { a: undefined, b: undefined });
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', { a: undefined, b: undefined });
      expectObservable(bridge.getStateChanges('someState'), '^ 1000ms !').toBe('a---b', { a: undefined, b: {x: 1} });
    });
  });

  test('Edit state restrictions', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setState('arg1', 1, 'restricted');
      });
      cold('---a').subscribe(() => {
        bridge.editState('arg1', 2);
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('abc', {a: {}, b: {}, c: {
        arg1: {
          assignedValue: 1,
          type: 'restricted',
        },
      }});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('abcd', {a: undefined, b: undefined, c: 1, d: 2});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Edit state RO restrictions', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, true);
      adapter.editState('arg1', 1);
      const bridge = new FuncCallInstancesBridge(io, [], true);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setState('arg1', 2, 'restricted');
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('abc', {a: {}, b: {}, c: {
        arg1: {
          assignedValue: 2,
          type: 'restricted',
        },
      }});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('--a', {a: [
        'arg1',
        {
          'assignedValue': 2,
          'type': 'restricted',
        },
      ]});
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: 1});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Runnable single validator', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setValidation('arg1', 'mock', ({errors: ['mock error']}));
      });
      cold('---a').subscribe(() => {
        bridge.setValidation('arg1', 'mock', undefined);
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a-bc', {
        a: {},
        b: {
          'mock': {
            'arg1': {
              'errors': [
                'mock error'
              ],
            },
          },
        },
        c: {
          'mock': {
            'arg1': undefined,
          },
        },
      });
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('abcd', {a: false, b: true, c: false, d: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

  test('Runnable miltiple validators', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, [], false);
      cold('-a').subscribe(() => {
        bridge.init({adapter, restrictions: {}, isOutputOutdated: true, initValues: false});
      });
      cold('--a').subscribe(() => {
        bridge.setValidation('arg1', 'mock1', ({errors: ['mock error 1']}));
      });
      cold('---a').subscribe(() => {
        bridge.setValidation('arg1', 'mock2', ({errors: ['mock error 2']}));
      });
      cold('----a').subscribe(() => {
        bridge.setValidation('arg1', 'mock2', undefined);
      });
      cold('-----a').subscribe(() => {
        bridge.setValidation('arg1', 'mock1', {});
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a-bcde', {
        a: {},
        b: {
          'mock1': {
            'arg1': {
              'errors': [
                'mock error 1'
              ],
            },
          },
        },
        c: {
          'mock1': {
            'arg1': {
              'errors': [
                'mock error 1'
              ],
            },
          },
          'mock2': {
            'arg1': {
              'errors': [
                'mock error 2'
              ],
            },
          },
        },
        d: {
          'mock1': {
            'arg1': {
              'errors': [
                'mock error 1'
              ],
            },
          },
          'mock2': {
            'arg1': undefined,
          },
        },
        e: {
          'mock1': {
            'arg1': {},
          },
          'mock2': {
            'arg1': undefined,
          },
        },
      });
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('abc--d', {a: false, b: true, c: false, d: true});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });
});
