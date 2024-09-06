import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { category, test, before } from '@datagrok-libraries/utils/src/test';
import { FuncCallMockAdapter } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallAdapters';
import { expectDeepEqual } from '@datagrok-libraries/utils/src/expect';
import {TestScheduler} from 'rxjs/testing';
import { map, take, takeUntil, toArray } from 'rxjs/operators';
import { Subject } from 'rxjs';
import { FuncallStateItem } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import { FuncCallInstancesBridge } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';

category('ComputeUtils: Driver instance bridge', async () => {
  let testScheduler: TestScheduler;

  const io: FuncallStateItem[] = [
    {direction: 'input', id: 'arg1', type: 'string'},
    {direction: 'input', id: 'arg2', type: 'string'},
    {direction: 'output', id: 'res', type: 'string'},
  ];

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Init', async() => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter(io, false);
      const bridge = new FuncCallInstancesBridge(io, false);
      cold('-a').subscribe(() => {
        bridge.init({ adapter, restrictions: {}, isOutputOutdated: true, initValues: false });
      });
      expectObservable(bridge.inputRestrictions$, '^ 1000ms !').toBe('ab', {a: {}, b: {}});
      expectObservable(bridge.isOutputOutdated$, '^ 1000ms !').toBe('a', {a: true});
      expectObservable(bridge.validations$, '^ 1000ms !').toBe('a', {a: {}});
      expectObservable(bridge.isRunning$, '^ 1000ms !').toBe('a', {a: false});
      expectObservable(bridge.isRunable$, '^ 1000ms !').toBe('ab', {a: false, b: false});
      expectObservable(bridge.inputRestrictionsUpdates$, '^ 1000ms !').toBe('');
      expectObservable(bridge.getStateChanges('arg1'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('arg2'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
      expectObservable(bridge.getStateChanges('res'), '^ 1000ms !').toBe('ab', {a: undefined, b: undefined});
    });
  });

});
