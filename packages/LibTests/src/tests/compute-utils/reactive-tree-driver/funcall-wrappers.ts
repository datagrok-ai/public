import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {FuncCallAdapter, FuncCallMockAdapter} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallAdapters';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {TestScheduler} from 'rxjs/testing';
import {map, take, takeUntil, toArray} from 'rxjs/operators';
import {Subject} from 'rxjs';

category('ComputeUtils: Driver mock wrapper', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  test('Set and get state', async () => {
    const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
    adapter.setState('state1', 1);
    const val = adapter.getState('state1');
    expectDeepEqual(val, 1);
  });

  test('Set and get df state', async () => {
    const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    adapter.setState('state1', df.clone());
    const val = adapter.getState('state1');
    expectDeepEqual(val, df);
  });

  test('Get state names', async () => {
    const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
    const names = adapter.getStateNames();
    expectDeepEqual(names, ['state1', 'state2']);
  });

  test('Track state changes', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
      cold('-a').subscribe(() => {
        adapter.setState('state1', 1);
      });
      cold('--a').subscribe(() => {
        adapter.setState('state2', 2);
      });
      expectObservable(adapter.getStateChanges('state1'), '^ 1000ms !').toBe('ab', {a: undefined, b: 1});
      expectObservable(adapter.getStateChanges('state2'), '^ 1000ms !').toBe('a-b', {a: undefined, b: 2});
    });
  });

  test('Ignore read only state changes', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], true);
      cold('-a').subscribe(() => {
        adapter.setState('state1', 1);
      });
      cold('--a').subscribe(() => {
        adapter.setState('state2', 2);
      });
      expectObservable(adapter.getStateChanges('state1'), '^ 1000ms !').toBe('a', {a: undefined});
      expectObservable(adapter.getStateChanges('state2'), '^ 1000ms !').toBe('a', {a: undefined});
    });
  });

  test('Track dataframe mutations', async () => {
    const df1 = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
    setTimeout(() => {
      adapter.setState('state1', df1);
    }, 1);
    setTimeout(() => {
      df1.set('col1', 1, 'val3');
    }, 2);
    const val = await adapter.getStateChanges('state1', true).pipe(map((x?: DG.DataFrame) =>x?.toJson()), take(3), toArray()).toPromise();
    expectDeepEqual(val, [
      null,
      [
        {
          'col1': 'val1',
        },
        {
          'col1': 'val2',
        },
      ],
      [
        {
          'col1': 'val1',
        },
        {
          'col1': 'val3',
        },
      ],
    ]);
  });

  test('Support mock runs', async () => {
    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const adapter = new FuncCallMockAdapter([{id: 'state1'}, {id: 'state2'}], false);
      expectObservable(adapter.getStateChanges('state2'), '^ 1000ms !').toBe('a 9ms b', {a: undefined, b: 1});
      adapter.run({'state2': 1}, 10).subscribe();
    });
  });
});

category('ComputeUtils: Driver FuncCall wrapper', async () => {
  before(async () => {});

  test('Set and get state', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestAdd2');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    adapter.setState('a', 1);
    const val = adapter.getState('a');
    expectDeepEqual(val, 1);
  });

  test('Set and get df state', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestDF1');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    adapter.setState('df', df.clone());
    const val = adapter.getState('df');
    expectDeepEqual(val, df);
  });

  test('Get state names', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestAdd2');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    const names = adapter.getStateNames();
    expectDeepEqual(names, ['a', 'b', 'res']);
  });

  test('Track state changes', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestAdd2');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    setTimeout(() => {
      adapter.setState('a', 1);
    }, 1);
    setTimeout(() => {
      adapter.setState('a', 2);
    }, 2);
    const val = await adapter.getStateChanges('a', true).pipe(take(3), toArray()).toPromise();
    expectDeepEqual(val, [
      null,
      1,
      2,
    ]);
  });

  test('Ignore read only state changes', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestAdd2');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, true);
    setTimeout(() => {
      adapter.setState('a', 1);
    }, 1);
    setTimeout(() => {
      adapter.setState('a', 2);
    }, 2);
    const end$ = new Subject<true>();
    setTimeout(() => {
      end$.next(true);
    }, 10);
    const val = await adapter.getStateChanges('a', true).pipe(takeUntil(end$), toArray()).toPromise();
    expectDeepEqual(val, [null]);
  });

  test('Track dataframe mutations', async () => {
    const df1 = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    const f: DG.Func = await grok.functions.eval('Libtests:TestDF1');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    setTimeout(() => {
      adapter.setState('df', df1);
    }, 1);
    setTimeout(() => {
      df1.set('col1', 1, 'val3');
    }, 2);
    const val = await adapter.getStateChanges('df', true).pipe(map((x?: DG.DataFrame) => x?.toJson()), take(3), toArray()).toPromise();
    expectDeepEqual(val, [
      null,
      [
        {
          'col1': 'val1',
        },
        {
          'col1': 'val2',
        },
      ],
      [
        {
          'col1': 'val1',
        },
        {
          'col1': 'val3',
        },
      ],
    ]);
  });

  test('Support run', async () => {
    const f: DG.Func = await grok.functions.eval('Libtests:TestAdd2');
    const call = f.prepare({});
    const adapter = new FuncCallAdapter(call, false);
    setTimeout(() => {
      adapter.setState('a', 1);
      adapter.setState('b', 2);
    }, 1);
    setTimeout(async () => {
      await adapter.run().toPromise();
    }, 2);
    const end$ = new Subject<true>();
    setTimeout(() => {
      end$.next(true);
    }, 10);
    const val = await adapter.getStateChanges('res', true).pipe(takeUntil(end$), toArray()).toPromise();
    expectDeepEqual(val, [null, 3]);
  });
});
