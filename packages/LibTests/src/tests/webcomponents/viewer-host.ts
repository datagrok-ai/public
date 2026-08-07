import {category, test, before} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {of, Observable} from 'rxjs';
import {TestScheduler} from 'rxjs/testing';
import {ViewerHost, ViewerLike} from '@datagrok-libraries/webcomponents/src/Viewer/ViewerHost';
import {createTestScheduler} from '../../test-utils';

// Virtual-time tests of the dg-viewer state core: type/dataFrame/options set within one
// tick coalesce into a single transition; a transition applies the frame first and the
// final options second — stale options are never applied (re-asserting them against a new
// frame crashes the Dart legend refresh, see the split select/deselect bug).

class MockViewer implements ViewerLike {
  private frame: any;
  constructor(
    readonly id: string,
    readonly type: string,
    private readonly log: string[],
    frame: any,
  ) {
    this.frame = frame;
  }

  get dataFrame() {
    return this.frame;
  }

  set dataFrame(df: any) {
    this.frame = df;
    this.log.push(`frame:${this.id}:${df.name}`);
  }

  setOptions(options: Record<string, any>) {
    this.log.push(`options:${this.id}:${JSON.stringify(options)}`);
  }
}

type Factory = (type: string, df: any) => Observable<MockViewer>;

function setup(createOverride?: (mk: (type: string, df: any) => MockViewer) => Factory) {
  const log: string[] = [];
  let counter = 0;
  const mk = (type: string, df: any) => {
    const viewer = new MockViewer(`v${++counter}`, type, log, df);
    log.push(`create:${viewer.id}:${type}:${df.name}`);
    return viewer;
  };
  const factory: Factory = createOverride ?
    createOverride(mk) :
    (type, df) => of(mk(type, df));
  const host = new ViewerHost<MockViewer>(factory);
  return {log, host};
}

const df1 = {name: 'df1'};
const df2 = {name: 'df2'};

category('WebComponents: ViewerHost', () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('creates a viewer once type and frame are both set', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      flush();
      expectDeepEqual(log, [], {prefix: 'type alone must not create'});
      host.dfSetted$.next(df1);
      flush();
      expectDeepEqual(log, ['create:v1:grid:df1', 'options:v1:{}']);
    });
  });

  test('creation applies the options of the same tick', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      flush();
      expectDeepEqual(log, ['create:v1:grid:df1', 'options:v1:{"a":1}']);
    });
  });

  test('property order within a tick is irrelevant', async () => {
    const run = (apply: (host: ViewerHost<MockViewer>) => void) => {
      let result: string[] = [];
      testScheduler.run(({flush}) => {
        const {log, host} = setup();
        apply(host);
        flush();
        result = log;
      });
      return result;
    };
    const frameFirst = run((host) => {
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      host.typeSetted$.next('grid');
    });
    const optionsFirst = run((host) => {
      host.optionsSetted$.next({a: 1});
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
    });
    expectDeepEqual(frameFirst, ['create:v1:grid:df1', 'options:v1:{"a":1}']);
    expectDeepEqual(optionsFirst, frameFirst);
  });

  test('multiple frame updates in one tick coalesce to the last', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.dfSetted$.next(df2);
      host.dfSetted$.next(df1);
      flush();
      expectDeepEqual(log, ['create:v1:grid:df1', 'options:v1:{}']);
    });
  });

  test('frame-only change swaps in place and reapplies options', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      flush();
      host.dfSetted$.next(df2);
      flush();
      expectDeepEqual(log.slice(2), ['frame:v1:df2', 'options:v1:{"a":1}']);
    });
  });

  test('frame and options together apply the frame first, final options once', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('line');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({split: ['species']});
      flush();
      host.dfSetted$.next(df2);
      host.optionsSetted$.next({split: []});
      flush();
      // the old {split: ['species']} must never be applied to df2 — that pairing
      // is exactly what crashed the Dart legend refresh
      expectDeepEqual(log.slice(2), ['frame:v1:df2', 'options:v1:{"split":[]}']);
    });
  });

  test('intermediate options within a tick are never applied', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('line');
      host.dfSetted$.next(df1);
      flush();
      host.dfSetted$.next(df2);
      host.optionsSetted$.next({step: 'intermediate'});
      host.optionsSetted$.next({step: 'final'});
      flush();
      expectDeepEqual(log.slice(2), ['frame:v1:df2', 'options:v1:{"step":"final"}']);
    });
  });

  test('options-only change applies once', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      flush();
      host.optionsSetted$.next({a: 1});
      flush();
      expectDeepEqual(log.slice(2), ['options:v1:{"a":1}']);
    });
  });

  test('referentially new but equal options are not reapplied', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      flush();
      host.optionsSetted$.next({a: 1});
      flush();
      expectDeepEqual(log.length, 2, {prefix: 'no effects for equal options'});
    });
  });

  test('in-place mutation of the same options object is detected', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      const options: Record<string, any> = {a: 1};
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next(options);
      flush();
      options.a = 2;
      host.optionsSetted$.next(options);
      flush();
      expectDeepEqual(log.slice(2), ['options:v1:{"a":2}']);
    });
  });

  test('re-setting the same frame reference is a no-op', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      flush();
      host.dfSetted$.next(df1);
      flush();
      expectDeepEqual(log.length, 2, {prefix: 'no effects for same frame'});
    });
  });

  test('type change recreates with the current frame and options', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      flush();
      host.typeSetted$.next('line');
      flush();
      expectDeepEqual(log.slice(2), ['create:v2:line:df1', 'options:v2:{"a":1}']);
      expectDeepEqual(host.viewer$.value!.id, 'v2');
    });
  });

  test('clearing the frame detaches the viewer', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      const attachLog: (string | undefined)[] = [];
      host.viewer$.subscribe((viewer) => attachLog.push(viewer?.id));
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      flush();
      host.dfSetted$.next(undefined);
      flush();
      expectDeepEqual(attachLog, [undefined, 'v1', undefined]);
      expectDeepEqual(log.length, 2, {prefix: 'detaching produces no viewer effects'});
    });
  });

  test('frame re-set after detach creates a fresh viewer', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      flush();
      host.dfSetted$.next(undefined);
      flush();
      host.dfSetted$.next(df2);
      flush();
      expectDeepEqual(log.slice(2), ['create:v2:grid:df2', 'options:v2:{}']);
    });
  });

  test('slow creation queues the next transition in order', async () => {
    testScheduler.run(({cold, flush}) => {
      const {log, host} = setup((mk) => (type, df) =>
        new Observable<MockViewer>((subscriber) => {
          const sub = cold('--x', {x: null}).subscribe(() => {
            subscriber.next(mk(type, df));
            subscriber.complete();
          });
          return () => sub.unsubscribe();
        }));
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      // a frame update arriving while the creation is in flight
      cold('-x').subscribe(() => host.dfSetted$.next(df2));
      flush();
      expectDeepEqual(log, [
        'create:v1:grid:df1', 'options:v1:{}',
        'frame:v1:df2', 'options:v1:{}',
      ]);
    });
  });

  test('creation applies options that changed during the flight', async () => {
    testScheduler.run(({cold, flush}) => {
      const {log, host} = setup((mk) => (type, df) =>
        new Observable<MockViewer>((subscriber) => {
          const sub = cold('--x', {x: null}).subscribe(() => {
            subscriber.next(mk(type, df));
            subscriber.complete();
          });
          return () => sub.unsubscribe();
        }));
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.optionsSetted$.next({a: 1});
      cold('-x').subscribe(() => host.optionsSetted$.next({a: 2}));
      flush();
      expectDeepEqual(log, ['create:v1:grid:df1', 'options:v1:{"a":2}']);
    });
  });

  test('frameApplied$ fires after an in-place swap only', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      const applied: string[] = [];
      host.frameApplied$.subscribe((viewer) => applied.push(`${viewer.id}:${viewer.dataFrame.name}`));
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      flush();
      expectDeepEqual(applied, [], {prefix: 'creation must not fire frameApplied$'});
      host.optionsSetted$.next({a: 1});
      flush();
      expectDeepEqual(applied, [], {prefix: 'options-only must not fire frameApplied$'});
      host.dfSetted$.next(df2);
      flush();
      expectDeepEqual(applied, ['v1:df2']);
      expectDeepEqual(log.slice(3), ['frame:v1:df2', 'options:v1:{"a":1}'],
        {prefix: 'frameApplied$ fires after the swap effects'});
    });
  });

  test('destroy cancels a pending transition', async () => {
    testScheduler.run(({flush}) => {
      const {log, host} = setup();
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      host.destroy();
      flush();
      expectDeepEqual(log, []);
    });
  });

  test('destroy during a slow creation prevents the attach', async () => {
    testScheduler.run(({cold, flush}) => {
      const {log, host} = setup((mk) => (type, df) =>
        new Observable<MockViewer>((subscriber) => {
          const sub = cold('----x', {x: null}).subscribe(() => {
            subscriber.next(mk(type, df));
            subscriber.complete();
          });
          return () => sub.unsubscribe();
        }));
      const attachLog: (string | undefined)[] = [];
      host.viewer$.subscribe((viewer) => attachLog.push(viewer?.id));
      host.typeSetted$.next('grid');
      host.dfSetted$.next(df1);
      cold('--x').subscribe(() => host.destroy());
      flush();
      expectDeepEqual(log, [], {prefix: 'creation must be cancelled'});
      expectDeepEqual(attachLog, [undefined]);
    });
  });
});
