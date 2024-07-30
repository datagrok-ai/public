import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import cloneDeepWith from 'lodash.clonedeepwith';
import {Observable, defer, from} from 'rxjs';
import {HandlerBase} from './config/PipelineConfiguration';

export function callHandler<R, P = any>(handler: HandlerBase<P, R>, params: P): Observable<R> {
  if (typeof handler === 'string') {
    return defer(async () => {
      const f: DG.Func = await grok.functions.eval(handler);
      const call = f.prepare({params});
      await call.call();
      const res = call.outputs as R;
      return res;
    });
  } else {
    return defer(() => {
      const res = handler(params);
      if (res instanceof Observable || res instanceof Promise)
        return res;
      else
        return from([res]);
    });
  }
}

export function cloneConfig<T>(config: T): T {
  return cloneDeepWith(config, (val) => {
    if (typeof val === 'function')
      return val;
  });
}

export function indexFromEnd<T>(arr: Readonly<T[]>, offset = 0): T | undefined {
  return arr[arr.length - offset - 1];
}

const PIPELINE_DEBUG = true;

export function debuglog(...args: any[]) {
  if (PIPELINE_DEBUG)
    console.log(...args);
}
