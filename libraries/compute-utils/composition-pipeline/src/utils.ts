import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {keyToPath} from './config/config-processing-utils';
import {Observable, defer, from} from 'rxjs';

const PIPELINE_DEBUG = true;

export function path(strings: TemplateStringsArray): string[] {
  const str = strings[0];
  return keyToPath(str);
}

type Result<R> = R | Observable<R> | Promise<R>;
type Args = Record<string, any>;

export function callHandler(fn: string | ((args: any) => Result<Args>), params: Args): Observable<Args> {
  if (typeof fn === 'string') {
    return defer(async () => {
      const f: DG.Func = await grok.functions.eval(fn);
      const call = f.prepare({params});
      await call.call();
      const res = call.outputs;
      return res;
    });
  } else {
    return defer(() => {
      const res = fn(params);
      if (res instanceof Observable || res instanceof Promise)
        return res;
      else
        return from([res]);
    });
  }
}

export function debuglog(...args: any[]) {
  if (PIPELINE_DEBUG)
    console.log(...args);
}
