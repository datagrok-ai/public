import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {keyToPath} from './config-processing-utils';

const PIPELINE_DEBUG = true;

export function path(strings: TemplateStringsArray): string[] {
  const str = strings[0];
  return keyToPath(str);
}

export async function callHandler(fn: string | Function, params: Record<string, any>) {
  if (typeof fn === 'string') {
    const f: DG.Func = await grok.functions.eval(fn);
    const call = f.prepare({params});
    await call.call();
    return call.outputs;
  } else {
    const res = await fn(params);
    return res;
  }
}

export function debuglog(...args: any[]) {
  if (PIPELINE_DEBUG)
    console.log(...args);
}
