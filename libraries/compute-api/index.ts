/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

type ConstructorTypeOf<T> = new (...args:any[]) => T;

declare global {
  interface Window {
    compute: {
      CFV: ConstructorTypeOf<CustomFunctionView>,
    },
  }
}

export type * from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';

import {
  CustomFunctionView, createCFV
} from './src/views';
export {
  createCFV
};

export async function initComputeApi() {
  const initFunc = DG.Func.find({package: 'Compute', name: 'init'})[0];
  if (initFunc)
    await initFunc.prepare().call();
}
