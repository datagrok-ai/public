/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';

type ConstructorTypeOf<T> = new (...args:any[]) => T;

declare global {
  interface Window {
    compute: {
      CFV: ConstructorTypeOf<CustomFunctionView>,
    },
  }
}

export type * from '@datagrok-libraries/compute-utils/reactive-tree-driver/index';

export async function initComputeApi() {
  const initFunc = DG.Func.find({package: 'Compute', name: 'init'})[0];
  if (initFunc)
    await initFunc.prepare().call();
}

export async function startWorkflow<T=any>(nqName: string, version: string, instanceConfig: PipelineInstanceConfig): Promise<T> {
  const startFunc = DG.Func.find({package: 'Compute2', name: 'StartWorkflow'})[0];
  const call = startFunc.prepare({nqName, version, instanceConfig});
  await call.call();
  return call.getOutputParamValue();
}

import type {
  CustomFunctionView,
} from '@datagrok-libraries/compute-utils';

export function createCFV(
  ...args: ConstructorParameters<typeof CustomFunctionView>
): CustomFunctionView {
  return new window.compute.CFV(...args);
};
