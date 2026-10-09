import {
  PipelineState,
  PipelineInstanceRuntimeData,
  PipelineStateStatic,
  PipelineStateDynamic,
  StepDynamicDescription,
  StepFunCallState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

export function mockFuncCall(uuid: string, opts?: {isReadonly?: boolean}): StepFunCallState {
  return {
    type: 'funccall',
    uuid,
    configId: uuid,
    isReadonly: opts?.isReadonly ?? false,
  };
}

const defaultPipelineRuntimeData: PipelineInstanceRuntimeData = {
  actions: undefined,
  disableHistory: false,
  customExports: undefined,
  forceNavigate: false,
  compactView: false,
};

export function mockStaticPipeline(
  uuid: string,
  steps: PipelineState[],
  opts?: {
    isReadonly?: boolean; isActionStep?: boolean; forceNavigate?: boolean; compactView?: boolean;
    nqName?: string; disableHistory?: boolean; actions?: ViewAction[];
  },
): PipelineStateStatic<StepFunCallState, PipelineInstanceRuntimeData> {
  return {
    type: 'static',
    uuid,
    configId: uuid,
    friendlyName: undefined,
    description: undefined,
    version: undefined,
    nqName: opts?.nqName,
    isReadonly: opts?.isReadonly ?? false,
    steps,
    isActionStep: opts?.isActionStep,
    ...defaultPipelineRuntimeData,
    forceNavigate: opts?.forceNavigate ?? false,
    compactView: opts?.compactView ?? false,
    disableHistory: opts?.disableHistory ?? false,
    actions: opts?.actions,
  };
}

export function mockDynamicPipeline(
  uuid: string,
  steps: PipelineState[],
  opts?: {
    isReadonly?: boolean; forceNavigate?: boolean; type?: 'dynamic' | 'parallel' | 'sequential';
    stepTypes?: StepDynamicDescription[]; actions?: ViewAction[];
  },
): PipelineStateDynamic<StepFunCallState, PipelineInstanceRuntimeData> {
  return {
    type: opts?.type ?? 'dynamic',
    uuid,
    configId: uuid,
    friendlyName: undefined,
    description: undefined,
    version: undefined,
    nqName: undefined,
    isReadonly: opts?.isReadonly ?? false,
    steps,
    stepTypes: opts?.stepTypes ?? [],
    ...defaultPipelineRuntimeData,
    forceNavigate: opts?.forceNavigate ?? false,
    actions: opts?.actions,
  };
}
