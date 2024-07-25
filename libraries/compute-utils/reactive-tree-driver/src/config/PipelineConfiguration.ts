import {Observable} from 'rxjs';
import {IRuntimeController} from '../RuntimeController';
import {ItemId, NqName, ItemPath, InputState, ItemType, ItemPathArray} from '../data/common-types';
import {PipelineInstanceConfig, StepParallelInitialConfig, StepSequentialInitialConfig} from './PipelineInstance';

//
// Pipeline public configuration
//

export type StateItem = {
  id: ItemId;
  type?: string;
};

export type PipelineGlobalId = {
  globalId: string,
}

export type PipelineSelfRef = {
  type: 'selfRef',
  selfRef: string,
}

// handlers

export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R>) | NqName;
export type Handler = HandlerBase<{ controller: IRuntimeController }, void>;
export type SelectorKeyExtractor = HandlerBase<{ controller: IRuntimeController }, string>;
export type PipelineProvider = HandlerBase<{ version?: string }, PipelineConfigurationInitial & PipelineGlobalId>;

// link-like

export type PipelineLinkConfigurationBase<P> = {
  from: P;
  to: P;
  dataFrameMutations?: boolean | P;
  inputState?: InputState | InputState[];
  allTypeNodes?: boolean | P;
  handler?: Handler;
}

export type PipelineLinkConfiguration<P> = {
  id: ItemId;
} & PipelineLinkConfigurationBase<P>;

export type PipelineDynamicLinkConfiguration<P> = {
  // TODO: better spec
  id: ItemId;
} & PipelineLinkConfigurationBase<P>;

export type PipelineHookConfiguration<P> = {
  id: ItemId;
  handler: Handler;
} & Partial<PipelineLinkConfigurationBase<P>>;

export type PipelineActionConfiguraion<P> = {
  id: ItemId;
  path: P;
  position: ActionPositions;
  friendlyName?: string;
} & PipelineLinkConfigurationBase<P>;

const actionPositions = ['buttons', 'menu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

// hooks config

export type PipelineHooks<P> = {
  onInit?: PipelineHookConfiguration<P>[];
  beforeLoadFuncCall?: PipelineHookConfiguration<P>[];
  afterLoadFuncCall?: PipelineHookConfiguration<P>[];
  beforeInputFormRender?: PipelineHookConfiguration<P>[];
  afterInputFormRender?: PipelineHookConfiguration<P>[];
  beforeViewerRender?: PipelineHookConfiguration<P>[];
  afterViewerRender?: PipelineHookConfiguration<P>[];
  beforeLoadRun?: PipelineHookConfiguration<P>[];
  afterLoadRun?: PipelineHookConfiguration<P>[];
  beforeSaveRun?: PipelineHookConfiguration<P>[];
  afterSaveRun?: PipelineHookConfiguration<P>[];
  onClose?: PipelineHookConfiguration<P>[];
};

export type HooksSpec = {
  path: ItemPathArray;
  hooks: PipelineHooks<ItemPathArray>;
}

// static steps config

export type PipelineStepConfiguration<S> = {
  id: ItemId;
  nqName: NqName;
  friendlyName?: string;
  io?: S;
};

export type PipelineConfigurationBase<P> = {
  id: ItemId;
  nqName?: NqName;
  friendlyName?: string;
  hooks?: PipelineHooks<P>;
  actions?: PipelineActionConfiguraion<P>[];
  states?: StateItem[];
};

// fixed pipeline

export type AbstractPipelineStaticConfiguration<P, S, R> = {
  links?: PipelineLinkConfiguration<P>[];
  steps: (PipelineStepConfiguration<S> | AbstractPipelineConfiguration<P, S, R> | R)[];
  type: 'static';
} & PipelineConfigurationBase<P>;

// parallel pipeline

export type ParallelItemContext<P> = {
  id: ItemId;
  allowAdding: boolean;
  selectorPath?: P;
  selectorExtractor?: SelectorKeyExtractor;
};

export type PipelineParallelItem<P, S, R> = (AbstractPipelineConfiguration<P, S, R> | R) & ParallelItemContext<P>;

export type AbstractPipelineParallelConfiguration<P, S, R> = {
  initialSteps?: StepParallelInitialConfig[];
  stepTypes: PipelineParallelItem<P, S, R>[];
  type: 'parallel';
} & PipelineConfigurationBase<P>;

// sequential pipeline

export type SequentialItemContext = {
  id: ItemId;
  allowAdding: boolean;
};

export type PipelineSequentialItem<P, S, R> = (AbstractPipelineConfiguration<P, S, R> | R) & SequentialItemContext;

export type AbstractPipelineSequentialConfiguration<P, S, R> = {
  initialSteps?: StepSequentialInitialConfig[];
  stepTypes: PipelineSequentialItem<P, S, R>[];
  links?: PipelineDynamicLinkConfiguration<P>[];
  type: 'sequential';
} & PipelineConfigurationBase<P>;

// pipeline config

export type AbstractPipelineConfiguration<P, S, R> =
AbstractPipelineStaticConfiguration<P, S, R> |
AbstractPipelineParallelConfiguration<P, S, R> |
AbstractPipelineSequentialConfiguration<P, S, R>;

export type PipelineRefIntial = {
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<ItemPath | ItemPath[], never, PipelineRefIntial>;
export type PipelineConfigurationParallelIntial = AbstractPipelineParallelConfiguration<ItemPath | ItemPath[], never, PipelineRefIntial>;
export type PipelineConfigurationSequentialInitial = AbstractPipelineSequentialConfiguration<ItemPath | ItemPath[], never, PipelineRefIntial>;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationParallelIntial | PipelineConfigurationSequentialInitial | PipelineRefIntial;

// extrenal instance config

export type ExternalInitialConfig = {
  provider: PipelineProvider | NqName;
  version?: string;
  config: PipelineInstanceConfig;
}
