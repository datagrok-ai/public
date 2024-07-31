import {Observable} from 'rxjs';
import {IRuntimeController} from '../RuntimeController';
import {ItemId, NqName, ItemPath, InputState, ItemPathArray} from '../data/common-types';
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
  // TODO: create spec
  id: ItemId;
} & PipelineLinkConfigurationBase<P>;

export type PipelineHookConfiguration<P> = {
  id: ItemId;
  handler: Handler;
} & Partial<PipelineLinkConfigurationBase<P>>;

export type PipelineActionConfiguraion<P> = {
  id: ItemId;
  path: P;
  friendlyName?: string;
  menuCategory?: string;
} & PipelineLinkConfigurationBase<P>;

export type StepActionConfiguraion<P> = {
  id: ItemId;
  path: P;
  position: ActionPositions;
  friendlyName?: string;
  menuCategory?: string;
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

export type PipelineStepConfiguration<P, S> = {
  id: ItemId;
  nqName: NqName;
  friendlyName?: string;
  io?: S;
  actions?: StepActionConfiguraion<P>[];
};

export type PipelineConfigurationBase<P, G> = {
  id: ItemId;
  globalId?: G;
  nqName?: NqName;
  friendlyName?: string;
  hooks?: PipelineHooks<P>;
  actions?: PipelineActionConfiguraion<P>[];
  states?: StateItem[];
};

// fixed pipeline

export type PipelineStaticItem<P, S, R, G> = {
  id: ItemId;
} & (PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R, G> | R);

export type AbstractPipelineStaticConfiguration<P, S, R, G> = {
  links?: PipelineLinkConfiguration<P>[];
  steps: PipelineStaticItem<P, S, R, G>[];
  type: 'static';
} & PipelineConfigurationBase<P, G>;

// parallel pipeline

export type ParallelItemContext<P> = {
  id: ItemId;
  allowAdding: boolean;
  selectorPath?: P;
  selectorExtractor?: SelectorKeyExtractor;
};

export type PipelineParallelItem<P, S, R, G> = (AbstractPipelineConfiguration<P, S, R, G> | R) & ParallelItemContext<P>;

export type AbstractPipelineParallelConfiguration<P, S, R, G> = {
  initialSteps?: StepParallelInitialConfig[];
  stepTypes: PipelineParallelItem<P, S, R, G>[];
  type: 'parallel';
} & PipelineConfigurationBase<P, G>;

// sequential pipeline

export type SequentialItemContext = {
  id: ItemId;
  allowAdding: boolean;
};

export type PipelineSequentialItem<P, S, R, G> = (AbstractPipelineConfiguration<P, S, R, G> | R) & SequentialItemContext;

export type AbstractPipelineSequentialConfiguration<P, S, R, G> = {
  initialSteps?: StepSequentialInitialConfig[];
  stepTypes: PipelineSequentialItem<P, S, R, G>[];
  links?: PipelineDynamicLinkConfiguration<P>[];
  type: 'sequential';
} & PipelineConfigurationBase<P, G>;

// pipeline config

export type AbstractPipelineConfiguration<P, S, R, G> =
AbstractPipelineStaticConfiguration<P, S, R, G> |
AbstractPipelineParallelConfiguration<P, S, R, G> |
AbstractPipelineSequentialConfiguration<P, S, R, G>;

export type PipelineRefInitial = {
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<ItemPath | ItemPath[], never, PipelineRefInitial, never>;
export type PipelineConfigurationParallelInitial = AbstractPipelineParallelConfiguration<ItemPath | ItemPath[], never, PipelineRefInitial, never>;
export type PipelineConfigurationSequentialInitial = AbstractPipelineSequentialConfiguration<ItemPath | ItemPath[], never, PipelineRefInitial, never>;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationParallelInitial | PipelineConfigurationSequentialInitial | PipelineRefInitial;

// extrenal instance config

export type ExternalInitialConfig = {
  provider: PipelineProvider | NqName;
  version?: string;
  config: PipelineInstanceConfig;
}
