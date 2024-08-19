import {Observable} from 'rxjs';
import {IRuntimeController} from '../IRuntimeController';
import {ItemId, NqName, FromSpec, ToSpec, RestrictionType} from '../data/common-types';
import {StepParallelInitialConfig, StepSequentialInitialConfig} from './PipelineInstance';

//
// Pipeline public configuration
//

export type StateItem = {
  id: ItemId;
  type?: string;
}

export type LoadedPipelineToplevelNode = {
  id: string,
  nqName: string,
}

export type PipelineSelfRef = {
  type: 'selfRef',
  selfRef: string,
  id: string,
}

// handlers

export type LoadedPipeline = (PipelineConfigurationStaticInitial | PipelineConfigurationParallelInitial | PipelineConfigurationSequentialInitial) & LoadedPipelineToplevelNode;

export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R>) | NqName;
export type Handler = HandlerBase<{ controller: IRuntimeController }, void>;
export type SelectorKeyExtractor = HandlerBase<{ controller: IRuntimeController }, string>;
export type PipelineProvider = HandlerBase<{ version?: string }, LoadedPipeline>;

// link-like

export type PipelineLinkConfigurationBase<F, T> = {
  from: F;
  to: T;
  dataFrameMutations?: boolean | string[] | Record<string, string>;
  inputState?: RestrictionType | RestrictionType[] | Record<string, RestrictionType>;
  handler?: Handler;
}

export type PipelineLinkConfiguration<F, T> = {
  id: ItemId;
} & PipelineLinkConfigurationBase<F, T>;

export type PipelineDynamicLinkConfiguration<F, T> = {
  // TODO: create spec
  id: ItemId;
} & PipelineLinkConfigurationBase<F, T>;

export type PipelineHookConfiguration<F, T> = {
  id: ItemId;
  handler: Handler;
} & Partial<PipelineLinkConfigurationBase<F, T>>;

export type PipelineActionConfiguraion<F, T> = {
  id: ItemId;
  path: F;
  friendlyName?: string;
  menuCategory?: string;
} & PipelineLinkConfigurationBase<F, T>;

export type StepActionConfiguraion<F, T> = {
  id: ItemId;
  path: F;
  position: ActionPositions;
  friendlyName?: string;
  menuCategory?: string;
} & PipelineLinkConfigurationBase<F, T>;


const actionPositions = ['buttons', 'menu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

// hooks config

export type PipelineHooks<F, T> = {
  onInit?: PipelineHookConfiguration<F, T>[];
  beforeLoadFuncCall?: PipelineHookConfiguration<F, T>[];
  afterLoadFuncCall?: PipelineHookConfiguration<F, T>[];
  beforeInputFormRender?: PipelineHookConfiguration<F, T>[];
  afterInputFormRender?: PipelineHookConfiguration<F, T>[];
  beforeViewerRender?: PipelineHookConfiguration<F, T>[];
  afterViewerRender?: PipelineHookConfiguration<F, T>[];
  beforeLoadRun?: PipelineHookConfiguration<F, T>[];
  afterLoadRun?: PipelineHookConfiguration<F, T>[];
  beforeSaveRun?: PipelineHookConfiguration<F, T>[];
  afterSaveRun?: PipelineHookConfiguration<F, T>[];
  onClose?: PipelineHookConfiguration<F, T>[];
};

// static steps config

export type PipelineStepConfiguration<F, T, S> = {
  id: ItemId;
  nqName: NqName;
  friendlyName?: string;
  io?: S;
  actions?: StepActionConfiguraion<F, T>[];
};

export type PipelineConfigurationBase<F, T> = {
  id: ItemId;
  nqName?: NqName;
  provider?: NqName;
  version?: string;
  friendlyName?: string;
  hooks?: PipelineHooks<F, T>;
  actions?: PipelineActionConfiguraion<F, T>[];
  states?: StateItem[];
};

// fixed pipeline

export type PipelineStaticItem<F, T, S, R> = {
  id: ItemId;
} & (PipelineStepConfiguration<F, T, S> | AbstractPipelineConfiguration<F, T, S, R> | R);

export type AbstractPipelineStaticConfiguration<F, T, S, R> = {
  links?: PipelineLinkConfiguration<F, T>[];
  steps: PipelineStaticItem<F, T, S, R>[];
  type: 'static';
} & PipelineConfigurationBase<F, T>;

// parallel pipeline

export type ParallelItemContext<F, T> = {
  disableUIAdding?: boolean;
  selectorPath?: F;
  selectorExtractor?: SelectorKeyExtractor;
};

export type PipelineParallelItem<F, T, S, R> = ((PipelineStepConfiguration<F, T, S> | AbstractPipelineConfiguration<F, T, S, R> | R) & ParallelItemContext<F, T>);

export type AbstractPipelineParallelConfiguration<F, T, S, R> = {
  initialSteps?: StepParallelInitialConfig[];
  stepTypes: PipelineParallelItem<F, T, S, R>[];
  type: 'parallel';
} & PipelineConfigurationBase<F, T>;

// sequential pipeline

export type SequentialItemContext = {
  disableUIAdding?: boolean;
};

export type PipelineSequentialItem<F, T, S, R> = ((PipelineStepConfiguration<F, T, S> | AbstractPipelineConfiguration<F, T, S, R> | R) & SequentialItemContext);


// SequentialItemContext;

export type AbstractPipelineSequentialConfiguration<F, T, S, R> = {
  initialSteps?: StepSequentialInitialConfig[];
  stepTypes: PipelineSequentialItem<F, T, S, R>[];
  links?: PipelineDynamicLinkConfiguration<F, T>[];
  type: 'sequential';
} & PipelineConfigurationBase<F, T>;

// pipeline config

export type AbstractPipelineConfiguration<F, T, S, R> =
AbstractPipelineStaticConfiguration<F, T, S, R> |
AbstractPipelineParallelConfiguration<F, T, S, R> |
AbstractPipelineSequentialConfiguration<F, T, S, R>;

export type PipelineRefInitial = {
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<FromSpec, ToSpec, never, PipelineRefInitial>;
export type PipelineConfigurationParallelInitial = AbstractPipelineParallelConfiguration<FromSpec, ToSpec, never, PipelineRefInitial>;
export type PipelineConfigurationSequentialInitial = AbstractPipelineSequentialConfiguration<FromSpec, ToSpec, never, PipelineRefInitial>;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationParallelInitial | PipelineConfigurationSequentialInitial | PipelineRefInitial;

export type PipelineConfiguration = PipelineConfigurationInitial;
