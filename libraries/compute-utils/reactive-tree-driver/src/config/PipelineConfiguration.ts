import {Observable} from 'rxjs';
import {IRuntimeLinkController, IRuntimePipelineController, IRuntimeValidatorController} from '../RuntimeControllers';
import {ItemId, NqName, RestrictionType, LinkSpecString} from '../data/common-types';
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

export type IRuntimeController = IRuntimeLinkController | IRuntimePipelineController | IRuntimeValidatorController;
export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R>) | NqName;
export type Handler = HandlerBase<{ controller: IRuntimeController }, void>;
export type PipelineProvider = HandlerBase<{ version?: string }, LoadedPipeline>;

// link-like

export type PipelineLinkConfigurationBase<P> = {
  from: P;
  to: P;
  base?: P,
  dataFrameMutations?: boolean | string[] | Record<string, string>;
  inputState?: RestrictionType | RestrictionType[] | Record<string, RestrictionType>;
  handler?: Handler;
}

export type PipelineLinkConfiguration<P> = {
  id: ItemId;
} & PipelineLinkConfigurationBase<P>;

export type PipelineHookConfiguration<P> = {
  id: ItemId;
  handler: Handler;
} & Partial<PipelineLinkConfigurationBase<P>>;

export type PipelineActionConfiguraion<P> = {
  id: ItemId;
  friendlyName?: string;
  menuCategory?: string;
} & PipelineLinkConfigurationBase<P>;

export type StepActionConfiguraion<P> = {
  id: ItemId;
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

// static steps config

export type PipelineStepConfiguration<P, S> = {
  id: ItemId;
  nqName: NqName;
  friendlyName?: string;
  io?: S;
  actions?: StepActionConfiguraion<P>[];
};

export type PipelineConfigurationBase<P> = {
  id: ItemId;
  nqName?: NqName;
  provider?: NqName;
  version?: string;
  friendlyName?: string;
  hooks?: PipelineHooks<P>;
  actions?: PipelineActionConfiguraion<P>[];
  states?: StateItem[];
};

// fixed pipeline

export type PipelineStaticItem<P, S, R> = {
  id: ItemId;
} & (PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R> | R);

export type AbstractPipelineStaticConfiguration<P, S, R> = {
  links?: PipelineLinkConfiguration<P>[];
  steps: PipelineStaticItem<P, S, R>[];
  type: 'static';
} & PipelineConfigurationBase<P>;

// parallel pipeline

export type ParallelItemContext<P> = {
  disableUIAdding?: boolean;
  selectorPath?: P;
};

export type PipelineParallelItem<P, S, R> = ((PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R> | R) & ParallelItemContext<P>);

export type AbstractPipelineParallelConfiguration<P, S, R> = {
  initialSteps?: StepParallelInitialConfig[];
  stepTypes: PipelineParallelItem<P, S, R>[];
  type: 'parallel';
} & PipelineConfigurationBase<P>;

// sequential pipeline

export type SequentialItemContext = {
  disableUIAdding?: boolean;
};

export type PipelineSequentialItem<P, S, R> = ((PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R> | R) & SequentialItemContext);


// SequentialItemContext;

export type AbstractPipelineSequentialConfiguration<P, S, R> = {
  initialSteps?: StepSequentialInitialConfig[];
  stepTypes: PipelineSequentialItem<P, S, R>[];
  links?: PipelineLinkConfiguration<P>[];
  type: 'sequential';
} & PipelineConfigurationBase<P>;

// pipeline config

export type AbstractPipelineConfiguration<P, S, R> =
AbstractPipelineStaticConfiguration<P, S, R> |
AbstractPipelineParallelConfiguration<P, S, R> |
AbstractPipelineSequentialConfiguration<P, S, R>;

export type PipelineRefInitial = {
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<LinkSpecString, never, PipelineRefInitial>;
export type PipelineConfigurationParallelInitial = AbstractPipelineParallelConfiguration<LinkSpecString, never, PipelineRefInitial>;
export type PipelineConfigurationSequentialInitial = AbstractPipelineSequentialConfiguration<LinkSpecString, never, PipelineRefInitial>;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationParallelInitial | PipelineConfigurationSequentialInitial | PipelineRefInitial;

export type PipelineConfiguration = PipelineConfigurationInitial;
