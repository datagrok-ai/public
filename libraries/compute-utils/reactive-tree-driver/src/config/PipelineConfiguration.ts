import {Observable} from 'rxjs';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimePipelineMutationController, IRuntimeValidatorController} from '../RuntimeControllers';
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

export type IRuntimeController = IRuntimeLinkController | IRuntimeValidatorController;
export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R> | R) | NqName;
export type Handler = HandlerBase<{ controller: IRuntimeLinkController }, void>;
export type Validator = HandlerBase<{ controller: IRuntimeValidatorController }, void>;
export type MetaHandler = HandlerBase<{ controller: IRuntimeMetaController }, void>;
export type MutationHandler = HandlerBase<{ controller: IRuntimePipelineMutationController }, void>;

export type PipelineProvider = HandlerBase<{ version?: string }, LoadedPipeline>;

// link-like

export type PipelineLinkConfigurationBase<P> = {
  id: ItemId;
  from: P;
  to: P;
  base?: P,
  actions?: P;
  dataFrameMutations?: boolean | string[];
  defaultRestrictions?: Record<string, RestrictionType>;
}

export type PipelineHandlerConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  isValidator?: false;
  isMeta?: false;
  isPipeline?: false;
  actions?: undefined;
  handler?: Handler;
};

export type PipelineValidatorConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  isValidator: true;
  isMeta?: false;
  isPipeline?: false;
  handler: Validator;
};

export type PipelineMetaConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  isValidator?: false;
  isMeta: true;
  isPipeline?: false;
  actions?: undefined;
  handler: MetaHandler;
};

export type PipelineHookConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  isValidator?: false;
  isMeta?: false;
  isPipeline?: false;
  base?: undefined,
  actions?: undefined;
  handler: Handler;
};

export type PipelineLinkConfiguration<P> = PipelineHandlerConfiguration<P> | PipelineValidatorConfiguration<P> | PipelineMetaConfiguration<P> | PipelineHookConfiguration<P>;

export type PipelineActionConfiguraion<P> = PipelineLinkConfigurationBase<P> & {
  position: ActionPositions;
  friendlyName?: string;
  menuCategory?: string;
  handler: Handler;
  isValidator?: false;
  isMeta?: false;
  isPipeline?: false;
};

export type PipelineMutationConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  position: ActionPositions;
  friendlyName?: string;
  menuCategory?: string;
  handler: MutationHandler;
  isValidator?: false;
  isMeta?: false;
  isPipeline: true;
};

export type StepActionConfiguraion<P> = PipelineActionConfiguraion<P>;

const actionPositions = ['buttons', 'menu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

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
  onInit?: PipelineHookConfiguration<P>;
  actions?: (PipelineActionConfiguraion<P> | PipelineMutationConfiguration<P>)[];
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
