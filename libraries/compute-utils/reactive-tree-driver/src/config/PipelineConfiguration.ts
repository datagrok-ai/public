import {Observable} from 'rxjs';
import {RuntimeController} from '../RuntimeController';
import {ItemId, NqName, ItemPath, InputState, ItemType} from '../data/common-types';

//
// Pipeline public configuration
//

export type StateItemConfiguration = {
  id: ItemId;
  type?: string;
};

// handlers

export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R>) | NqName;
export type Handler = HandlerBase<{ controller: RuntimeController }, void>;
export type SelectorKeyExtractor = HandlerBase<{ controller: RuntimeController }, string>;
export type PipelineProvider = HandlerBase<{ version?: string }, PipelineConfigurationProvided>;

// link-like

export type PipelineLinkConfigurationBase<P> = {
  from: P;
  to: P;
  dataFrameMutations?: boolean | P;
  inputState?: InputState | InputState[];
  handler?: Handler;
}

export type PipelineLinkConfiguration<P> = {
  id: ItemId;
} & PipelineLinkConfigurationBase<P>;

export type PipelineDynamicLinkConfiguration<P> = {
  type: ItemType;
  inputType: ItemType;
  outputType: ItemType;
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

// static steps config

export type PipelineStepConfiguration<S> = {
  id: ItemId;
  nqName: NqName;
  friendlyName?: string;
  io?: S;
};

export type PipelineConfigurationBase<P> = {
  id: ItemId;
  nqName: NqName;
  dynamic?: 'sequential' | 'parallel';
  hooks?: PipelineHooks<P>;
  actions?: PipelineActionConfiguraion<P>[];
  states?: StateItemConfiguration[];
};

// fixed pipeline

export type AbstractPipelineStaticConfiguration<P, S, R> = {
  links?: PipelineLinkConfiguration<P>[];
  steps: (PipelineStepConfiguration<S> | AbstractPipelineConfiguration<P, S, R> | R)[];
  dynamic?: undefined;
} & PipelineConfigurationBase<P>;

// parallel pipeline

export type ParallelItemContext<P> = {
  nqName?: NqName;
  type: ItemType;
  selectorPath?: P;
  selectorExtractor?: SelectorKeyExtractor;
};

export type PipelineParallelItem<P, S, R> = ({
  config: AbstractPipelineConfiguration<P, S, R> | R;
}) & ParallelItemContext<P>;

export type AbstractPipelineParallelConfiguration<P, S, R> = {
  items: PipelineParallelItem<P, S, R>[];
  dynamic: 'parallel';
} & PipelineConfigurationBase<P>;

// sequential pipeline

export type SequentialItemContext = {
  nqName?: NqName;
  type: ItemType;
  inputTypeId: ItemId;
  outputTypeId: ItemId;
};

export type PipelineSequentialItem<P, S, R> = ({
  config: AbstractPipelineConfiguration<P, S, R> | R;
}) & SequentialItemContext;

export type AbstractPipelineSequentialConfiguration<P, S, R> = {
  items: PipelineSequentialItem<P, S, R>[];
  links?: PipelineDynamicLinkConfiguration<P>[];
  dynamic: 'sequential';
} & PipelineConfigurationBase<P>;

// global config

export type PipelineRef = {
  id: ItemId;
  version?: string;
  provider: PipelineProvider | NqName;
}

export type AbstractPipelineConfiguration<P, S, R> =
AbstractPipelineStaticConfiguration<P, S, R> |
AbstractPipelineParallelConfiguration<P, S, R> |
AbstractPipelineSequentialConfiguration<P, S, R>;

export type PipelineConfigurationProvided = AbstractPipelineConfiguration<ItemPath | ItemPath[], never, PipelineRef>;
export type PipelineConfiguration = PipelineConfigurationProvided | PipelineRef;
