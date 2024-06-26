import {Observable} from 'rxjs';
import {RuntimeController} from '../RuntimeController';
import {ItemName, NqName, ItemPath, InputState, TypeKey} from './CommonTypes';

//
// CompositionPipeline public configuration
//

export type StateItemConfiguration = {
  id: ItemName;
};

export type Handler = ((params: { controller: RuntimeController }) => Promise<void> | Observable<void>) | NqName;
export type SelectorKeyExtractor = ((params: { controller: RuntimeController }) => Promise<string> | Observable<string>) | NqName;

export type PipelineLinkConfigurationBase = {
  includeDataFrameMutations?: boolean | ItemPath[];
  inputState?: InputState | InputState[];
  handler?: Handler;
}

export type PipelineLinkConfiguration = {
  id: ItemName;
  from: ItemPath | ItemPath[];
  to: ItemPath | ItemPath[];
} & PipelineLinkConfigurationBase;

export type PipelineDynamicLinkConfiguration = {
  inputType: TypeKey;
  outputType: TypeKey;
} & PipelineLinkConfigurationBase;

export type PipelineHookConfiguration = {
  id: ItemName;
  from?: ItemPath | ItemPath[];
  to?: ItemPath | ItemPath[];
  handler: Handler;
};

const actionPositions = ['buttons', 'menu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

export type PipelineActionConfiguraion = PipelineHookConfiguration & {
  position: ActionPositions;
  friendlyName?: string;
};

export type PipelineStepConfiguration = {
  id: ItemName;
  nqName: NqName;
  friendlyName?: string;
  actions?: PipelineActionConfiguraion[];
};

export type PipelineHooks = {
  onInit?: PipelineHookConfiguration[];
  beforeLoadFuncCall?: PipelineHookConfiguration[];
  afterLoadFuncCall?: PipelineHookConfiguration[];
  beforeInputFormRender?: PipelineHookConfiguration[];
  afterInputFormRender?: PipelineHookConfiguration[];
  beforeViewerRender?: PipelineHookConfiguration[];
  afterViewerRender?: PipelineHookConfiguration[];
  beforeLoadRun?: PipelineHookConfiguration[];
  afterLoadRun?: PipelineHookConfiguration[];
  beforeSaveRun?: PipelineHookConfiguration[];
  afterSaveRun?: PipelineHookConfiguration[];
  onClose?: PipelineHookConfiguration[];
};

export type PipelineConfigurationBase = {
  id: ItemName;
  nqName: NqName;
  dynamic?: 'sequential' | 'parallel';
  hooks?: PipelineHooks;
  actions?: PipelineActionConfiguraion[];
  states?: StateItemConfiguration[];
};

// fixed pipeline

export type FixedPipelineConfiguration = {
  links?: PipelineLinkConfiguration[];
  steps: (PipelineStepConfiguration | NestedPipelineConfiguration)[];
  dynamic: undefined;
} & PipelineConfigurationBase;

// parallel pipeline

export type PipelineParallelItem = ({
  nqName: NqName;
} | {
  config: PipelineConfiguration;
}) & {
  type: TypeKey;
  typeName: string;
  selectorPath: ItemPath[];
  selectorHandler?: SelectorKeyExtractor;
}

export type PipelineParallelConfiguration = {
  items: PipelineParallelItem[];
  dynamic: 'parallel';
} & PipelineConfigurationBase;

// sequential pipeline

export type PipelineSequentialItem = ({
  nqName: NqName;
} | {
  config: PipelineConfiguration;
}) & {
  type: TypeKey;
  typeName: string;
  inputTypeId: ItemName;
  outputTypeId: ItemName;
};

export type PipelineSequentialHandlers = Record<ItemName, Record<ItemName, Handler>>;

export type PipelineSequentialConfiguration = {
  items: PipelineSequentialItem[];
  handlers: PipelineSequentialHandlers;
  links?: PipelineDynamicLinkConfiguration[];
  dynamic: 'sequential';
} & PipelineConfigurationBase;

// global config

export type PipelineOptions = {
  actions?: [PipelineActionConfiguraion, ItemPath][];
};

export type PipelineRef = {
  provider: NqName;
}

export type PipelineConfiguration = FixedPipelineConfiguration | PipelineParallelConfiguration | PipelineSequentialConfiguration;

export type NestedPipelineConfiguration = (PipelineConfiguration | PipelineRef) & PipelineOptions
