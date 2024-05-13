import {Observable} from 'rxjs';
import {RichFunctionView} from '../../function-views';
import {ActionItem, ValidationResult} from '../../shared-utils/validation';

//
// State optional spec
//
const primitiveTypes = ['bool', 'int', 'number', 'string'] as const;
type SupportedPrimitives = (typeof primitiveTypes)[number];
const containerTypes = ['dataframe', 'object'] as const;
type SupportedContainers = (typeof primitiveTypes)[number];
const supportedTypes = [...primitiveTypes, ...containerTypes] as const;
type SupportedTypes = (typeof supportedTypes)[number];

//
// Common interfaces
//

export type TypeSpec = {
  [key: string]: SupportedTypes | { type: SupportedContainers; spec?: TypeSpec; };
};

export type ItemName = string;
export type ItemPath = ItemName[];
export type NqName = string;

export type InputState = 'disabled' | 'restricted' | 'user input';

export interface RuntimeController {
  enableLink(path: ItemPath): void;
  disableLink(path: ItemPath): void;
  triggerLink(path: ItemPath): void;
  isLinkEnabled(path: ItemPath): void;
  enableStep(path: ItemPath): void;
  disableStep(path: ItemPath): void;
  isStepEnabled(path: ItemPath): boolean;
  enablePipeline(path: ItemPath): void;
  disablePipeline(path: ItemPath): void;
  loadNestedPipeline(path: ItemPath, runId: string): void;
  getState<T = any>(path: ItemPath): T | void;
  setState<T = any>(path: ItemPath, state: T, inputState?: InputState): void;
  updateStateInputState(path: ItemPath, inputState?: InputState): void;
  getValidationAction(path: ItemPath, name?: string): ActionItem;
  setValidation(path: ItemPath, validation?: ValidationResult | undefined): void;
  getView(path: ItemPath): RichFunctionView | void;
  goToStep(path: ItemPath): void;
  getCloseNotifier(): Observable<boolean>;
}

//
// CompositionPipeline public configuration
//

export type ExportConfig = {
  export?: ((format: string) => Promise<Blob>);
  filename?: ((format: string) => string);
  supportedFormats?: string[];
  supportedExtensions?: Record<string, string>;
};

export type StateType = 'input' | 'output' | 'state';

export type StateItemConfiguration = {
  id: ItemName;
  type?: TypeSpec;
  stateType?: StateType; // system field
};

export type Handler = ((params: { controller: RuntimeController; }) => Promise<void>) | NqName;

export type PipelineLinkConfiguration = {
  id: ItemName;
  from: ItemPath | ItemPath[];
  to: ItemPath | ItemPath[];
  handler?: Handler;
  ignoreNotifier?: boolean;
  includeDataFrameMutations?: boolean;
  inputState?: InputState;
};

export type PipelineHookConfiguration = {
  id: ItemName;
  from?: ItemPath | ItemPath[];
  to?: ItemPath | ItemPath[];
  handler: Handler;
};

export type PipelineActionConfiguraion = PipelineHookConfiguration & {
  position: 'buttons' | 'menu' | 'none';
  friendlyName: string;
};

export type PipelinePopupConfiguration = {
  id: ItemName;
  nqName: NqName;
  position: 'buttons' | 'menu' | 'none';
  friendlyName: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: PipelineActionConfiguraion[];
};

export type PipelineStepConfiguration = {
  id: ItemName;
  nqName: NqName;
  friendlyName?: string;
  helpUrl?: string;
  states?: StateItemConfiguration[];
  actions?: PipelineActionConfiguraion[];
  popups?: PipelinePopupConfiguration[];
};

export type PipelineHooks = {
  beforeInit?: PipelineHookConfiguration[];
  afterInit?: PipelineHookConfiguration[];
  beforeFuncCallReady?: PipelineHookConfiguration[];
  afterFuncCallReady?: PipelineHookConfiguration[];
  beforeLoadRun?: PipelineHookConfiguration[];
  afterLoadRun?: PipelineHookConfiguration[];
  beforeSaveRun?: PipelineHookConfiguration[];
  afterSaveRun?: PipelineHookConfiguration[];
  onViewReady?: PipelineHookConfiguration[];
};

export type PipelineConfiguration = {
  id: ItemName;
  nqName: NqName;
  steps: PipelineStepConfiguration[];
  hooks?: PipelineHooks;
  links?: PipelineLinkConfiguration[];
  actions?: PipelineActionConfiguraion[];
  states?: StateItemConfiguration[];
  exportConfig?: ExportConfig;
  helpUrl?: string;
};

//
// Composition Pipeline configuration
//

export type ItemsToAdd = {
  popupsToAdd?: [PipelinePopupConfiguration, ItemPath][];
  actionsToAdd?: [PipelineActionConfiguraion, ItemPath][];
};

export type ItemsToRemove = {
  itemsToRemove?: ItemPath[];
};

export type NestedPipelineConfig = {
  insertBeforeStep?: ItemName;
  friendlyPrefix?: string;
};

export type ComposedPipelinesConfig = {
  nestedPipelinesConfig?: {
    [id: string]: NestedPipelineConfig;
  };
};

export type PipelineCompositionConfiguration = ItemsToAdd & PipelineConfiguration & ComposedPipelinesConfig & ItemsToRemove;
