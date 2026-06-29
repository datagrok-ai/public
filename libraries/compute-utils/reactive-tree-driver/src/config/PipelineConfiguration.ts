import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable} from 'rxjs';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimePipelineMutationController, INameSelectorController, IRuntimeValidatorController, IFuncallActionController, IRuntimeReturnController, IRuntimePipelineValidatorController} from '../RuntimeControllers';
import {DynamicPipelineType, ItemId, NqName, RestrictionType, LinkSpecString, ValidationResult} from '../data/common-types';
import {PipelineState, StepDynamicInitialConfig} from './PipelineInstance';
import {LinkIOParsed} from './LinkSpec';
import type ExcelJS from 'exceljs';
import {ConsistencyInfo} from '../runtime/StateTreeNodes';
import {Zippable} from 'fflate';
import {BehaviorSubject} from 'rxjs';

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
  nqName: string,
  version: string | undefined,
  id: string,
}

// handlers

export type LoadedPipeline = (PipelineConfigurationStaticInitial | PipelineConfigurationDynamicInitial) & LoadedPipelineToplevelNode;

export type IRuntimeController = IRuntimeLinkController | IRuntimeValidatorController;
export type HandlerBase<P, R> = ((params: P) => Promise<R> | Observable<R> | R) | NqName;
export type Handler = HandlerBase<{ controller: IRuntimeLinkController }, void>;
export type Validator = HandlerBase<{ controller: IRuntimeValidatorController }, void>;
export type PipelineValidator = HandlerBase<{ controller: IRuntimePipelineValidatorController }, void>;
export type MetaHandler = HandlerBase<{ controller: IRuntimeMetaController }, void>;
export type MutationHandler = HandlerBase<{ controller: IRuntimePipelineMutationController }, void>;
export type SelectorHandler = HandlerBase<{ controller: INameSelectorController }, void>;
export type FunccallActionHandler = HandlerBase<{ controller: IFuncallActionController }, void>;
export type PipelineProvider = HandlerBase<{ version?: string }, LoadedPipeline>;
export type ReturnHandler = HandlerBase<{ controller: IRuntimeReturnController }, void>;

export interface ExportCbInput {
  fc: DG.FuncCall,
  wb: ExcelJS.Workbook,
  archive: Zippable,
  path: string[],
  fileName: string,
  isOutputOutdated?: boolean,
  runError?: string,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
  meta: Record<string, BehaviorSubject<any>>,
  description?: Record<string, string | string[]>,
}

export type ExportUtils = {
  reportStateExcel: (pipelineState: PipelineState, cb?: <T>(input: ExportCbInput) => Promise<T>) => Promise<readonly [Blob, Zippable, string]>;
  reportFuncCallExcel: (fc: DG.FuncCall, uuid: string) => Promise<readonly [Blob, ExcelJS.Workbook]>;
  getFuncCallCustomExports: (fc: DG.FuncCall) => string[];
  runFuncCallCustomExport: (fc: DG.FuncCall, uuid: string, exportName: string) => Promise<any>;

}
export type PipelineExport = (pipelineState: PipelineState, utils: ExportUtils) => Promise<any>;
export type ViewersHook = (ioName: string, type: string, viewer?: DG.Viewer, meta?: any) => void;

// link-like

export type PipelineLinkConfigurationBase<P> = {
  id: ItemId;
  from: P;
  to: P;
  not?: P;
  base?: P,
  actions?: P;
  dataFrameMutations?: boolean | string[];
  defaultRestrictions?: Record<string, RestrictionType> | RestrictionType;
  nodePriority?: number;
}

export type PipelineHandlerConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type?: 'data',
  actions?: undefined;
  handler?: Handler;
  runOnInit?: boolean;
};

export type PipelineValidatorConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'validator'
  handler: Validator;
  runOnInit?: undefined;
  sequential?: boolean;
  debounce?: number;
};

export type PipelineMetaConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'meta'
  actions?: undefined;
  handler: MetaHandler;
  runOnInit?: undefined;
  sequential?: boolean;
};

export type PipelineInitConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type?: 'data',
  base?: undefined,
  actions?: undefined;
  handler: Handler;
  runOnInit?: undefined;
};

export type PipelineReturnConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'return',
  base?: undefined,
  actions?: undefined;
  handler: ReturnHandler;
  runOnInit?: undefined;
};

export type PipelineSelectorConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'nodemeta' | 'selector', // selector for API compat
  actions?: undefined;
  handler: SelectorHandler;
  runOnInit?: undefined;
  sequential?: boolean;
};

export type PipelinePipelineValidatorConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'pipelineValidator';
  handler: PipelineValidator;
  runOnInit?: undefined;
  sequential?: boolean;
  debounce?: number;
};

export type PipelineLinkConfiguration<P> = PipelineHandlerConfiguration<P> | PipelineValidatorConfiguration<P> | PipelineMetaConfiguration<P> | PipelineInitConfiguration<P> | PipelineReturnConfiguration<P> | PipelineSelectorConfiguration<P> | PipelinePipelineValidatorConfiguration<P>;

/** Action fields shared between config-time (ActionInfo<P>) and the UI-facing ViewAction.
 *  Excludes runtime-only matcher fields (showWhen/hideWhen) and UI-only fields (uuid/visible). */
export type ActionInfoBase = {
  id: string;
  position: ActionPositions;
  friendlyName?: string;
  description?: string;
  menuCategory?: string;
  confirmationMessage?: string;
  icon?: string;
  /** Show this action on a specific child step instead of where it's defined.
   *  Uses configId of the target step. The action's from/to still resolve at definition scope. */
  visibleOn?: string;
};

export type ActionInfo<P> = ActionInfoBase & {
  runOnInit?: undefined;
  /** LQL string or array. Action's ViewAction.visible is true only when all
   *  required entries match. (optional) entries are ignored. Empty / absent => no positive gate. */
  showWhen?: P;
  /** LQL string or array. Action's ViewAction.visible is false when any entry matches.
   *  Mirrors the link-level `not` field. Empty / absent => no negative gate. */
  hideWhen?: P;
};

export type DataActionConfiguraion<P> = PipelineLinkConfigurationBase<P> & {
  type?: 'data',
  handler: Handler;
} & ActionInfo<P>;

export type PipelineMutationConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'pipeline',
  handler: MutationHandler;
} & ActionInfo<P>;

export type FuncCallActionConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'funccall',
  handler: FunccallActionHandler;
} & ActionInfo<P>;

const actionPositions = ['buttons', 'menu', 'globalmenu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

type LinkOf<S> = [S] extends [never] ? LinkSpecString : LinkIOParsed[];
type RefOf<S> = [S] extends [never] ? PipelineRefInitial : PipelineSelfRef;
type StatesOf<S> = [S] extends [never] ? Array<ItemId | StateItem> : StateItem[];

// static steps config
export type PipelineStepConfiguration<S> = {
  id: ItemId;
  type?: 'step',
  nqName: NqName;
  friendlyName?: string;
  links?: PipelineLinkConfiguration<LinkOf<S>>[];
  actions?: (DataActionConfiguraion<LinkOf<S>> | FuncCallActionConfiguration<LinkOf<S>>)[];
  states?: StatesOf<S>;
  tags?: string[];
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
  viewersHook?: ViewersHook;
  // Per-step opt-in: when true, TreeWizard shows save-to-history and a history panel for this step.
  enableHistory?: boolean;
  io?: S;
};

export interface CustomExport {
  id: string,
  friendlyName?: string,
  handler: PipelineExport,
}

export type PipelineConfigurationBase<S> = {
  id: ItemId;
  nqName?: NqName;
  version?: string;
  friendlyName?: string;
  description?: string;
  links?: PipelineLinkConfiguration<LinkOf<S>>[];
  actions?: (DataActionConfiguraion<LinkOf<S>> | PipelineMutationConfiguration<LinkOf<S>> | FuncCallActionConfiguration<LinkOf<S>>)[];
  onInit?: PipelineInitConfiguration<LinkOf<S>>;
  onReturn?: PipelineReturnConfiguration<LinkOf<S>>;
  states?: StatesOf<S>;
  tags?: string[];
  forceNavigate?: boolean;
  customExports?: CustomExport[];
  disableHistory?: boolean;
  disableDefaultExport?: boolean;
  approversGroup?: string; // not used rn
};

export type NestedItemContext = {
  disableUIControlls?: boolean;
  disableUIAdding?: boolean;
  disableUIRemoving?: boolean;
  disableUIDragging?: boolean;
};

// action step (lightweight placeholder for displaying actions via visibleOn)

export type AbstractPipelineActionConfiguration = {
  type: 'action';
  id: ItemId;
  nqName?: NqName;
  friendlyName?: string;
  description?: string;
  tags?: string[];
};

// fixed pipeline

export type PipelineStaticItem<S> =
PipelineStepConfiguration<S> | AbstractPipelineConfiguration<S> | AbstractPipelineActionConfiguration | RefOf<S>;

export type AbstractPipelineStaticConfiguration<S> = {
  steps: PipelineStaticItem<S>[];
  type: 'static';
  isActionStep?: boolean;
} & PipelineConfigurationBase<S>;

// dynamic pipeline (unified type for parallel and sequential)

export type PipelineDynamicItem<S> = ((PipelineStepConfiguration<S> | AbstractPipelineConfiguration<S> | AbstractPipelineActionConfiguration | RefOf<S>) & NestedItemContext);

export type AbstractPipelineDynamicConfiguration<S> = {
  initialSteps?: Array<ItemId | StepDynamicInitialConfig>;
  stepTypes: PipelineDynamicItem<S>[];
  type: DynamicPipelineType;
} & PipelineConfigurationBase<S>;

// pipeline config

export type AbstractPipelineConfiguration<S> =
AbstractPipelineStaticConfiguration<S> |
AbstractPipelineDynamicConfiguration<S>;

export type PipelineRefInitial = {
  id?: string;
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<never>;
export type PipelineConfigurationDynamicInitial = AbstractPipelineDynamicConfiguration<never>;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationDynamicInitial | PipelineRefInitial;

export type PipelineConfiguration = PipelineConfigurationInitial;
