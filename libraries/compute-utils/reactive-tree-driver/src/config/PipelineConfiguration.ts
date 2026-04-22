import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable} from 'rxjs';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimePipelineMutationController, INameSelectorController, IRuntimeValidatorController, IFuncallActionController, IRuntimeReturnController} from '../RuntimeControllers';
import {ItemId, NqName, RestrictionType, LinkSpecString, ValidationResult} from '../data/common-types';
import {PipelineOutline, PipelineState, StepDynamicInitialConfig, StepParallelInitialConfig, StepSequentialInitialConfig} from './PipelineInstance';
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
export type StructureCheckHook = (data: PipelineOutline) => ValidationResult | undefined;

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

export type PipelineLinkConfiguration<P> = PipelineHandlerConfiguration<P> | PipelineValidatorConfiguration<P> | PipelineMetaConfiguration<P> | PipelineInitConfiguration<P> | PipelineReturnConfiguration<P> | PipelineSelectorConfiguration<P>;

export type ActionInfo = {
  id: string;
  position: ActionPositions;
  friendlyName?: string;
  description?: string;
  menuCategory?: string;
  confirmationMessage?: string;
  icon?: string;
  runOnInit?: undefined;
  /** Show this action on a specific child step instead of where it's defined.
   *  Uses configId of the target step. The action's from/to still resolve at definition scope. */
  visibleOn?: string;
};

export type DataActionConfiguraion<P> = PipelineLinkConfigurationBase<P> & {
  type?: 'data',
  handler: Handler;
} & ActionInfo;

export type PipelineMutationConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'pipeline',
  handler: MutationHandler;
} & ActionInfo;

export type FuncCallActionConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'funccall',
  handler: FunccallActionHandler;
} & ActionInfo;

const actionPositions = ['buttons', 'menu', 'globalmenu', 'none'] as const;
export type ActionPositions = typeof actionPositions[number];

// static steps config
export type PipelineStepConfiguration<P, S> = {
  id: ItemId;
  type?: 'step',
  nqName: NqName;
  friendlyName?: string;
  actions?: (DataActionConfiguraion<P> | FuncCallActionConfiguration<P>)[];
  states?: StateItem[];
  tags?: string[];
  initialValues?: Record<string, any>;
  inputRestrictions?: Record<string, RestrictionType>;
  viewersHook?: ViewersHook;
  io?: S;
};

export interface CustomExport {
  id: string,
  friendlyName?: string,
  handler: PipelineExport,
}

export type PipelineConfigurationBase<P> = {
  id: ItemId;
  nqName?: NqName;
  version?: string;
  friendlyName?: string;
  links?: PipelineLinkConfiguration<P>[];
  actions?: (DataActionConfiguraion<P> | PipelineMutationConfiguration<P> | FuncCallActionConfiguration<P>)[];
  onInit?: PipelineInitConfiguration<P>;
  onReturn?: PipelineReturnConfiguration<P>;
  states?: StateItem[];
  tags?: string[];
  structureCheck?: StructureCheckHook;
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
  friendlyName?: string;
  description?: string;
  tags?: string[];
};

// fixed pipeline

export type PipelineStaticItem<P, S, R> =
PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R> | AbstractPipelineActionConfiguration | R;

export type AbstractPipelineStaticConfiguration<P, S, R> = {
  steps: PipelineStaticItem<P, S, R>[];
  type: 'static';
  isActionStep?: boolean;
} & PipelineConfigurationBase<P>;

// dynamic pipeline (unified type for parallel and sequential)

export type PipelineDynamicItem<P, S, R> = ((PipelineStepConfiguration<P, S> | AbstractPipelineConfiguration<P, S, R> | AbstractPipelineActionConfiguration | R) & NestedItemContext);

export type AbstractPipelineDynamicConfiguration<P, S, R> = {
  initialSteps?: StepDynamicInitialConfig[];
  stepTypes: PipelineDynamicItem<P, S, R>[];
  type: 'dynamic' | 'parallel' | 'sequential';
} & PipelineConfigurationBase<P>;

/** @deprecated Use PipelineDynamicItem */
export type PipelineParallelItem<P, S, R> = PipelineDynamicItem<P, S, R>;
/** @deprecated Use PipelineDynamicItem */
export type PipelineSequentialItem<P, S, R> = PipelineDynamicItem<P, S, R>;
/** @deprecated Use AbstractPipelineDynamicConfiguration */
export type AbstractPipelineParallelConfiguration<P, S, R> = AbstractPipelineDynamicConfiguration<P, S, R>;
/** @deprecated Use AbstractPipelineDynamicConfiguration */
export type AbstractPipelineSequentialConfiguration<P, S, R> = AbstractPipelineDynamicConfiguration<P, S, R>;

// pipeline config

export type AbstractPipelineConfiguration<P, S, R> =
AbstractPipelineStaticConfiguration<P, S, R> |
AbstractPipelineDynamicConfiguration<P, S, R>;

export type PipelineRefInitial = {
  id?: string;
  version?: string;
  provider: PipelineProvider | NqName;
  type: 'ref';
}

export type PipelineConfigurationStaticInitial = AbstractPipelineStaticConfiguration<LinkSpecString, never, PipelineRefInitial>;
export type PipelineConfigurationDynamicInitial = AbstractPipelineDynamicConfiguration<LinkSpecString, never, PipelineRefInitial>;
/** @deprecated Use PipelineConfigurationDynamicInitial */
export type PipelineConfigurationParallelInitial = PipelineConfigurationDynamicInitial;
/** @deprecated Use PipelineConfigurationDynamicInitial */
export type PipelineConfigurationSequentialInitial = PipelineConfigurationDynamicInitial;

export type PipelineConfigurationInitial = PipelineConfigurationStaticInitial | PipelineConfigurationDynamicInitial | PipelineRefInitial;

export type PipelineConfiguration = PipelineConfigurationInitial;
