import type {CheckOptions, CheckSeverity} from './checks';
import * as DG from 'datagrok-api/dg';
import {Observable} from 'rxjs';
import type {RulesLogic} from 'json-logic-js';
import {IRuntimeLinkController, IRuntimeMetaController, IRuntimePipelineMutationController, INameSelectorController, IRuntimeValidatorController, IFuncallActionController, IRuntimeReturnController, IRuntimePipelineValidatorController} from '../RuntimeControllers';
import {
  AnnotationLinkKind, DynamicPipelineType, ItemId, NqName, RestrictionType, LinkSpecString, ValidationResult,
} from '../data/common-types';
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

export type StepStatus = 'next' | 'next warn' | 'next error' |
  'pending' | 'pending executed' |
  'running' |
  'succeeded' | 'succeeded info' | 'succeeded warn' | 'succeeded inconsistent' |
  'failed';

export interface ExportSummaryRollup {
  steps: number,
  notLoaded: number,
  failed: number,
  outdated: number,
  withErrors: number,
  withWarnings: number,
  inconsistent: number,
}

export interface ExportSummaryItem {
  kind: 'step' | 'workflow',
  /** Folders in the exported zip */
  path: string[],
  name: string,
  title?: string,
  /** The step workbook in the zip; absent for workflows and for steps that are not loaded */
  fileName?: string,
  /** The status the tree shows; absent for workflows and for steps that are not loaded */
  status?: StepStatus,
  runError?: string,
  errors: string[],
  warnings: string[],
  notifications: string[],
  inconsistentInputs: string[],
  /** Workflows only: counts over all steps below */
  rollup?: ExportSummaryRollup,
}

export interface ExportCbInput {
  fc: DG.FuncCall,
  wb: ExcelJS.Workbook,
  archive: Zippable,
  path: string[],
  fileName: string,
  status?: StepStatus,
  isOutputOutdated?: boolean,
  runError?: string,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
  meta: Record<string, BehaviorSubject<any>>,
  description?: Record<string, string | string[]>,
}

export type ExportUtils = {
  reportStateExcel: (pipelineState: PipelineState, cb?: <T>(input: ExportCbInput) => Promise<T>) => Promise<readonly [Blob, Zippable, string, ExportSummaryItem[]]>;
  getExportSummary: (pipelineState: PipelineState) => ExportSummaryItem[];
  reportSummaryExcel: (pipelineState: PipelineState) => Promise<readonly [Blob, ExcelJS.Workbook]>;
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
  params?: Record<string, any>;
  /** Set by config processing on links generated from function annotations. */
  annotation?: [P] extends [LinkSpecString] ? never : AnnotationLinkKind;
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

// rule links (expanded into meta/validator/data links at config processing)

type RuleOps =
  | {literal: any}
  | {columns: RuleExpr | [RuleExpr] | [RuleExpr, string]}
  | {columnsMissing: [RuleExpr, RuleExpr]}
  | {columnIs: [RuleExpr, string]}
  | {nulls: RuleExpr}
  | {column: [RuleExpr, string]}
  | {row: [RuleExpr, string, RuleExpr]}
  | {regex: [RuleExpr, string] | [RuleExpr, string, string]}
  | {script: string}
  | {scriptVerdict: string}
  | {len: RuleExpr};

/** JSON Logic expression, or formula text such as `'gt(m, 0)'`. */
export type RuleLogic = RulesLogic<RuleOps>;
/** JSON Logic expression or a plain literal; in value fields a string starting with `=` is a formula. */
export type RuleExpr = RuleLogic | RuleExpr[];
export type RuleTargets = string | string[];

/** An effect's own condition, combined with the rule's `when`. */
type RuleEffectWhen = {when?: RuleLogic};

export type RuleMetaEffect = RuleEffectWhen & (
  | {effect: 'hide' | 'show', targets: RuleTargets}
  | {effect: 'items', targets: RuleTargets, items: RuleExpr}
  | {effect: 'meta', targets: RuleTargets, meta: Record<string, RuleExpr>});

export type RuleValidatorEffect = RuleEffectWhen & (
  | {effect: 'error' | 'warning' | 'notification', targets: RuleTargets, message: RuleExpr}
  /** Writes a source's verdicts: `isError` items as errors, the rest as warnings. */
  | {effect: 'verdicts', targets: RuleTargets, source: string});

export type RuleDataEffect = RuleEffectWhen & (
  | {effect: 'set', targets: RuleTargets, value: RuleExpr, restriction?: RestrictionType}
  | {effect: 'clear', targets: RuleTargets, restriction?: RestrictionType}
  /** Writes each key of the `values` object to the target alias of the same name;
   *  keys without a target are ignored, targets without a key are left as they are.
   *  Without `targets` every `to` alias is a target. */
  | {effect: 'assign', targets?: RuleTargets, values: RuleExpr, restriction?: RestrictionType, ignoreCase?: boolean});

export type RuleEffect = RuleMetaEffect | RuleValidatorEffect | RuleDataEffect;

/** Values resolved before a rule evaluates. `validators` runs the named functions
 *  on the io behind `input`; without `names` it runs the io's own annotation validators
 *  through the platform, using the `(call)` input that expansion adds as `call`.
 *  `js` calls `fn` with the values of the `args` input aliases on every run of the
 *  rule; a returned promise is awaited. */
export type RuleSource =
  {validators: {input: string, names?: string[], call?: string}} |
  /** The annotation `choices` of the io behind `input`, evaluated by the platform (1.28+) through the step's
   *  FuncCall: `{items, values, inList, row}`, `row` being the `propagateChoice` lookup row of the current value. */
  {choices: {input: string, call?: string}} |
  {js: {args: string[], fn: (...values: any[]) => any}} |
  /** Calls the platform function `name`; `args` maps its parameters to expressions over the inputs. */
  {func: {name: string, args?: Record<string, RuleExpr>}} |
  /** Runs `sql` on the connection; `args` binds the query's `@name` parameters. */
  {query: {connection: string, sql: string, args?: Record<string, RuleExpr>}} |
  /** Loads a table from a file share path or a URL, once per link. */
  {file: string} |
  /** A table given in the config: a dataframe, or CSV text parsed with optional import options, once per link. */
  {table: DG.DataFrame | string | {csv: string, options?: DG.CsvImportOptions}};

export type PipelineRuleConfiguration<P> = PipelineLinkConfigurationBase<P> & {
  type: 'rule';
  when?: RuleLogic;
  /** Source objects, or source calls such as `'file("System:AppData/Pkg/presets.csv")'`. */
  sources?: Record<string, RuleSource | string>;
  /** Effect objects, or effect calls such as `'set(t, m, restriction: "restricted")'`. */
  effects: (RuleEffect | string)[];
  debounce?: number;
  runOnInit?: boolean;
  handler?: undefined;
  actions?: undefined;
  params?: undefined;
};

/** Annotation-style checks on one io, expanded into validator links at config processing. */
export type PipelineCheckConfiguration<P> = {
  id: ItemId;
  type: 'check';
  /** LQL query of the checked io, without an alias. */
  io: P;
  /** The annotation options to check; `table` is an LQL query of the table io, without an alias. */
  check: CheckOptions;
  /** Inputs a GrokScript expression reads, as variable name to io query; `value` is the checked io. */
  vars?: Record<string, P>;
  when?: RuleLogic;
  message?: RuleExpr;
  severity?: CheckSeverity;
  not?: P;
  base?: P;
  nodePriority?: number;
  debounce?: number;
};

export type PipelineLinkConfigurationInput<P> = PipelineLinkConfiguration<P> | PipelineRuleConfiguration<P> | PipelineCheckConfiguration<P>;

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
type LinksOf<S> = [S] extends [never] ?
  PipelineLinkConfigurationInput<LinkSpecString>[] :
  PipelineLinkConfiguration<LinkIOParsed[]>[];
type RefOf<S> = [S] extends [never] ? PipelineRefInitial : PipelineSelfRef;
type StatesOf<S> = [S] extends [never] ? Array<ItemId | StateItem> : StateItem[];

// static steps config
export type PipelineStepConfiguration<S> = {
  id: ItemId;
  type?: 'step',
  nqName: NqName;
  friendlyName?: string;
  links?: LinksOf<S>;
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
  links?: LinksOf<S>;
  actions?: (DataActionConfiguraion<LinkOf<S>> | PipelineMutationConfiguration<LinkOf<S>> | FuncCallActionConfiguration<LinkOf<S>>)[];
  onInit?: PipelineInitConfiguration<LinkOf<S>>;
  onReturn?: PipelineReturnConfiguration<LinkOf<S>>;
  states?: StatesOf<S>;
  tags?: string[];
  forceNavigate?: boolean;
  customExports?: CustomExport[];
  disableHistory?: boolean;
  disableDefaultExport?: boolean;
  compactView?: boolean;
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
