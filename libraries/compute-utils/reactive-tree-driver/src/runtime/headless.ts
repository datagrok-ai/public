import * as DG from 'datagrok-api/dg';
import {Observable, defer, of} from 'rxjs';
import {concatMap, finalize, map, mapTo} from 'rxjs/operators';
import {TreeNode} from '../data/BaseTree';
import {PipelineConfiguration, RunSummaryItem, StepStatus} from '../config/PipelineConfiguration';
import {PipelineConfigurationProcessed, getProcessedConfig} from '../config/config-processing-utils';
import {PipelineInstanceConfig} from '../config/PipelineInstance';
import {ItemId, NqName, ValidationResult} from '../data/common-types';
import {DriverLogger} from '../data/Logger';
import {ItemMetadata} from '../view/ViewCommunication';
import {callHandler, copyDataFrame, mergeValidationResults} from '../utils';
import {MockRun, StateTree} from './StateTree';
import {ConsistencyInfo, StateTreeNode, isFuncCallNode} from './StateTreeNodes';
import {RunStates, getRunSummary, statesToStatus} from './run-summary';
import {toInstanceConfigRec} from './StateTreeSerializer';

type NodeRunState = {
  id: ItemId;
  uuid: string;
  /** Inputs and, for a step that has run, outputs: the instance config meaning */
  initialValues: Record<string, any>;
  isReadonly?: boolean;
  sourceRunId?: string;
};

export type StepRunState = NodeRunState & {
  status?: StepStatus;
  isOutputOutdated: boolean;
  runError?: string;
  validations: Record<string, ValidationResult>;
  consistency: Record<string, ConsistencyInfo>;
};

export type WorkflowRunState = NodeRunState & {
  skipOnInit?: boolean;
  validations?: ValidationResult;
  steps: RunState[];
};

/** A workflow or step with its values and run status; a valid instance config, so it can start a workflow again. */
export type RunState = StepRunState | WorkflowRunState;

export type HeadlessResult = {
  /** false when a step ready to run had validation errors, then nothing ran */
  isRun: boolean;
  state: RunState;
  summary: RunSummaryItem[];
};

const valuesOf = <T>(subjects: Record<string, {value: T}>) =>
  Object.fromEntries(Object.entries(subjects).map(([uuid, subj]) => [uuid, subj.value]));

function readStates(tree: StateTree): RunStates {
  const pipelineValidations: Record<string, ValidationResult | undefined> = {};
  const descriptions: Record<string, Record<string, string>> = {};
  tree.nodeTree.traverse(tree.nodeTree.root, (acc, node) => {
    const item = node.getItem();
    const title = item.nodeDescription.getState('title');
    if (title != null)
      descriptions[item.uuid] = {title};
    if (!isFuncCallNode(item)) {
      const results = Object.values(item.pipelineValidations$.value);
      pipelineValidations[item.uuid] = results.length ? mergeValidationResults(...results) : undefined;
    }
    return acc;
  }, null);
  return {
    callInfoStates: valuesOf(tree.getFuncCallStates()),
    validationStates: valuesOf(tree.getValidations()),
    consistencyStates: valuesOf(tree.getConsistency()),
    pipelineValidations,
    descriptions,
  };
}

function toRunState(node: TreeNode<StateTreeNode>, conf: PipelineInstanceConfig, states: RunStates): RunState {
  const item = node.getItem();
  const base = {...conf, id: item.config.id, uuid: item.uuid, initialValues: conf.initialValues ?? {}};
  if (isFuncCallNode(item)) {
    const callInfo = states.callInfoStates?.[item.uuid];
    const validations = states.validationStates?.[item.uuid] ?? {};
    const consistency = states.consistencyStates?.[item.uuid] ?? {};
    return {
      ...base,
      status: callInfo && statesToStatus(callInfo, validations, consistency),
      isOutputOutdated: callInfo?.isOutputOutdated ?? true,
      runError: callInfo?.runError,
      validations,
      consistency,
    };
  }
  const steps = node.getChildren().map((child, idx) => toRunState(child.item, conf.steps![idx], states));
  return {...base, validations: states.pipelineValidations?.[item.uuid], steps};
}

// a step may edit a table in place, which must not reach the caller's table
function copyTables<C extends {initialValues?: Record<string, any>, steps?: C[]}>(conf: C): C {
  const initialValues = conf.initialValues && Object.fromEntries(Object.entries(conf.initialValues)
    .map(([name, value]) => [name, value instanceof DG.DataFrame ? copyDataFrame(value) : value]));
  return {...conf, initialValues, steps: conf.steps?.map(copyTables)};
}

// steps that wait for earlier ones get their inputs during the run
function hasReadyStepErrors(tree: StateTree) {
  return tree.nodeTree.traverse(tree.nodeTree.root, (acc, node) => {
    const item = node.getItem();
    return acc || (isFuncCallNode(item) && !item.isReadonly && !item.pendingDependencies$.value.length &&
      !item.instancesWrapper.isRunable$.value);
  }, false);
}

const treeOptions = (config: PipelineConfigurationProcessed, mockMode = false, logger?: DriverLogger) =>
  ({config, isReadonly: false, annotationChecks: true, batchLinks: true, mockMode, logger});

/** Creates the workflow from the instance config, then, unless a step ready to run has validation errors,
 * runs every step that can run, as the Run ready steps button does. */
export function runWorkflowConfig({config, instanceConfig, mockMode, mockRun, logger}: {
  config: PipelineConfigurationProcessed,
  instanceConfig?: PipelineInstanceConfig,
  mockMode?: boolean,
  mockRun?: MockRun,
  logger?: DriverLogger,
}): Observable<HeadlessResult> {
  return defer(() => {
    const options = treeOptions(config, mockMode, logger);
    const tree = instanceConfig ?
      StateTree.fromInstanceConfig({...options, instanceConfig: copyTables(instanceConfig)}) :
      StateTree.fromPipelineConfig(options);
    return tree.init().pipe(
      concatMap(() => {
        if (hasReadyStepErrors(tree))
          return of(false);
        const rootUuid = tree.nodeTree.root.getItem().uuid;
        return tree.runSequence(rootUuid, false, true, false, mockRun).pipe(mapTo(true));
      }),
      map((isRun) => {
        const states = readStates(tree);
        const root = tree.nodeTree.root;
        return {
          isRun,
          state: toRunState(root, toInstanceConfigRec(root), states),
          summary: getRunSummary(tree.toState(), states),
        };
      }),
      finalize(() => tree.close()),
    );
  });
}

/** Saves a run state as a new workflow run, the way Duplicate builds a copy; a failed step is saved as not run. */
export function saveWorkflowConfig({config, state, metadata}: {
  config: PipelineConfigurationProcessed,
  state: RunState,
  metadata?: ItemMetadata,
}): Observable<DG.FuncCall> {
  return defer(() => {
    const tree = StateTree.fromInstanceConfig({...treeOptions(config), instanceConfig: state});
    return tree.init().pipe(
      concatMap(() => tree.save(undefined, metadata)),
      map((call) => call!),
      finalize(() => tree.close()),
    );
  });
}

const loadConfig = (provider: NqName, version?: string) =>
  callHandler<PipelineConfiguration>(provider, {version}).pipe(concatMap((conf) => getProcessedConfig(conf)));

export function runWorkflow(provider: NqName, version?: string, instanceConfig?: PipelineInstanceConfig) {
  return loadConfig(provider, version).pipe(
    concatMap((config) => runWorkflowConfig({config, instanceConfig})),
  ).toPromise() as Promise<HeadlessResult>;
}

export function saveWorkflowRun(
  provider: NqName, version: string | undefined, state: RunState, metadata?: ItemMetadata,
) {
  return loadConfig(provider, version).pipe(
    concatMap((config) => saveWorkflowConfig({config, state, metadata})),
  ).toPromise() as Promise<DG.FuncCall>;
}
