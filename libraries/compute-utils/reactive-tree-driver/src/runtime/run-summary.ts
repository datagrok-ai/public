import * as DG from 'datagrok-api/dg';
import {ExportSummaryRollup, RunSummaryItem, StepStatus} from '../config/PipelineConfiguration';
import {PipelineStateRec, StepFunCallStateBase, isFuncCallState} from '../config/PipelineInstance';
import {ValidationItem, ValidationResult} from '../data/common-types';
import {ConsistencyInfo, FuncCallStateInfo} from './StateTreeNodes';

/** Driver states by node uuid */
export type RunStates = {
  callInfoStates?: Record<string, FuncCallStateInfo | undefined>,
  validationStates?: Record<string, Record<string, ValidationResult> | undefined>,
  consistencyStates?: Record<string, Record<string, ConsistencyInfo> | undefined>,
  pipelineValidations?: Record<string, ValidationResult | undefined>,
  descriptions?: Record<string, Record<string, string | string[]> | undefined>,
};

type SummaryNode = PipelineStateRec<StepFunCallStateBase & {funcCall?: DG.FuncCall}, {}>;

// per-input conditions behind a step status
export const inputHasWarning = (val: ValidationResult) => !!val.warnings?.length;
export const inputHasError = (val: ValidationResult) => !!val.errors?.length;
export const isFixableInconsistency = (val: ConsistencyInfo) =>
  !!val.inconsistent && (val.restriction === 'disabled' || val.restriction === 'restricted');
export const isChange = (val: ConsistencyInfo) => !!val.inconsistent && val.restriction === 'info';
export const isInconsistent = (val: ConsistencyInfo) => !!val.inconsistent;

export const anyInput = <T>(states: Record<string, T> | undefined, cond: (state: T) => boolean) =>
  Object.values(states ?? {}).some(cond);

export const inputsWhere = <T>(states: Record<string, T> | undefined, cond: (state: T) => boolean) =>
  Object.entries(states ?? {}).filter(([, state]) => cond(state)).map(([name]) => name);

export const statesToStatus = (
  callState: FuncCallStateInfo,
  validationsState?: Record<string, ValidationResult>,
  consistencyStates?: Record<string, ConsistencyInfo>,
): StepStatus => {
  const hasErrors = anyInput(validationsState, inputHasError);
  const hasWarnings = anyInput(validationsState, inputHasWarning);
  const hasInconsistencies = anyInput(consistencyStates, isFixableInconsistency);
  if (callState.isRunning) return 'running';
  if (callState.runError)
    return 'failed';
  if (callState.pendingDependencies?.length)
    return callState.isOutputOutdated ? 'pending' : 'pending executed';
  if (!callState.isOutputOutdated) {
    if (hasInconsistencies)
      return 'succeeded inconsistent';
    if (hasWarnings || hasErrors)
      return 'succeeded warn';
    if (anyInput(consistencyStates, isChange))
      return 'succeeded info';
    return 'succeeded';
  }
  if (hasErrors)
    return 'next error';
  if (hasWarnings || hasInconsistencies)
    return 'next warn';

  return 'next';
};

export const friendlyIoName = (funcCall: DG.FuncCall | undefined, ioName: string): string => {
  const prop =
    funcCall?.func?.inputs?.find((p: DG.Property) => p.name === ioName) ??
    funcCall?.func?.outputs?.find((p: DG.Property) => p.name === ioName);
  return prop?.friendlyName ?? prop?.caption ?? ioName;
};

// as in the tree, a title given as a list is not used
export const textTitle = (description?: Record<string, string | string[]>) =>
  typeof description?.title === 'string' ? description.title : undefined;

const resultMessages = (result: ValidationResult | undefined, prefix = '') => {
  const texts = (items?: ValidationItem[]) =>
    (items ?? []).map((item) => prefix + (typeof item === 'string' ? item : item.description));
  return {
    errors: texts(result?.errors),
    warnings: texts(result?.warnings),
    notifications: texts(result?.notifications),
  };
};

const stepSummary = (state: SummaryNode & {type: 'funccall'}, states: RunStates): RunSummaryItem => {
  const base = {kind: 'step' as const, uuid: state.uuid, name: state.friendlyName ?? state.configId,
    title: textTitle(states.descriptions?.[state.uuid])};
  const callInfo = states.callInfoStates?.[state.uuid];
  if (!callInfo)
    return {...base, errors: [], warnings: [], notifications: [], inconsistentInputs: []};
  const validation = states.validationStates?.[state.uuid];
  const consistency = states.consistencyStates?.[state.uuid];
  const ios = Object.entries(validation ?? {})
    .map(([io, result]) => resultMessages(result, `${friendlyIoName(state.funcCall, io)}: `));
  return {
    ...base,
    status: statesToStatus(callInfo, validation, consistency),
    runError: callInfo.runError,
    errors: ios.flatMap((m) => m.errors),
    warnings: ios.flatMap((m) => m.warnings),
    notifications: ios.flatMap((m) => m.notifications),
    inconsistentInputs: inputsWhere(consistency, isInconsistent).map((io) => friendlyIoName(state.funcCall, io)),
  };
};

const OUTDATED_STATUSES: StepStatus[] = ['next', 'next warn', 'next error', 'pending'];

const summaryRollup = (steps: RunSummaryItem[]): ExportSummaryRollup => {
  const count = (pred: (item: RunSummaryItem) => boolean) => steps.filter(pred).length;
  return {
    steps: steps.length,
    notLoaded: count((item) => !item.status),
    failed: count((item) => item.status === 'failed'),
    outdated: count((item) => !!item.status && OUTDATED_STATUSES.includes(item.status)),
    withErrors: count((item) => item.errors.length > 0),
    withWarnings: count((item) => item.warnings.length > 0),
    inconsistent: count((item) => item.inconsistentInputs.length > 0),
  };
};

/** Summary items in tree order, a workflow before its steps, with counts over all steps below each workflow. */
export function getRunSummary(treeState: SummaryNode, states: RunStates): RunSummaryItem[] {
  const items: RunSummaryItem[] = [];
  const visit = (state: SummaryNode): RunSummaryItem[] => {
    if (isFuncCallState(state)) {
      const item = stepSummary(state, states);
      items.push(item);
      return [item];
    }
    const workflow: RunSummaryItem = {
      kind: 'workflow',
      uuid: state.uuid,
      name: state.friendlyName ?? state.configId,
      title: textTitle(states.descriptions?.[state.uuid]),
      ...resultMessages(states.pipelineValidations?.[state.uuid]),
      inconsistentInputs: [],
    };
    items.push(workflow);
    const steps = state.steps.flatMap((step) => visit(step));
    workflow.rollup = summaryRollup(steps);
    return steps;
  };
  visit(treeState);
  return items;
}
