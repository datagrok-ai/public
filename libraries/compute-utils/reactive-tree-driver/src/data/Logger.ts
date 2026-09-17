import * as grok from 'datagrok-api/grok';
import {v4 as uuidv4} from 'uuid';
import {NodePath} from './BaseTree';
import {BehaviorSubject, Subject} from 'rxjs';

export type DebugLogType =
 'linkRunScheduled' | 'linkRunFinished' | 'linkAdded' | 'linkRemoved' | 'actionAdded' | 'actionRemoved' |
 'treeUpdateStarted' | 'treeUpdateMutation' | 'treeUpdateFinished'

export interface DebugLogBase {
  uuid: string,
  timestamp: Date,
}

export interface LinkLogPayload {
  linkUUID: string,
  prefix: Readonly<NodePath>,
  basePath?: Readonly<NodePath>,
  id: string,
  isDefaultValidator?: boolean,
}

export interface LinkLogItem extends DebugLogBase, LinkLogPayload {}

export interface LinkRunStartedLogItem extends LinkLogItem {
  type: 'linkRunStarted',
}

export interface LinkRunFinishedLogItem extends LinkLogItem {
  type: 'linkRunFinished',
}

export interface LinkAddedLogItem extends LinkLogItem {
  type: 'linkAdded',
}

export interface ActionAddedLogItem extends LinkLogItem {
  type: 'actionAdded',
}

export interface LinkRemoveLogItem extends LinkLogItem {
  type: 'linkRemoved',
}

export interface ActionRemoveLogItem extends LinkLogItem {
  type: 'actionRemoved',
}

export interface TreeUpdateStartedLogItem extends DebugLogBase {
  type: 'treeUpdateStarted',
}

export interface TreeUpdateMutationPayload {
  mutationRootPath?: NodePath,
  addIdx?: number,
  removeIdx?: number,
  id: string,
}

export interface TreeUpdateMutationLogItem extends DebugLogBase, TreeUpdateMutationPayload {
  type: 'treeUpdateMutation',
}

export interface TreeUpdateFinishedLogItem extends DebugLogBase {
  type: 'treeUpdateFinished',
}

export type ErrorSeverity = 'fatal' | 'recoverable' | 'warning';

export interface ErrorLogItem extends DebugLogBase {
  type: 'error';
  severity: ErrorSeverity;
  context: string;
  message: string;
  error?: Error;
  links?: string[];
}

export type LogItem =
TreeUpdateStartedLogItem | TreeUpdateMutationLogItem | TreeUpdateFinishedLogItem |
LinkRunStartedLogItem | LinkRunFinishedLogItem | LinkAddedLogItem | LinkRemoveLogItem | ActionAddedLogItem | ActionRemoveLogItem |
ErrorLogItem;

const MAX_LOG_ENTRIES = 5000;

// emitted arrays are immutable snapshots: append copies into a fresh array,
// so consumers (vue refs, logs$.value readers) never observe mutation
function capAppend<T>(arr: readonly T[], item: T): T[] {
  const next = arr.length >= MAX_LOG_ENTRIES ? arr.slice(arr.length - MAX_LOG_ENTRIES + 1) : arr.slice();
  next.push(item);
  return next;
}

export class DriverLogger {
  private log: LogItem[] = [];
  public logs$ = new BehaviorSubject<LogItem[]>([]);
  private _errors: ErrorLogItem[] = [];
  public errors$ = new Subject<ErrorLogItem>();

  get errors(): readonly ErrorLogItem[] { return this._errors; }

  logLink(type: 'linkRunStarted' | 'linkRunFinished' | 'linkAdded' | 'linkRemoved' | 'actionAdded' | 'actionRemoved', data: LinkLogPayload) {
    const uuid = uuidv4();
    const timestamp = new Date();
    this.append({type, uuid, timestamp, ...data});
  }

  logTreeUpdates(type: 'treeUpdateStarted' | 'treeUpdateFinished') {
    const uuid = uuidv4();
    const timestamp = new Date();
    this.append({type, uuid, timestamp});
  }

  logMutations(data: TreeUpdateMutationPayload) {
    const uuid = uuidv4();
    const timestamp = new Date();
    // no emission: mutation entries are batched until the next event flushes them
    this.log = capAppend(this.log, {type: 'treeUpdateMutation', uuid, timestamp, ...data});
  }

  logError(severity: ErrorSeverity, context: string, error: Error, links?: string[]) {
    const entry: ErrorLogItem = {
      type: 'error',
      uuid: uuidv4(),
      timestamp: new Date(),
      severity,
      context,
      message: error.message,
      error,
      links,
    };
    this._errors = capAppend(this._errors, entry);
    this.errors$.next(entry);
    this.append(entry);
    reportErrorToUI(severity, context, error);
  }

  private append(item: LogItem) {
    this.log = capAppend(this.log, item);
    this.logs$.next(this.log);
  }
}

export function reportErrorToUI(severity: ErrorSeverity, context: string, error: Error | string) {
  const msg = error instanceof Error ? error.message : error;
  console.error(`[RTD:${context}] ${msg}`, error);
  if (severity === 'fatal' || severity === 'recoverable')
    grok.shell.error(msg);
  else
    grok.shell.warning(msg);
}

export function reportError(severity: ErrorSeverity, context: string, error: Error | string, logger?: DriverLogger, links?: string[]) {
  const err = error instanceof Error ? error : new Error(error);
  if (logger)
    logger.logError(severity, context, err, links);
  else
    reportErrorToUI(severity, context, err);
}
