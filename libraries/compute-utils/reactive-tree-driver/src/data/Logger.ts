import {v4 as uuidv4} from 'uuid';
import {NodeAddress} from './BaseTree';
import {BehaviorSubject} from 'rxjs';

export type DebugLogType =
 'linkRunScheduled' | 'linkRunFinished' |
 'treeUpdateStarted' | 'treeUpdateMutation' | 'treeUpdateFinished'

export interface DebugLogBase {
  uuid: string,
  timestamp: Date,
}

export interface LinkLogPayload {
  linkUUID: string,
  prefix: NodeAddress,
  id: string,
}

export interface LinkLogItem extends DebugLogBase, LinkLogPayload {}

export interface LinkRunStartedLogItem extends LinkLogItem {
  type: 'linkRunStarted',
}

export interface LinkRunFinishedLogItem extends LinkLogItem {
  type: 'linkRunFinished',
}

export interface TreeUpdateStartedLogItem extends DebugLogBase {
  type: 'treeUpdateStarted',
}

export interface TreeUpdateMutationPayload {
  mutationRootPath?: NodeAddress,
  addIdx?: number,
  removeIdx?: number,
}

export interface TreeUpdateMutationLogItem extends DebugLogBase, TreeUpdateMutationPayload {
  type: 'treeUpdateMutation',
}

export interface TreeUpdateFinishedLogItem extends DebugLogBase {
  type: 'treeUpdateFinished',
}

export type LogItem =
TreeUpdateStartedLogItem | TreeUpdateMutationLogItem | TreeUpdateFinishedLogItem |
LinkRunStartedLogItem | LinkRunFinishedLogItem;

export class DriverLogger {
  private log: LogItem[] = [];
  public logs$ = new BehaviorSubject<LogItem[]>([]);

  logLink(type: 'linkRunStarted' | 'linkRunFinished', data: LinkLogPayload) {
    const uuid = uuidv4();
    const timestamp = new Date();
    const logItem = {type, uuid, timestamp, ...data};
    this.log = [...this.log, logItem];
    this.logs$.next(this.log);
  }

  logTreeUpdates(type: 'treeUpdateStarted' | 'treeUpdateFinished') {
    const uuid = uuidv4();
    const timestamp = new Date();
    this.log = [...this.log, {type, uuid, timestamp}];
    this.logs$.next(this.log);
  }

  logMutations(data: TreeUpdateMutationPayload = {}) {
    const uuid = uuidv4();
    const timestamp = new Date();
    this.log = [...this.log, {type: 'treeUpdateMutation', uuid, timestamp, ...data}];
  }
}
