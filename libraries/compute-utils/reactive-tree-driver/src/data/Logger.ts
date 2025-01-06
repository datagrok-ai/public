import {v4 as uuidv4} from 'uuid';
import {NodeAddress} from './BaseTree';
import {BehaviorSubject} from 'rxjs';

export type DebugLogType =
 'linkRunScheduled' | 'linkRunFinished' | 'linkAdded' | 'linkRemoved' | 'actionAdded' | 'actionRemoved' |
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
LinkRunStartedLogItem | LinkRunFinishedLogItem | LinkAddedLogItem | LinkRemoveLogItem | ActionAddedLogItem | ActionRemoveLogItem;

export class DriverLogger {
  private log: LogItem[] = [];
  public logs$ = new BehaviorSubject<LogItem[]>([]);

  logLink(type: 'linkRunStarted' | 'linkRunFinished' | 'linkAdded' | 'linkRemoved' | 'actionAdded' | 'actionRemoved', data: LinkLogPayload) {
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
