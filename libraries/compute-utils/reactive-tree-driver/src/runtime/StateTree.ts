import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, defer, of, merge, Subject, BehaviorSubject, from, combineLatest} from 'rxjs';
import {finalize, map, mapTo, toArray, concatMap, tap, takeUntil, debounceTime, scan, withLatestFrom, filter} from 'rxjs/operators';
import {NodePath, BaseTree, TreeNode} from '../data/BaseTree';
import {PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {PipelineInstanceConfig, PipelineSerializedState, PipelineState} from '../config/PipelineInstance';
import {isPipelineStepConfig} from '../config/config-utils';
import {FuncCallAdapter} from './FuncCallAdapters';
import {saveFuncCall, saveInstanceState} from './funccall-utils';
import {ConsistencyInfo, FuncCallNode, FuncCallStateInfo, isFuncCallNode, StateTreeNode, StateTreeSerializationOptions} from './StateTreeNodes';
import {indexFromEnd} from '../utils';
import {LinksState} from './LinksState';
import {GranularMutationOp, ValidationResult} from '../data/common-types';
import {DriverLogger, TreeUpdateMutationPayload} from '../data/Logger';
import {ItemMetadata} from '../view/ViewCommunication';
import * as Factory from './StateTreeFactory';
import * as Serializer from './StateTreeSerializer';

const MAX_CONCURENT_SAVES = 5;

export interface TreeUpdateData {
  isMutation: boolean;
  details?: TreeUpdateMutationPayload[],
}

export class LockError extends Error {}

export class StateTree {
  private closed$ = new Subject<true>();

  public nodeTree: BaseTree<StateTreeNode>;
  public nodesMap: Map<string, StateTreeNode> = new Map();
  public linksState: LinksState;

  public makeStateRequests$ = new Subject<true>();
  public globalROLocked$ = new BehaviorSubject(false);
  public treeMutationsLocked$ = new BehaviorSubject(false);
  public result$ = new Subject<any>();

  constructor(
    item: StateTreeNode,
    public config: PipelineConfigurationProcessed,
    private mockMode = false,
    private defaultValidators = false,
    private logger?: DriverLogger,
    private batchLinks = false,
  ) {
    this.nodeTree = new BaseTree(item);
    this.linksState = new LinksState(defaultValidators, this.logger, batchLinks);

    this.linksState.runningLinks$.pipe(
      map((links) => !!links?.length),
      takeUntil(this.closed$),
    ).subscribe(this.treeMutationsLocked$);
  }

  //
  // tree state getting (delegates to StateTreeSerializer)
  //

  public toSerializedState(options: StateTreeSerializationOptions = {}): PipelineSerializedState {
    return StateTree.toSerializedStateRec(this.nodeTree.root, options);
  }

  public toState(options: StateTreeSerializationOptions = {}): PipelineState {
    return StateTree.toStateRec(this.nodeTree.root, options, this.linksState);
  }

  // Static re-exports for backward compatibility
  public static toStateRec = Serializer.toStateRec;
  public static toSerializedStateRec = Serializer.toSerializedStateRec;

  public getValidations() {
    return Serializer.getValidations(this.nodeTree);
  }

  public getConsistency() {
    return Serializer.getConsistency(this.nodeTree);
  }

  public getMeta() {
    return Serializer.getMeta(this.nodeTree);
  }

  public getFuncCallStates() {
    return Serializer.getFuncCallStates(this.nodeTree);
  }

  public getNodesDescriptions() {
    return Serializer.getNodesDescriptions(this.nodeTree);
  }

  public getIOMutations() {
    return Serializer.getIOMutations(this.nodeTree);
  }

  //
  // events handling
  //

  public init() {
    return this.mutateTree(() => {
      return StateTree.loadOrCreateCalls(this, this.mockMode).pipe(mapTo([]));
    }, false).pipe(
      concatMap(() => this.linksState.waitForLinks()),
      mapTo(this),
    );
  }

  // for testing only
  public runMutateTree() {
    return this.mutateTree(() => of([{isMutation: true}]));
  }

  public save(uuid?: string, metaData?: ItemMetadata) {
    return this.mutateTree(() => {
      const [root, nqName] = StateTree.findPipelineNode(this, uuid);
      if (this.mockMode)
        throw new Error(`Cannot save mock tree`);
      if (nqName == null)
        throw new Error(`Attempting to save pipeline with no nqName`);
      const pendingFuncCallSaves = this.nodeTree.traverse(root, (acc, node) => {
        const item = node.getItem();
        if (isFuncCallNode(item)) {
          const obs$ = defer(() => saveFuncCall(item.instancesWrapper)).pipe(
            map((nextCall) => {
              const nadapter = new FuncCallAdapter(nextCall, item.isReadonly);
              item.changeAdapter(nadapter);
            }),
          );
          return [...acc, obs$];
        }
        return acc;
      }, [] as Array<Observable<void>>);
      return merge(...pendingFuncCallSaves, MAX_CONCURENT_SAVES).pipe(
        toArray(),
        concatMap(() => this.saveMetaCall(root, nqName, metaData)),
        map((call) => [undefined, call]),
      );
    });
  }

  // --- Internal unguarded primitives (caller must hold the tree lock) ---

  private addSubTreeInternal(parentPath: NodePath, id: string, pos: number): TreeUpdateMutationPayload {
    const subConfig = StateTree.getSubConfig(this.config, parentPath, id);
    if (isPipelineStepConfig(subConfig)) {
      const fcNode = new FuncCallNode(subConfig, false, this.logger);
      fcNode.initState(subConfig);
      this.nodeTree.attachBrunch(parentPath, new TreeNode(fcNode), id, pos);
    } else {
      StateTree.fromPipelineConfig({
        config: this.config,
        startNode: subConfig,
        startPath: [...parentPath, {id, idx: pos}],
        startState: this,
        isReadonly: false,
        defaultValidators: this.defaultValidators,
        batchLinks: this.batchLinks,
        mockMode: this.mockMode,
        logger: this.logger,
      });
    }
    return {mutationRootPath: parentPath, addIdx: pos, id};
  }

  private removeSubtreeInternal(uuid: string): TreeUpdateMutationPayload {
    const data = this.nodeTree.find((item) => item.uuid === uuid);
    if (data == null)
      throw new Error(`Node uuid ${uuid} not found`);
    const [node, path] = data;
    this.nodeTree.removeBrunch(path);
    const [mutationRootPath, removeIdx] = this.getMutationSlice(path);
    return {mutationRootPath, removeIdx, id: node.getItem().config.id};
  }

  private moveSubtreeInternal(uuid: string, newIdx: number): TreeUpdateMutationPayload {
    const data = this.nodeTree.find((item) => item.uuid === uuid);
    if (data == null)
      throw new Error(`Node uuid ${uuid} not found`);
    const [node, path] = data;
    this.nodeTree.removeBrunch(path);
    const [ppath, oldIdx] = this.getMutationSlice(path);
    this.nodeTree.attachBrunch(ppath, node, node.getItem().config.id, newIdx);
    return {mutationRootPath: ppath, addIdx: newIdx, removeIdx: oldIdx, id: node.getItem().config.id};
  }

  private replaceSubtreeInternal(path: NodePath, initConfig: PipelineInstanceConfig): {detail: TreeUpdateMutationPayload, subTree: StateTree} {
    const ppath = path.slice(0, -1);
    const last = indexFromEnd(path);
    const subConfig = last ? StateTree.getSubConfig(this.config, ppath, last.id) : this.config;
    if (isPipelineStepConfig(subConfig))
      throw new Error(`FuncCall node ${JSON.stringify(path)}, but pipeline is expected`);
    const subTree = StateTree.fromInstanceConfig({
      instanceConfig: initConfig,
      config: subConfig,
      isReadonly: false,
      defaultValidators: this.defaultValidators,
      batchLinks: this.batchLinks,
      mockMode: this.mockMode,
      logger: this.logger,
    });
    if (last) {
      this.nodeTree.removeBrunch(path);
      this.nodeTree.attachBrunch(ppath, subTree.nodeTree.root, last.id, last.idx);
    } else
      this.nodeTree.replaceRoot(subTree.nodeTree.root);

    const detail: TreeUpdateMutationPayload = last ? {
      mutationRootPath: ppath,
      addIdx: last.idx,
      id: last.id,
    } : {mutationRootPath: [], id: subTree.nodeTree.root.getItem().config.id};
    return {detail, subTree};
  }

  // --- Public guarded methods ---

  public moveSubtree(uuid: string, newIdx: number) {
    return this.mutateTree(() => {
      const detail = this.moveSubtreeInternal(uuid, newIdx);
      return of([{isMutation: true, details: [detail]}]);
    });
  }

  public removeSubtree(uuid: string) {
    return this.mutateTree(() => {
      const detail = this.removeSubtreeInternal(uuid);
      return of([{isMutation: true, details: [detail]}]);
    });
  }

  public addSubTree(puuid: string, id: string, pos: number) {
    return this.mutateTree(() => {
      const [_root, _nqName, ppath] = StateTree.findPipelineNode(this, puuid);
      const detail = this.addSubTreeInternal(ppath, id, pos);
      return StateTree.loadOrCreateCalls(this, this.mockMode).pipe(mapTo([{isMutation: true, details: [detail]}]));
    });
  }

  public loadSubTree(puuid: string, dbId: string, id: string, pos: number, isReadonly: boolean, isReplace = false) {
    return this.mutateTree(() => {
      const [_root, _nqName, ppath] = StateTree.findPipelineNode(this, puuid);
      const subConfig = StateTree.getSubConfig(this.config, ppath, id);
      if (isPipelineStepConfig(subConfig))
        throw new Error(`FuncCall node ${JSON.stringify(ppath)}, but pipeline is expected`);
      const tree = StateTree.load(
        {
          dbId,
          config: subConfig,
          mockMode: this.mockMode,
          defaultValidators: this.defaultValidators,
          batchLinks: this.batchLinks,
          isReadonly,
        },
      ).pipe(
        concatMap((tree) => StateTree.loadOrCreateCalls(tree, this.mockMode)),
      ).pipe(
        map((tree) => isReplace ?
          this.nodeTree.replaceBrunch(ppath, tree.nodeTree.root, id, pos) :
          this.nodeTree.attachBrunch(ppath, tree.nodeTree.root, id, pos),
        ),
      );
      const details: TreeUpdateMutationPayload[] =[{
        mutationRootPath: ppath,
        addIdx: pos,
        id,
      }];
      return tree.pipe(mapTo([{isMutation: true, details}]));
    });
  }

  private applyGranularOps(parentPath: NodePath, ops: GranularMutationOp[]): Observable<TreeUpdateMutationPayload> {
    const details: TreeUpdateMutationPayload[] = [];
    for (const op of ops) {
      if (op.op === 'add') {
        const pos = op.position ?? this.nodeTree.getNode(parentPath).getChildren().length;
        details.push(this.addSubTreeInternal(parentPath, op.configId, pos));
      } else if (op.op === 'remove')
        details.push(this.removeSubtreeInternal(op._uuid));
      else if (op.op === 'move')
        details.push(this.moveSubtreeInternal(op._uuid, op.position));
    }
    return StateTree.loadOrCreateCalls(this, this.mockMode).pipe(
      concatMap(() => from(details)),
    );
  }

  public runStep(uuid: string, mockResults?: Record<string, any>, mockDelay?: number) {
    return this.withTreeLock(() => {
      if (!this.mockMode && mockResults)
        throw new Error(`Mock results passed to runStep while mock mode is off`);
      const res = this.nodeTree.find((item) => item.uuid === uuid);
      if (!res)
        throw new Error(`Step uuid ${uuid} not found`);
      const [snode] = res;
      const node = snode.getItem();
      if (!isFuncCallNode(node))
        throw new Error(`Step uuid ${uuid} is not FuncCall`);
      return node.instancesWrapper.run(mockResults, mockDelay).pipe(mapTo(undefined));
    });
  }

  public runSequence(startUuid: string, rerunWithConsistent?: boolean, includeNonNested?: boolean) {
    return this.withTreeLock(() => {
      const startNode = includeNonNested ? this.nodeTree.root : this.nodeTree.find((item) => item.uuid === startUuid)?.[0];
      if (startNode == null)
        return of(undefined);
      const nodesSeq = this.nodeTree.traverse(startNode, (acc, node) => [...acc, node.getItem()], [] as StateTreeNode[]);
      const startIdx = nodesSeq.findIndex((node) => node.uuid === startUuid);
      if (startIdx < 0)
        return of(undefined);
      return from(nodesSeq.slice(startIdx)).pipe(
        concatMap((node) => {
          if (!isFuncCallNode(node) || node.pendingDependencies$.value?.length)
            return of(undefined);
          if (rerunWithConsistent) {
            return node.getStateStore().overrideToConsistent().pipe(
              concatMap(() => this.linksState.waitForLinks()),
              withLatestFrom(node.getStateStore().isRunable$, node.getStateStore().isOutputOutdated$),
              filter(([, runable, outdated]) => runable && outdated),
              concatMap(() => node.getStateStore().run()),
              concatMap(() => this.linksState.waitForLinks()),
            );
          } else {
            if (!node.getStateStore().isRunable$.value || !node.getStateStore().isOutputOutdated$.value)
              return of(undefined);
            return node.getStateStore().run().pipe(
              concatMap(() => this.linksState.waitForLinks()),
            );
          }
        }),
        toArray(),
        mapTo(undefined),
      );
    });
  }

  public runAction(uuid: string, additionalParams: Record<string, any> = {}) {
    const action = this.linksState.actions.get(uuid);
    if (!action)
      throw new Error(`Action ${uuid} not found`);
    if (action.spec.type === 'pipeline') {
      return this.mutateTree(() => {
        return action.execPipelineMutations(additionalParams).pipe(
          concatMap(({pipelineMutations, granularMutations}) => {
            const replacements$ = from(pipelineMutations ?? []).pipe(
              concatMap((data) => {
                const {detail, subTree} = this.replaceSubtreeInternal(data.path, data.initConfig);
                return StateTree.loadOrCreateCalls(subTree, this.mockMode).pipe(mapTo(detail));
              }),
            );

            const granular$ = from(granularMutations ?? []).pipe(
              concatMap((data) => this.applyGranularOps(data.path, data.ops)),
            );

            return merge(replacements$, granular$);
          }),
          toArray(),
          map((details) => [{isMutation: true, details}]),
        );
      });
    } else if (action.spec.type === 'funccall') {
      return this.withTreeLock(() => action.exec(additionalParams).pipe(
        finalize(() => this.makeStateRequests$.next(true)),
      ));
    } else
      return action.exec(additionalParams);
  }

  public returnResult() {
    return this.withTreeLock(() => {
      return this.linksState.runReturnResults(this.nodeTree).pipe(
        tap((res) => {
          this.result$.next(res);
        }),
      );
    });
  }

  public updateFuncCall(uuid: string, call: DG.FuncCall) {
    return this.withTreeLock(() => {
      const data = this.nodeTree.find((item) => item.uuid === uuid);
      if (!data)
        throw new Error(`No FuncCall node found ${uuid}`);
      const [node, path] = data;
      const item = node.getItem();
      if (!isFuncCallNode(item) || item.instancesWrapper.isReadonly)
        throw new Error(`FuncCall writable node is expected on path ${JSON.stringify(path)}`);
      const adapter = new FuncCallAdapter(call, false);
      item.changeAdapter(adapter, true);
      item.setOutdatedStatus(false);
      return of(item);
    }).pipe(
      finalize(() => this.makeStateRequests$.next(true)),
    );
  }

  public resetToConsistent(uuid: string, ioName: string) {
    return this.withTreeLock(() => {
      const data = this.nodeTree.find((item) => item.uuid === uuid);
      if (!data)
        return of(undefined);
      const [node] = data;
      const item = node.getItem();
      if (!isFuncCallNode(item))
        return of(undefined);
      item.instancesWrapper.setToConsistent(ioName);
      return of(undefined);
    });
  }

  public close() {
    this.nodeTree.traverse(this.nodeTree.root, (acc, node) => {
      const item = node.getItem();
      if (isFuncCallNode(item))
        item.close();
      return acc;
    }, undefined);
    this.closed$.next(true);
  }

  //
  // Static factory re-exports (backward compatibility)
  //

  static load = Factory.loadStateTree;
  static fromPipelineConfig = Factory.fromPipelineConfig;
  static fromInstanceState = Factory.fromPipelineInstanceState;
  static fromInstanceConfig = Factory.fromPipelineInstanceConfig;
  static loadOrCreateCalls = Factory.loadOrCreateCalls;
  static findPipelineNode = Factory.findPipelineNode;
  static getSubConfig = Factory.getSubConfig;

  //
  // locking, tree mutation and deps tracking
  //

  private mutateTree<R>(fn: () => Observable<readonly [TreeUpdateData?, R?]>, waitForLinks = true, disableLinks = true) {
    return defer(() => {
      if (this.logger)
        this.logger.logTreeUpdates('treeUpdateStarted');
      this.treeLock();
      return (waitForLinks ? this.linksState.waitForLinks() : of(undefined)).pipe(
        tap(() => disableLinks ? this.linksState.disableLinks() : undefined),
        concatMap(() => fn()),
      );
    }).pipe(
      tap(() => this.updateNodesMap()),
      concatMap((data) => {
        const [mutationData, res] = data ?? [];
        if (this.logger && mutationData?.details?.length) {
          for (const item of mutationData?.details)
            this.logger.logMutations(item);
        }
        const isMutation = mutationData?.isMutation ?? false;
        return this.linksState.update(this.nodeTree, isMutation).pipe(mapTo(res));
      }),
      tap(() => {
        this.removeOrphanedIOMetadata();
        this.setDepsTracker();
      }),
      finalize(() => {
        this.treeUnlock();
        if (this.logger)
          this.logger.logTreeUpdates('treeUpdateFinished');
        this.makeStateRequests$.next(true);
      }),
    );
  }

  private withTreeLock<T = undefined>(fn: () => Observable<T>) {
    return defer(() => {
      this.treeLock();
      return this.linksState.waitForLinks().pipe(concatMap(() => fn()));
    }).pipe(
      finalize(() => {
        this.treeUnlock();
      }),
    );
  }

  private updateNodesMap() {
    this.nodesMap = new Map();
    this.nodeTree.traverse(this.nodeTree.root, (acc, node) => {
      const item = node.getItem();
      this.nodesMap.set(item.uuid, item);
      return acc;
    }, undefined);
  }

  private removeOrphanedIOMetadata() {
    const currentLinkIds = new Set(this.linksState.links.keys());

    this.nodeTree.traverse(this.nodeTree.root, (acc, node) => {
      const item = node.getItem();

      if (!isFuncCallNode(item)) {
        item.clearOldTags(currentLinkIds);
        return acc;
      }

      item.clearOldValidations(currentLinkIds);
      item.clearOldMetas(currentLinkIds);

      const ioDeps = this.linksState.ioDependencies.get(item.uuid) ?? {};
      const ioNames = item.instancesWrapper.getStateNames();

      for (const ioName of ioNames) {
        const deps = ioDeps[ioName];
        if (!deps?.data)
          item.clearIORestriction(ioName);
      }

      return acc;
    }, null);
  }

  private setDepsTracker() {
    this.nodeTree.traverse(this.nodeTree.root, (acc, node) => {
      const item = node.getItem();
      if (!isFuncCallNode(item))
        return acc;

      const deps = this.linksState.stepsDependencies.get(item.uuid);
      const depsStates = [...(deps?.nodes ?? [])].filter((depId) => {
        const depItem = this.nodesMap.get(depId);
        return depItem && isFuncCallNode(depItem);
      }).map((depId) => {
        const depItem = this.nodesMap.get(depId)! as FuncCallNode;
        const hasPending$ = combineLatest([
          depItem.instancesWrapper.isOutputOutdated$,
          depItem.pendingDependencies$.pipe(
            map((d) => d.length !== 0),
          ),
        ]).pipe(map(([isOutdated, hasPending]) => isOutdated || hasPending));
        return [depId, hasPending$] as const;
      });
      item.setDeps(depsStates);
      return acc;
    }, null);
  }

  private treeLock() {
    if (this.globalROLocked$.value)
      throw new LockError(`Changes double lock`);
    this.globalROLocked$.next(true);
  }

  private treeUnlock() {
    if (this.globalROLocked$.value)
      this.globalROLocked$.next(false);
  }

  private getMutationSlice(path: NodePath) {
    if (path.length === 0)
      return [path, 0] as const;

    const ppath = path.slice(0, -1);
    const endSegment = indexFromEnd(path)!;
    return [ppath, endSegment.idx] as const;
  }

  // meta call saving

  private saveMetaCall(
    root: TreeNode<StateTreeNode>,
    nqName: string,
    metaData?: ItemMetadata,
  ) {
    return defer(() => {
      if (this.mockMode)
        return of(undefined);
      const state = StateTree.toSerializedStateRec(root, {disableNodesUUID: true});
      const item = root.getItem();
      const version = !isFuncCallNode(item) ? item.config.version : undefined;
      return saveInstanceState(nqName, state, metaData, version);
    });
  }
}
