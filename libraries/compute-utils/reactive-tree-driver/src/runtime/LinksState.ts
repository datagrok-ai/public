import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BaseTree, NodePath, NodePathSegment} from '../data/BaseTree';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, LinkSpec, MatchInfo, matchNodeLink} from './link-matching';
import {Action, Link} from './Link';
import {BehaviorSubject, concat, merge, Subject, of, Observable, defer, combineLatest, identity, EMPTY} from 'rxjs';
import {takeUntil, map, scan, switchMap, filter, mapTo, toArray, take, tap, debounceTime, delay, concatMap, finalize} from 'rxjs/operators';
import {DriverLogger} from '../data/Logger';
import {getLinksDiff} from './links-diff';
import {ViewAction} from '../config/PipelineInstance';
import {calculateStepsDependencies, calculateIoDependencies, createDefaultValidators, DependenciesData, IoDeps} from './links-dependencies';

export interface LinksData {
  uuid: string;
  id: string;
  prefix: Readonly<NodePath>;
  basePath?: Readonly<NodePath>;
  isAction: boolean;
  matchInfo: MatchInfo;
}

export class LinksState {
  private closed$ = new Subject<true>();
  private linksUpdates = new Subject<true>();
  private runnedInit = new Set<string>();

  public links: Map<string, Link> = new Map();
  public actions: Map<string, Action> = new Map();
  public nodesActions: Map<string, Action[]> = new Map();
  public stepsDependencies: Map<string, DependenciesData> = new Map();
  public ioDependencies: Map<string, IoDeps> = new Map();

  public runningLinks$ = new BehaviorSubject<undefined | string[]>(undefined);
  public forceInitialMetaRun = false;
  public batchLinks: boolean;

  private batchPending = new Set<string>();
  private batchTrigger$ = new Subject<void>();

  constructor(
    private defaultValidators: boolean = false,
    private logger?: DriverLogger,
    batchLinks: boolean = false,
  ) {
    this.batchLinks = batchLinks;

    this.linksUpdates.pipe(
      switchMap(() => this.getRunningLinks()),
      takeUntil(this.closed$),
    ).subscribe(this.runningLinks$);

    this.batchTrigger$.pipe(
      debounceTime(0),
      takeUntil(this.closed$),
    ).subscribe(() => {
      const pending = [...this.batchPending];
      this.batchPending.clear();
      for (const uuid of pending)
        this.links.get(uuid)?.trigger();
    });
  }

  public scheduleBatch(linkUUID: string) {
    this.batchPending.add(linkUUID);
    this.batchTrigger$.next();
  }

  public update(state: BaseTree<StateTreeNode>, isMutation: boolean) {
    this.disableLinks();

    const oldLinks = [...this.links.values()];
    const oldActions = [...this.actions.values()];

    const [links, addedLinks] = this.updateLinks(state, oldLinks);
    this.links = new Map(links.map((link) => [link.uuid, link] as const));

    const [actions, nodesActions] = this.updateActions(state, oldActions);
    this.actions = new Map(actions.map((link) => [link.uuid, link] as const));
    this.nodesActions = nodesActions;

    this.stepsDependencies = calculateStepsDependencies(state, links);
    this.ioDependencies = calculateIoDependencies(state, links, this.logger);

    if (isMutation) {
      const newData = addedLinks.filter((link) => this.isDataLink(link));
      const newMeta = addedLinks.filter((link) => !this.isDataLink(link));
      this.linksUpdates.next(true);
      return concat(
        of(this.wireLinks(state)),
        this.runNewInits(state),
        this.runLinks(state, this.toLinksMap(newData), false),
        this.runLinks(state, this.toLinksMap(newMeta), true),
      ).pipe(toArray(), mapTo(undefined));
    } else {
      this.linksUpdates.next(true);
      const metaMap = this.toLinksMap(links.filter((link) => !this.isDataLink(link)));
      const initLinks = this.toLinksMap(links.filter((link) => this.isOnInitDataLink(link)));
      return concat(
        of(this.wireLinks(state)),
        this.runNewInits(state),
        this.runLinks(state, initLinks, false),
        (this.defaultValidators || this.forceInitialMetaRun) ? of(null).pipe(delay(0), concatMap(() => this.runLinks(state, metaMap, true))) : of(null),
      ).pipe(toArray(), mapTo(undefined));
    }
  }

  public getNodeActionsData(uuid: string): ViewAction[] | undefined {
    const actions = this.nodesActions.get(uuid);
    if (actions) {
      return actions.map(({uuid, spec: {id, position, description, menuCategory, friendlyName, icon, confirmationMessage}}) =>
        ({uuid, id, position, description, menuCategory, friendlyName, icon, confirmationMessage}));
    }
    return actions;
  }

  public updateLinks(state: BaseTree<StateTreeNode>, oldLinks: Link[]) {
    const newLinks = this.createStateLinks(state);
    if (this.defaultValidators) {
      const validators = createDefaultValidators(state, this.logger);
      newLinks.push(...validators);
    }
    return this.mergeLinks(oldLinks, newLinks, 'link');
  }

  public updateActions(state: BaseTree<StateTreeNode>, oldActions: Action[]) {
    const newActions = this.createStateActions(state);
    const [mergedActions] = this.mergeLinks(oldActions, newActions, 'action');
    const nodeActions = new Map<string, Action[]>;
    for (const action of mergedActions) {
      const visibleOn = action.spec.visibleOn;
      let targetUuid: string;
      if (visibleOn) {
        // Route action to the descendant node matching visibleOn configId
        const found = state.getNode(action.prefix).find((item) => item.config.id === visibleOn);
        targetUuid = found ? found[0].getItem().uuid : state.getNode(action.prefix).getItem().uuid;
      } else {
        targetUuid = state.getNode(action.prefix).getItem().uuid;
      }
      const acts = nodeActions.get(targetUuid) ?? [];
      acts.push(action);
      nodeActions.set(targetUuid, acts);
    }
    return [mergedActions, nodeActions] as const;
  }

  private mergeLinks<L extends Link>(oldLinks: L[], newLinks: L[], prefix: 'link' | 'action'): [L[], L[]] {
    const addedLinks: L[] = [];
    const {toRemove, toAdd} = getLinksDiff(oldLinks, newLinks);
    const mergedLinks: L[] = [];
    for (const oldLink of oldLinks) {
      if (!toRemove.has(oldLink.uuid))
        mergedLinks.push(oldLink);
      else {
        if (this.logger && !oldLink.matchInfo.isDefaultValidator)
          this.logger.logLink(`${prefix}Removed`, {linkUUID: oldLink.uuid, prefix: oldLink.prefix, basePath: oldLink.matchInfo.basePath, id: oldLink.matchInfo.spec.id});
        oldLink.destroy();
      }
    }
    for (const newLink of newLinks) {
      if (toAdd.has(newLink.uuid)) {
        if (this.logger && !newLink.matchInfo.isDefaultValidator)
          this.logger.logLink(`${prefix}Added`, {linkUUID: newLink.uuid, prefix: newLink.prefix, basePath: newLink.matchInfo.basePath, id: newLink.matchInfo.spec.id});
        mergedLinks.push(newLink);
        addedLinks.push(newLink);
      }
    }
    return [mergedLinks, addedLinks] as const;
  }

  public createStateLinks(state: BaseTree<StateTreeNode>) {
    const links = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (!isFuncCallNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => {
          const spec = minfo.spec;
          const debounce = spec.type === 'validator'
            ? (spec.debounce ?? (this.batchLinks ? 0 : undefined))
            : undefined;
          const link = new Link(path, minfo, debounce, this.logger);
          return link;
        });
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }

  public createStateActions(state: BaseTree<StateTreeNode>) {
    const actions = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      const {config} = item;
      const matchedLinks = (config.actions ?? [])
        .map((link) => matchNodeLink(node, link))
        .filter((x) => !!x)
        .flat();
      const links = matchedLinks.map((minfo) => {
        const spec = minfo.spec as ActionSpec;
        const action = new Action(path, minfo, spec, this.logger);
        return action;
      });
      return [...acc, ...links];
    }, [] as Action[]);
    return actions;
  }

  public runLinks(state: BaseTree<StateTreeNode>, links: Map<string, Link>, isMeta: boolean) {
    const scheduledLinks = new Set<string>();
    const batchableLinks: Link[] = [];
    const sequentialObs: Observable<undefined>[] = [];

    state.traverse(state.root, (_acc, node) => {
      const item = node.getItem();
      const deps = this.stepsDependencies.get(item.uuid);
      const toRunLinks = [...(deps?.links ?? [])]
        .filter((linkUUID) => {
          const link = links.get(linkUUID);
          return link && !scheduledLinks.has(linkUUID);
        });
      const priorityOrderedLinks = toRunLinks.sort((a, b) => {
        const l1 = links.get(a)!;
        const l2 = links.get(b)!;
        return (l2.matchInfo.spec.nodePriority ?? 0) - (l1.matchInfo.spec.nodePriority ?? 0);
      });

      for (const linkUUID of priorityOrderedLinks) {
        scheduledLinks.add(linkUUID);
        const link = links.get(linkUUID)!;
        if (this.batchLinks && link.isBatchable) {
          batchableLinks.push(link);
        } else {
          sequentialObs.push(
            combineLatest([
              this.runningLinks$,
              this.getLinkRunObs(linkUUID),
            ]).pipe(
              isMeta ? identity : debounceTime(0),
              filter(([running]) => !running || running.length === 0),
              take(1),
              mapTo(undefined),
            ),
          );
        }
      }
      return _acc;
    }, undefined as void);

    const batchObs = this.runLinksBatched(batchableLinks);
    return concat(batchObs, ...sequentialObs);
  }

  private runLinksBatched(links: Link[]) {
    if (links.length === 0) return EMPTY;
    return defer(() => {
      for (const link of links)
        link.trigger();
      return this.runningLinks$.pipe(
        filter((x) => !x || x.length === 0),
        take(1),
        mapTo(undefined),
      );
    });
  }

  public waitForLinks() {
    return this.runningLinks$.pipe(
      filter((links) => links == null || links?.length === 0),
      debounceTime(0),
      take(1),
    );
  }

  public runNewInits(state: BaseTree<StateTreeNode>) {
    const obs = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (!isFuncCallNode(item) && !this.runnedInit.has(item.uuid) && item.config.onInit) {
        const minfo = matchNodeLink(node, item.config.onInit, path);
        if (!minfo)
          return acc;
        const initLink = new Link(path, minfo[0], undefined, this.logger);
        const obs$ = defer(() => {
          initLink.wire(state, this);
          this.runnedInit.add(item.uuid);
          initLink.trigger();
          return initLink.isRunning$.pipe(
            filter((x) => !x),
            take(1),
            tap(() => initLink.destroy()),
            mapTo(undefined),
          );
        });
        return [...acc, obs$];
      }
      return acc;
    }, [] as Observable<undefined>[]);
    return concat(...obs);
  }

  public runReturnResults(state: BaseTree<StateTreeNode>) {
    const obs$ = defer(() => {
      const item = state.root.getItem();
      if (isFuncCallNode(item) || !item.config.onReturn)
        return of(undefined);
      const minfo = matchNodeLink(state.root, item.config.onReturn, []);
      if (!minfo)
        return of(undefined);
      const returnLink = new Link([], minfo[0], undefined, this.logger);
      returnLink.wire(state, this);
      returnLink.trigger();
      return returnLink.isRunning$.pipe(
        filter((x) => !x),
        take(1),
        map(() => returnLink.returnResult),
        finalize(() => returnLink.destroy()),
      );
    });
    return obs$;
  }


  public wireLinks(state: BaseTree<StateTreeNode>) {
    for (const [, link] of this.links) {
      link.wire(state, this);
      link.setActive();
    }
    for (const [, action] of this.actions)
      action.wire(state, this);
  }

  public disableLinks() {
    for (const [, link] of this.links)
      link.setInactive();
  }

  public close() {
    this.closed$.next(true);
  }

  public getLinksInfo(): LinksData[] {
    const links = [...this.links.values()].filter((l) => !l.matchInfo.isDefaultValidator).map((l) => ({id: l.matchInfo.spec.id, uuid: l.uuid, prefix: l.prefix, basePath: l.matchInfo.basePath, isAction: false, matchInfo: l.matchInfo}));
    const actions = [...this.actions.values()].map((l) => ({id: l.matchInfo.spec.id, uuid: l.uuid, prefix: l.prefix, basePath: l.matchInfo.basePath, isAction: true, matchInfo: l.matchInfo}));
    return [...links, ...actions];
  }

  private getRunningLinks() {
    const obs = [...this.links.values(), ...[...this.nodesActions.values()].flat()].map(
      (link) => link.isRunning$.pipe(map((isRunning) => [link.uuid, isRunning] as const)));
    return merge(...obs).pipe(
      scan((acc, [uuid, isRunning]) => {
        acc[uuid] = isRunning;
        return acc;
      }, {} as Record<string, boolean>),
      map((data) => Object.entries(data).filter(
        ([, isRunning]) => isRunning).map(([uuid]) => uuid)),
    );
  }

  public isDataLink(link: Link) {
    return !link.matchInfo.spec.type || link.matchInfo.spec.type === 'data';
  }

  public isOnInitDataLink(link: Link) {
    return (!link.matchInfo.spec.type || link.matchInfo.spec.type === 'data') && link.matchInfo.spec.runOnInit;
  }

  public isDefaultValidatorLink(link: Link) {
    return link.matchInfo.isDefaultValidator;
  }

  private getLinkRunObs(linkUUID: string) {
    const link = this.links.get(linkUUID)!;
    const obs$ = defer(() => {
      link.trigger();
      return this.runningLinks$.pipe(filter((x) => x?.length === 0), mapTo(true as const));
    });
    return obs$;
  }

  private toLinksMap(links: Link[]) {
    return new Map(links.map((link) => [link.uuid, link]));
  }
}
