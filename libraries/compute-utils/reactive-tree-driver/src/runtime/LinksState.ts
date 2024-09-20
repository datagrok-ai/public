import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BaseTree, NodeAddress, NodePath} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {isFuncCallNode, isSequentialPipelineNode, isStaticPipelineNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, isActionSpec, MatchedNodePaths, matchNodeLink} from './link-matching';
import {Action, Link} from './Link';
import {BehaviorSubject, concat, merge, Subject, of, Observable, defer} from 'rxjs';
import {takeUntil, map, scan, switchMap, filter, mapTo, toArray, take, tap} from 'rxjs/operators';

export class DependenciesData {
  nodes: Set<string> = new Set();
  links: Set<string> = new Set();
}

export class LinksState {
  private closed$ = new Subject<true>();
  private linksUpdates = new Subject<true>();
  private runnedInit = new Set<string>();

  public links: Map<string, Link> = new Map();
  public actions: Map<string, Action> = new Map();
  public nodesActions: Map<string, Action[]> = new Map();
  public deps: Map<string, DependenciesData> = new Map();

  public runningLinks$ = new BehaviorSubject<undefined | string[]>(undefined);

  constructor() {
    this.linksUpdates.pipe(
      switchMap(() => this.getRunningLinks()),
      takeUntil(this.closed$),
    ).subscribe(this.runningLinks$);
  }

  // TODO: try keeping unaffected links (?)
  public update(state: BaseTree<StateTreeNode>, mutationPath?: NodePath, childOffset?: number) {
    this.destroyLinks();
    const links = this.createAutoLinks(state);
    this.links = new Map(links.map((link) => [link.uuid, link] as const));
    [this.actions, this.nodesActions] = this.createActionLinks(state);
    const deps = this.calculateDependencies(state, links);
    this.deps = deps;
    if (mutationPath) {
      const inbound = links.filter((link) => this.isInbound(mutationPath, link, childOffset));
      const outgoing = links.filter((link) => this.isOutgoing(mutationPath, link, childOffset));
      this.checkDisjoint(inbound, outgoing, mutationPath, childOffset);
      this.linksUpdates.next(true);
      return concat(
        this.runNewInits(state),
        this.runReadyInbound(state, mutationPath, childOffset),
        this.runReadyOutgoing(state, mutationPath, childOffset),
        this.runAffectedMetaLinks(state, mutationPath, childOffset),
        this.wireLinks(state),
      ).pipe(toArray(), mapTo(undefined));
    } else {
      this.linksUpdates.next(true);
      return this.wireLinks(state);
    }
  }

  public createAutoLinks(state: BaseTree<StateTreeNode>) {
    const links = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (isStaticPipelineNode(item) || isSequentialPipelineNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => {
          const link = new Link(path, minfo);
          return link;
        });
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }

  public createActionLinks(state: BaseTree<StateTreeNode>) {
    const actionEntries = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      const {config} = item;
      const matchedLinks = (config.actions ?? [])
        .map((link) => matchNodeLink(node, link))
        .filter((x) => !!x)
        .flat();
      const links = matchedLinks.map((minfo) => {
        const spec = minfo.spec as ActionSpec;
        const action = new Action(path, minfo, spec.position, spec.friendlyName, spec.menuCategory);
        return [item.uuid, action] as const;
      });
      return [...acc, ...links];
    }, [] as (readonly [string, Action])[]);
    const nodeActions = new Map<string, Action[]>;
    for (const [uuid, action] of actionEntries) {
      const acts = nodeActions.get(uuid) ?? [];
      acts.push(action);
      nodeActions.set(uuid, acts);
    }
    const actionsMap = new Map(actionEntries.map(([,action]) => [action.uuid, action]));
    return [actionsMap, nodeActions] as const;
  }

  // TODO: detection of cycles and double io deps, including self deps
  public calculateDependencies(state: BaseTree<StateTreeNode>, links: Link[]) {
    const deps: Map<string, DependenciesData> = new Map();
    for (const link of links) {
      for (const infosIn of Object.values(link.matchInfo.inputs)) {
        for (const infoIn of infosIn) {
          const inPathFull = [...link.prefix, ...infoIn.path];
          const nodeIn = state.getNode(inPathFull);
          // pipeline memory state should be immutable and set in
          // onInit, so they are always ready
          if (!isFuncCallNode(nodeIn.getItem()))
            continue;
          for (const infosOut of Object.values(link.matchInfo.outputs)) {
            for (const infoOut of infosOut) {
              const outPathFull = [...link.prefix, ...infoOut.path];
              const nodeOut = state.getNode(outPathFull);
              const depData = deps.get(nodeOut.getItem().uuid) ?? new DependenciesData();
              if (!BaseTree.isNodeAddressEq(inPathFull, outPathFull))
                depData.nodes.add(nodeIn.getItem().uuid);
              depData.links.add(link.uuid);
              deps.set(nodeOut.getItem().uuid, depData);
            }
          }
        }
      }
    }
    return deps;
  }

  public isLinkReady(state: BaseTree<StateTreeNode>, linkUUID: string) {
    const link = this.links.get(linkUUID);
    if (!link)
      return false;
    for (const infosIn of Object.values(link.matchInfo.inputs)) {
      for (const infoIn of infosIn) {
        const inPathFull = [...link.prefix, ...infoIn.path];
        const nodeIn = state.getNode(inPathFull);
        const item = nodeIn.getItem();
        if (!isFuncCallNode(item))
          continue;
        if (item.instancesWrapper.isOutputOutdated$.value)
          return false;
      }
    }
    return true;
  }

  public isMetaLink(link: Link) {
    return !isActionSpec(link.matchInfo.spec) && (link.matchInfo.spec.isValidator || link.matchInfo.spec.isMeta);
  }

  public runReadyInbound(state: BaseTree<StateTreeNode>, mutationPath: NodePath, childOffset?: number) {
    const root = state.getNode(mutationPath);
    const processedLinks = new Set<string>();
    const obs = state.traverse(root, (acc, node) => {
      const item = node.getItem();
      const deps = this.deps.get(item.uuid);
      const linkRuns = [...(deps?.links ?? [])].filter((linkUUID) => {
        const link = this.links.get(linkUUID);
        return link && !processedLinks.has(linkUUID) && this.isInbound(mutationPath, link, childOffset) && this.isLinkReady(state, linkUUID);
      }).map((linkUUID) => {
        processedLinks.add(linkUUID);
        return this.getLinkRunObs(linkUUID, mutationPath, childOffset);
      });
      return [...acc, ...linkRuns];
    }, [] as Observable<true>[]);
    return concat(...obs);
  }

  public runReadyOutgoing(state: BaseTree<StateTreeNode>, mutationPath: NodePath, childOffset?: number) {
    const processedLinks = new Set<string>();
    const obs = state.traverse(state.root, (acc, node, currentPath) => {
      const item = node.getItem();
      if (BaseTree.isNodeChildOffseted(mutationPath, currentPath, childOffset))
        return acc;

      const deps = this.deps.get(item.uuid);
      const linkRuns = [...(deps?.links ?? [])].filter((linkUUID) => {
        const link = this.links.get(linkUUID);
        return link && !processedLinks.has(linkUUID) && this.isOutgoing(mutationPath, link, childOffset) && this.isLinkReady(state, linkUUID);
      }).map((linkUUID) => {
        processedLinks.add(linkUUID);
        return this.getLinkRunObs(linkUUID, mutationPath, childOffset);
      });
      return [...acc, ...linkRuns];
    }, [] as Observable<true>[]);
    return concat(...obs);
  }

  public runAffectedMetaLinks(state: BaseTree<StateTreeNode>, mutationPath: NodePath, childOffset?: number) {
    const processedLinks = new Set<string>();
    const obs = state.traverse(state.root, (acc, node) => {
      const item = node.getItem();
      const deps = this.deps.get(item.uuid);
      const linkRuns = [...(deps?.links ?? [])].filter((linkUUID) => {
        const link = this.links.get(linkUUID);
        return link && !processedLinks.has(linkUUID) && this.isMetaLink(link) && this.isAffected(mutationPath, link, childOffset) && this.isLinkReady(state, linkUUID);
      }).map((linkUUID) => {
        processedLinks.add(linkUUID);
        return this.getLinkRunObs(linkUUID, mutationPath, childOffset);
      });
      return [...acc, ...linkRuns];
    }, [] as Observable<true>[]);
    return concat(...obs);
  }

  public runNewInits(state: BaseTree<StateTreeNode>) {
    const obs = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (!isFuncCallNode(item) && !this.runnedInit.has(item.uuid) && item.config.onInit) {
        const minfo = matchNodeLink(node, item.config.onInit, path);
        if (!minfo)
          return acc;
        const initLink = new Link(path, minfo[0]);
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

  public wireLinks(state: BaseTree<StateTreeNode>) {
    for (const [, link] of this.links) {
      link.wire(state, this);
      link.setActive();
    }
    for (const [, action] of this.actions) {
      action.wire(state, this);
    }
    return of(undefined);
  }

  public destroyLinks() {
    for (const [, link] of this.links)
      link.destroy();
  }

  public close() {
    this.closed$.next(true);
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

  public isAffected(rootPath: Readonly<NodeAddress>, link: Link, childOffset?: number) {
    return this.hasNested(rootPath, link.prefix, link.matchInfo.inputs, childOffset) ||
      this.hasNested(rootPath, link.prefix, link.matchInfo.outputs, childOffset);
  }

  public isInbound(rootPath: Readonly<NodeAddress>, link: Link, childOffset?: number) {
    return this.hasNested(rootPath, link.prefix, link.matchInfo.outputs, childOffset) &&
      this.hasNonNested(rootPath, link.prefix, link.matchInfo.inputs, childOffset);
  }

  public isOutgoing(rootPath: Readonly<NodeAddress>, link: Link, childOffset?: number) {
    return this.hasNested(rootPath, link.prefix, link.matchInfo.inputs, childOffset) &&
      this.hasNonNested(rootPath, link.prefix, link.matchInfo.outputs, childOffset);
  }

  private hasNested(rootPath: Readonly<NodeAddress>, prefix: Readonly<NodeAddress>, minfos: Record<string, MatchedNodePaths>, childOffset?: number) {
    return Object.entries(minfos).some(
      ([, minfo]) => minfo.some((io) => BaseTree.isNodeChildOffseted(rootPath, [...prefix, ...io.path], childOffset)));
  }

  private hasNonNested(rootPath: Readonly<NodeAddress>, prefix: Readonly<NodeAddress>, minfos: Record<string, MatchedNodePaths>, childOffset?: number) {
    return Object.entries(minfos).some(
      ([, minfo]) => minfo.some((io) => !BaseTree.isNodeChildOffseted(rootPath, [...prefix, ...io.path], childOffset)));
  }

  private checkDisjoint(inbound: Link[], outgoing: Link[], mutationPath: NodePath, childOffset?: number) {
    const inboundSet = new Set(inbound.map((item) => item.uuid));
    const outgoingSet = new Set(outgoing.map((item) => item.uuid));
    for (const out of outgoingSet) {
      if (inboundSet.has(out)) {
        const link = inbound.find((l) => l.uuid === out);
        throw new Error(`Link ${link?.matchInfo.spec.id} uuid ${link?.uuid} mutationPath ${JSON.stringify(mutationPath)} childOffset ${childOffset} is both inbound and outgoing`);
      }
    }
  }

  private getLinkRunObs(linkUUID: string, mutationPath: NodePath, childOffset?: number) {
    const link = this.links.get(linkUUID)!;
    const obs$ = defer(() => {
      link.trigger(mutationPath, childOffset);
      return this.runningLinks$.pipe(filter((x) => x?.length === 0), mapTo(true as const));
    });
    return obs$;
  }
}
