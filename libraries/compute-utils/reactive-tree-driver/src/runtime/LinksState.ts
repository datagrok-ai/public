import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BaseTree, NodeAddress, NodePath} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {isFuncCallNode, isSequentialPipelineNode, isStaticPipelineNode} from './StateTreeNodes';
import {ActionSpec, isActionSpec, matchNodeLink} from './link-matching';
import {Action, Link} from './Link';
import {BehaviorSubject, merge, Subject} from 'rxjs';
import {takeUntil, map, scan, switchMap} from 'rxjs/operators';

export class DependenciesData {
  nodes: Set<string> = new Set();
  links: Set<string> = new Set();
}

export class LinksState {
  private closed$ = new Subject<true>();
  private linksUpdates = new Subject<true>();

  public links: Map<string, Link> = new Map();
  public actions: Map<string, Action> = new Map();
  public baseNodeActions: Map<string, Action[]> = new Map();
  public deps: Map<string, DependenciesData> = new Map();

  public runningLinks$ = new BehaviorSubject<undefined | string[]>(undefined);

  constructor() {
    this.linksUpdates.pipe(
      switchMap(() => this.getRunningLinks()),
      takeUntil(this.closed$),
    ).subscribe(this.runningLinks$);
  }

  public update(state: StateTree, mutationPath?: NodePath, childOffset?: number) {
    this.destroyLinks();
    const links = this.createAutoLinks(state);
    this.links = new Map(links.map((link) => [link.uuid, link] as const));
    const tree = this.createActionLinksTree(state);
    [this.actions, this.baseNodeActions] = tree;
    const deps = this.calculateDependencies(state, links);
    this.deps = deps;
    if (mutationPath) {
      const inbound = links.filter((link) => this.isInbound(mutationPath, link, childOffset));
      const outgoing = links.filter((link) => this.isOutgoing(mutationPath, link, childOffset));
      this.checkDisjoint(inbound, outgoing, mutationPath, childOffset);
      this.linksUpdates.next(true);
      return {regionData: {inbound, outgoing}};
    }
    this.linksUpdates.next(true);
    return {};
  }

  public createAutoLinks(state: StateTree) {
    const links = state.nodeTree.traverse(state.nodeTree.root, (acc, node, path) => {
      const item = node.getItem();
      if (isStaticPipelineNode(item) || isSequentialPipelineNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => new Link(path, minfo));
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }

  public createActionLinksTree(state: StateTree) {
    const links = state.nodeTree.traverse(state.nodeTree.root, (acc, node, path) => {
      const item = node.getItem();
      const {config} = item;
      const matchedLinks = (config.actions ?? [])
        .map((link) => matchNodeLink(node, link))
        .filter((x) => !!x)
        .flat();
      const links = matchedLinks.map((minfo) => {
        const spec = minfo.spec as ActionSpec;
        return new Action(path, minfo, spec.position, spec.friendlyName, spec.menuCategory);
      });
      return [...acc, ...links];
    }, [] as Action[]);
    const actionsMap = new Map<string, Action[]>;
    for (const link of links) {
      if (link.matchInfo.basePath) {
        const item = state.nodeTree.getItem(link.matchInfo.basePath);
        const actions = actionsMap.get(item.uuid) ?? [];
        actions.push(link);
        actionsMap.set(item.uuid, actions);
      }
    }
    const actions = new Map(links.map((link) => [link.uuid, link] as const));
    return [actions, actionsMap] as const;
  }

  // TODO: detection of cycles and double io deps, including self deps
  public calculateDependencies(state: StateTree, links: Link[]) {
    const deps: Map<string, DependenciesData> = new Map();
    for (const link of links) {
      for (const infosIn of Object.values(link.matchInfo.inputs)) {
        for (const infoIn of infosIn) {
          const nodeIn = state.nodeTree.getNode(infoIn.path);
          // pipeline memory state should be immutable and set in onInit
          if (!isFuncCallNode(nodeIn.getItem()))
            continue;
          for (const infosOut of Object.values(link.matchInfo.outputs)) {
            for (const infoOut of infosOut) {
              const nodeOut = state.nodeTree.getNode(infoOut.path);
              const depData = deps.get(nodeOut.getItem().uuid) ?? new DependenciesData();
              if (!BaseTree.isNodeAddressEq(infoIn.path, infoOut.path))
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

  public isLinkReady(state: StateTree, linkUUID: string) {
    const link = this.links.get(linkUUID);
    if (!link)
      return false;
    for (const infosIn of Object.values(link.matchInfo.inputs)) {
      for (const infoIn of infosIn) {
        const nodeIn = state.nodeTree.getNode(infoIn.path);
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

  public close() {
    this.closed$.next(true);
  }

  private getRunningLinks() {
    const obs = [...this.links.values(), ...[...this.baseNodeActions.values()].flat()].map(
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

  private isInbound(rootPath: Readonly<NodeAddress>, link: Link, childOffset?: number) {
    return Object.entries(link.matchInfo.outputs).some(
      ([, ios]) => ios.some((io) => BaseTree.isNodeChildOffseted(rootPath, io.path, childOffset)));
  }

  private isOutgoing(rootPath: Readonly<NodeAddress>, link: Link, childOffset?: number) {
    return Object.entries(link.matchInfo.inputs).some(
      ([, ios]) => ios.some((io) => BaseTree.isNodeChildOffseted(rootPath, io.path, childOffset)));
  }

  private checkDisjoint(inbound: Link[], outgoing: Link[], mutationPath: NodePath, childOffset?: number) {
    const inboundSet = new Set(inbound.map((item) => item.uuid));
    const outgoingSet = new Set(outgoing.map((item) => item.uuid));
    for (const out of outgoingSet) {
      if (inboundSet.has(out)) {
        const link = inbound.find((l) => l.uuid === out);
        throw new Error(`Link ${link?.matchInfo.spec.id} uuid ${link?.uuid} mutationPath ${JSON.stringify(mutationPath)} childOffset ${childOffset}`);
      }
    }
  }

  private destroyLinks() {
    for (const [, link] of this.links)
      link.destroy();
  }
}
