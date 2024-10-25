import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BaseTree, NodeAddress, NodePath} from '../data/BaseTree';
import {StateTree} from './StateTree';
import {isFuncCallNode, isSequentialPipelineNode, isStaticPipelineNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, isActionSpec, LinkSpec, MatchedNodePaths, MatchInfo, matchNodeLink} from './link-matching';
import {Action, Link} from './Link';
import {BehaviorSubject, concat, merge, Subject, of, Observable, defer, combineLatest} from 'rxjs';
import {takeUntil, map, scan, switchMap, filter, mapTo, toArray, take, tap, debounceTime, delay, concatMap} from 'rxjs/operators';
import {parseLinkIO} from '../config/LinkSpec';
import {makeValidationResult} from '../utils';
import {DriverLogger} from '../data/Logger';

export interface NestedMutationData {
  mutationRootPath: NodeAddress,
  addIdx?: number,
  removeIdx?: number,
}

export interface LinksData {
  prefix: NodePath;
  isAction: boolean;
  matchInfo: MatchInfo;
}

class DependenciesData {
  nodes: Set<string> = new Set();
  links: Set<string> = new Set();
}

interface IoDep {
  data?: string;
  meta?: string;
  validation?: string[];
}

type IoDeps = Record<string, IoDep>;

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

  constructor(
    private defaultValidators: boolean = false,
    private logger?: DriverLogger,
  ) {
    this.linksUpdates.pipe(
      switchMap(() => this.getRunningLinks()),
      takeUntil(this.closed$),
    ).subscribe(this.runningLinks$);
  }

  public update(state: BaseTree<StateTreeNode>, nestedMutationData?: NestedMutationData) {
    this.destroyLinks();
    const links = this.createLinks(state);
    if (this.defaultValidators) {
      const validators = this.createDefaultValidators(state);
      links.push(...validators);
    }
    this.links = new Map(links.map((link) => [link.uuid, link] as const));
    [this.actions, this.nodesActions] = this.createActions(state);
    this.stepsDependencies = this.calculateStepsDependencies(state, links);
    this.ioDependencies = this.calculateIoDependencies(state, links);

    if (nestedMutationData) {
      const {mutationRootPath, addIdx, removeIdx} = nestedMutationData;
      const bound = this.getLowerBound(addIdx, removeIdx);
      const inbound = links.filter((link) => this.isDataLink(link) && this.isInbound(mutationRootPath, link, addIdx, removeIdx));
      const outgoing = links.filter((link) => this.isDataLink(link) && this.isOutgoing(mutationRootPath, link, addIdx));
      const affectedMeta = links
        .filter((link) => link.isMeta && this.isAffected(mutationRootPath, link, addIdx, removeIdx));
      const validatorLinks = links.filter((link) => link.isValidator);
      const inboundMap = new Map(inbound.map((link) => [link.uuid, link]));
      const outgoingMap = new Map(outgoing.map((link) => [link.uuid, link]));
      const validatorsMap = new Map(validatorLinks.map((link) => [link.uuid, link]));
      const affectedMap = new Map(affectedMeta.map((link) => [link.uuid, link]));
      this.linksUpdates.next(true);
      return concat(
        of(this.wireLinks(state)),
        this.runNewInits(state),
        this.runLinks(state, inboundMap, mutationRootPath, bound),
        this.runLinks(state, outgoingMap),
        this.runLinks(state, affectedMap),
        this.runLinks(state, validatorsMap),
      ).pipe(toArray(), mapTo(undefined));
    } else {
      this.linksUpdates.next(true);
      const metaMap = new Map(links.filter((link) => !this.isDataLink(link)).map((link) => [link.uuid, link]));
      return concat(
        of(this.wireLinks(state)),
        this.runNewInits(state),
        (this.defaultValidators || this.forceInitialMetaRun) ? of(null).pipe(delay(0), concatMap(() => this.runLinks(state, metaMap))) : of(null),
      ).pipe(toArray(), mapTo(undefined));
    }
  }

  public createLinks(state: BaseTree<StateTreeNode>) {
    const links = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (isStaticPipelineNode(item) || isSequentialPipelineNode(item)) {
        const {config} = item;
        const matchedLinks = (config.links ?? [])
          .map((link) => matchNodeLink(node, link))
          .filter((x) => !!x)
          .flat();
        const links = matchedLinks.map((minfo) => {
          const link = new Link(path, minfo, undefined, this.logger);
          return link;
        });
        return [...acc, ...links];
      }
      return acc;
    }, [] as Link[]);
    return links;
  }

  public createActions(state: BaseTree<StateTreeNode>) {
    const actionEntries = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      const {config} = item;
      const matchedLinks = (config.actions ?? [])
        .map((link) => matchNodeLink(node, link))
        .filter((x) => !!x)
        .flat();
      const links = matchedLinks.map((minfo) => {
        const spec = minfo.spec as ActionSpec;
        const action = new Action(path, minfo, spec.position, !!spec.isPipeline, spec.friendlyName, spec.menuCategory);
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
    const actionsMap = new Map(actionEntries.map(([, action]) => [action.uuid, action]));
    return [actionsMap, nodeActions] as const;
  }

  public createDefaultValidators(state: BaseTree<StateTreeNode>) {
    const defaultValidators = state.traverse(state.root, (acc, node, path) => {
      const item = node.getItem();
      if (!isFuncCallNode(item))
        return acc;
      const validators = item.config.io?.map((io) => {
        if (io.nullable || io.direction === 'output')
          return;
        const spec: LinkSpec = {
          id: uuidv4(),
          from: [parseLinkIO(`in:${io.id}`, io.direction)],
          to: [parseLinkIO(`out:${io.id}`, io.direction)],
          isValidator: true,
          handler({controller}) {
            const val = controller.getFirst('in');
            if (val == null)
              controller.setValidation('out', makeValidationResult({errors: ['Missing value']}));
            else
              controller.setValidation('out', undefined);
          },
        };
        // don't need to really match anything, just set to current node
        const minfo: MatchInfo = {
          spec,
          inputs: {
            'in': [{
              path: [],
              ioName: io.id,
            }],
          },
          outputs: {
            'out': [{
              path: [],
              ioName: io.id,
            }],
          },
          actions: {},
        };
        return new Link(path, minfo, 0);
      }).filter((x) => !!x);
      return [...acc, ...(validators ?? [])];
    }, [] as Link[]);
    return defaultValidators;
  }

  // TODO: cycles detection
  public calculateStepsDependencies(state: BaseTree<StateTreeNode>, links: Link[]) {
    const deps: Map<string, DependenciesData> = new Map();
    for (const link of links) {
      for (const infosIn of Object.values(link.matchInfo.inputs)) {
        for (const infoIn of infosIn) {
          const inPathFull = [...link.prefix, ...infoIn.path];
          const nodeIn = state.getNode(inPathFull);
          // pipeline memory states should be immutable and set in
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

  public calculateIoDependencies(state: BaseTree<StateTreeNode>, links: Link[]) {
    const deps = new Map<string, IoDeps>();
    for (const link of links) {
      const linkId = link.uuid;
      for (const infosOut of Object.values(link.matchInfo.outputs)) {
        for (const infoOut of infosOut) {
          const stepPath = [...link.prefix, ...infoOut.path];
          const stepNode = state.getNode(stepPath);
          const ioName = infoOut.ioName!;
          const depsData = deps.get(stepNode.getItem().uuid) ?? {};
          const depType = link.isValidator ? 'validation' : (link.isMeta ? 'meta' : 'data');
          if (depsData[ioName] == null)
            depsData[ioName] = {};
          if (depType === 'meta' || depType === 'data') {
            if (depsData[depType])
              grok.shell.warning(`Duplicate deps path ${JSON.stringify(stepPath)} io ${ioName} $`);

            depsData[ioName][depType] = linkId;
          } else {
            const currentDeps = depsData[ioName]['validation'] ?? [];
            currentDeps.push(linkId);
            depsData[ioName]['validation'] = currentDeps;
          }
          deps.set(stepNode.getItem().uuid, depsData);
        }
      }
    }
    return deps;
  }

  public runLinks(state: BaseTree<StateTreeNode>, links: Map<string, Link>, mutationRootPath?: NodeAddress, childOffset?: number) {
    const scheduledLinks = new Set<string>();
    const obs = state.traverse(state.root, (acc, node) => {
      const item = node.getItem();
      const deps = this.stepsDependencies.get(item.uuid);
      const linkRuns = [...(deps?.links ?? [])]
        .filter((linkUUID) => {
          const link = links.get(linkUUID);
          return link && !scheduledLinks.has(linkUUID);
        })
        .map((linkUUID) => {
          scheduledLinks.add(linkUUID);
          return combineLatest([
            this.runningLinks$,
            this.getLinkRunObs(linkUUID, mutationRootPath, childOffset),
          ]).pipe(
            debounceTime(0),
            filter(([running]) => !running || running.length === 0),
            take(1),
            mapTo(undefined),
          );
        });
      return [...acc, ...linkRuns];
    }, [] as Observable<undefined>[]);
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
    for (const [, action] of this.actions)
      action.wire(state, this);
  }

  public destroyLinks() {
    for (const [, link] of this.links)
      link.destroy();
  }

  public close() {
    this.closed$.next(true);
  }

  public getLinksInfo(): LinksData[] {
    const links = [...this.links.values()].map((l) => ({ prefix: l.prefix, isAction: false, matchInfo: l.matchInfo}));
    const actions = [...this.actions.values()].map((l) => ({ prefix: l.prefix, isAction: true, matchInfo: l.matchInfo}));
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
    return !isActionSpec(link.matchInfo.spec) && !link.matchInfo.spec.isValidator && !link.matchInfo.spec.isMeta;
  }

  public isAffected(rootPath: Readonly<NodeAddress>, link: Link, addIdx?: number, removeIdx?: number) {
    const bound = this.getLowerBound(addIdx, removeIdx);
    return this.hasNested(rootPath, link.prefix, link.matchInfo.inputs, bound) ||
      this.hasNested(rootPath, link.prefix, link.matchInfo.outputs, bound);
  }

  public isInbound(rootPath: Readonly<NodeAddress>, link: Link, addIdx?: number, removeIdx?: number) {
    if (addIdx != null) {
      const addedNodePath = [...rootPath, {idx: addIdx}];
      const isNodeInbound = (this.hasNested(addedNodePath, link.prefix, link.matchInfo.outputs) &&
        this.hasNonNested(addedNodePath, link.prefix, link.matchInfo.inputs));
      if (isNodeInbound)
        return isNodeInbound;
    }

    const bound = this.getLowerBound(addIdx, removeIdx);
    return (this.hasNested(rootPath, link.prefix, link.matchInfo.outputs, bound) &&
      this.hasNonNested(rootPath, link.prefix, link.matchInfo.inputs, bound));
  }

  public isOutgoing(rootPath: Readonly<NodeAddress>, link: Link, addIdx?: number) {
    const addedNodePath = addIdx == null ? rootPath : [...rootPath, {idx: addIdx}];
    return this.hasNested(addedNodePath, link.prefix, link.matchInfo.inputs) &&
      this.hasNonNested(addedNodePath, link.prefix, link.matchInfo.outputs);
  }

  private getLowerBound(addIdx?: number, removeIdx?: number) {
    if (addIdx == null && removeIdx == null)
      return undefined;
    if (addIdx != null && removeIdx != null)
      return Math.min(addIdx, removeIdx);
    if (addIdx != null)
      return addIdx;
    if (removeIdx != null)
      return removeIdx;
  }

  private hasNested(rootPath: Readonly<NodeAddress>, prefix: Readonly<NodeAddress>, minfos: Record<string, MatchedNodePaths>, childOffset?: number) {
    return Object.entries(minfos).some(
      ([, minfo]) => minfo.some((io) => BaseTree.isNodeChildOffseted(rootPath, [...prefix, ...io.path], childOffset)));
  }

  private hasNonNested(rootPath: Readonly<NodeAddress>, prefix: Readonly<NodeAddress>, minfos: Record<string, MatchedNodePaths>, childOffset?: number) {
    return Object.entries(minfos).some(
      ([, minfo]) => minfo.some((io) => !BaseTree.isNodeChildOffseted(rootPath, [...prefix, ...io.path], childOffset)));
  }

  private getLinkRunObs(linkUUID: string, mutationPath?: NodeAddress, childOffset?: number) {
    const link = this.links.get(linkUUID)!;
    const obs$ = defer(() => {
      link.trigger(mutationPath, childOffset);
      return this.runningLinks$.pipe(filter((x) => x?.length === 0), mapTo(true as const));
    });
    return obs$;
  }
}
