import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {v4 as uuidv4} from 'uuid';
import {BaseTree, NodePath} from '../data/BaseTree';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';
import {ActionSpec, LinkSpec, MatchInfo, matchNodeLink} from './link-matching';
import {Action, Link} from './Link';
import {BehaviorSubject, concat, merge, Subject, of, Observable, defer, combineLatest} from 'rxjs';
import {takeUntil, map, scan, switchMap, filter, mapTo, toArray, take, tap, debounceTime, delay, concatMap} from 'rxjs/operators';
import {parseLinkIO} from '../config/LinkSpec';
import {makeValidationResult} from '../utils';
import {DriverLogger} from '../data/Logger';
import {getLinksDiff} from './links-diff';
import {ViewAction} from '../config/PipelineInstance';

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

  public update(state: BaseTree<StateTreeNode>, isMutation: boolean) {
    this.disableLinks();

    const oldLinks = [...this.links.values()];
    const [links, addedLinks] = this.createLinks(state, oldLinks);
    this.links = new Map(links.map((link) => [link.uuid, link] as const));

    [this.actions, this.nodesActions] = this.createActions(state);
    this.stepsDependencies = this.calculateStepsDependencies(state, links);
    this.ioDependencies = this.calculateIoDependencies(state, links);

    if (isMutation) {
      const newData = addedLinks.filter((link) => this.isDataLink(link));
      const newMeta = addedLinks.filter((link) => !this.isDataLink(link));
      this.linksUpdates.next(true);
      return concat(
        of(this.wireLinks(state)),
        this.runNewInits(state),
        this.runLinks(state, this.toLinksMap(newData)),
        this.runLinks(state, this.toLinksMap(newMeta)),
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

  public getNodeActionsData(uuid: string): ViewAction[] | undefined {
    const actions = this.nodesActions.get(uuid);
    if (actions) {
      return actions.map(({uuid, spec: {position, description, menuCategory, friendlyName, icon, confirmationMessage}}) =>
        ({uuid, position, description, menuCategory, friendlyName, icon, confirmationMessage}));
    }
    return actions;
  }

  public createLinks(state: BaseTree<StateTreeNode>, oldLinks: Link[]) {
    const addedLinks: Link[] = [];
    const newLinks = this.createStateLinks(state);
    if (this.defaultValidators) {
      const validators = this.createDefaultValidators(state);
      newLinks.push(...validators);
    }
    const {toRemove, toAdd} = getLinksDiff(oldLinks, newLinks);
    const mergedLinks: Link[] = [];
    for (const oldLink of oldLinks) {
      if (!toRemove.has(oldLink.uuid))
        mergedLinks.push(oldLink);
      else {
        if (this.logger && !oldLink.matchInfo.isDefaultValidator)
          this.logger.logLink('linkRemoved', {linkUUID: oldLink.uuid, prefix: oldLink.prefix, id: oldLink.matchInfo.spec.id});
        oldLink.destroy();
      }
    }
    for (const newLink of newLinks) {
      if (toAdd.has(newLink.uuid)) {
        if (this.logger && !newLink.matchInfo.isDefaultValidator)
          this.logger.logLink('linkAdded', {linkUUID: newLink.uuid, prefix: newLink.prefix, id: newLink.matchInfo.spec.id});
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
        const action = new Action(path, minfo, spec, this.logger);
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
          type: 'validator',
          handler({controller}) {
            const val = controller.getFirst('in');
            if (val == null || val === '')
              controller.setValidation('out', makeValidationResult({errors: ['Missing value']}));
            else
              controller.setValidation('out', undefined);
          },
        };
        // don't really need to match anything, just set to the current node
        const minfo: MatchInfo = {
          spec,
          inputs: {
            'in': [{
              path: [],
              ioName: io.id,
            }],
          },
          inputsUUID: new Map(),
          outputs: {
            'out': [{
              path: [],
              ioName: io.id,
            }],
          },
          actions: {},
          outputsUUID: new Map(),
          isDefaultValidator: true,
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
          const node = state.getNode(stepPath);
          const ioName = infoOut.ioName!;
          const depsData = deps.get(node.getItem().uuid) ?? {};
          const depType = link.matchInfo.spec.type ?? 'data';
          if (depsData[ioName] == null)
            depsData[ioName] = {};
          if (depType === 'meta' || depType === 'data') {
            if (depsData[depType])
              grok.shell.warning(`Duplicate deps path ${JSON.stringify(stepPath)} io ${ioName} $`);

            depsData[ioName][depType] = linkId;
          }
          deps.set(node.getItem().uuid, depsData);
        }
      }
    }
    return deps;
  }

  public runLinks(state: BaseTree<StateTreeNode>, links: Map<string, Link>) {
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
            this.getLinkRunObs(linkUUID),
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

  public disableLinks() {
    for (const [, link] of this.links)
      link.setInactive();
  }

  public close() {
    this.closed$.next(true);
  }

  public getLinksInfo(): LinksData[] {
    const links = [...this.links.values()].filter((l) => !l.matchInfo.isDefaultValidator).map((l) => ({prefix: l.prefix, isAction: false, matchInfo: l.matchInfo}));
    const actions = [...this.actions.values()].map((l) => ({prefix: l.prefix, isAction: true, matchInfo: l.matchInfo}));
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
