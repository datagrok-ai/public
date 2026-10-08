import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {NodePath, TreeNode} from '../data/BaseTree';
import {formatLinkIO, formatLinkSegment, LinkIOParsed} from '../config/LinkSpec';
import {formatNodePath} from '../utils';
import {explainLinkMatch, LinkMatchExplanation, MatchedNodePaths, matchNodeLink} from './link-matching';
import {LinksData} from './LinksState';
import {StateTree} from './StateTree';
import {flattenTags} from './StateTreeSerializer';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';

export type InspectedNodeRef = {
  uuid: string;
  configId: string;
  path: string;
};

export type InspectedLink = {
  id: string;
  uuid: string;
  type: string;
  isAction: boolean;
  node: InspectedNodeRef;
  basePath?: string;
  resolvedInputs: Record<string, string[]>;
  resolvedOutputs: Record<string, string[]>;
  // aliases that resolved to nothing: optional matches, or all targets dropped by $linked
  emptyAliases?: Record<string, 'optional' | 'linked'>;
  // targets dropped from $linked aliases because data links write them
  linkedTargets?: Record<string, string[]>;
  defaultRestrictions?: any;
  dataFrameMutations?: any;
  hasHandler?: boolean;
};

export type NotMatchedLink = {
  id: string;
  type: string;
  isAction: boolean;
  node: InspectedNodeRef;
  explanation: LinkMatchExplanation;
};

export type LinksInspection = {
  matched: InspectedLink[];
  notMatched: NotMatchedLink[];
};

export type InspectedIO = {name: string, type?: string, nullable?: boolean};

export type InspectedNode = InspectedNodeRef & {
  friendlyName?: string;
  type: string;
  isReadonly: boolean;
  inputs?: InspectedIO[];
  outputs?: InspectedIO[];
  callState?: any;
  validations?: any;
  consistency?: any;
  meta?: Record<string, any>;
  descriptions?: Record<string, any>;
  states?: Record<string, any>;
  pipelineValidations?: any;
};

export function inspectLinks(state: StateTree): LinksInspection {
  const {linksState, nodeTree} = state;
  const links = linksState.getLinksInfo();
  const matched = links.map((link) => inspectLink(state, link));
  const instanceKeys = new Set(links.map((link) => `${formatNodePath(link.prefix)}#${link.id}`));
  const notMatched = nodeTree.traverse(nodeTree.root, (acc, node, path) => {
    const {config} = node.getItem();
    const specs = [
      ...(config.links ?? []).map((spec) => [spec, false] as const),
      ...(config.actions ?? []).map((spec) => [spec, true] as const),
    ];
    for (const [spec, isAction] of specs) {
      if (instanceKeys.has(`${formatNodePath(path)}#${spec.id}`))
        continue;
      acc.push({
        id: spec.id, type: spec.type ?? 'data', isAction, node: nodeRef(node, path),
        explanation: explainLinkMatch(node, spec),
      });
    }
    return acc;
  }, [] as NotMatchedLink[]);
  return {matched, notMatched};
}

function inspectLink(state: StateTree, link: LinksData): InspectedLink {
  const {matchInfo, prefix, isAction} = link;
  const spec = matchInfo.spec;
  const node = state.nodeTree.getNode(prefix);
  const res: InspectedLink = {
    id: spec.id, uuid: link.uuid, type: spec.type ?? 'data', isAction, node: nodeRef(node, prefix),
    resolvedInputs: formatMatched(matchInfo.inputs),
    resolvedOutputs: formatMatched(matchInfo.outputs),
  };
  if (matchInfo.basePath)
    res.basePath = formatNodePath(matchInfo.basePath);
  if ('defaultRestrictions' in spec && spec.defaultRestrictions)
    res.defaultRestrictions = spec.defaultRestrictions;
  if ('dataFrameMutations' in spec && spec.dataFrameMutations)
    res.dataFrameMutations = spec.dataFrameMutations;
  if ('handler' in spec && spec.handler)
    res.hasHandler = true;

  const emptyAliases: Record<string, 'optional' | 'linked'> = {};
  // $linked pruning keeps no record, so a fresh match shows what was dropped
  const fresh = spec.to.some((io) => io.unlinked) ? matchNodeLink(node, spec, matchInfo.basePath)?.[0] : undefined;
  const freshOutputs = fresh ? formatMatched(fresh.outputs) : {};
  const linkedTargets: Record<string, string[]> = {};
  for (const io of [...spec.from, ...spec.to]) {
    const current = res.resolvedInputs[io.name] ?? res.resolvedOutputs[io.name];
    const dropped = io.unlinked ? (freshOutputs[io.name] ?? []).filter((p) => !current?.includes(p)) : [];
    if (dropped.length)
      linkedTargets[io.name] = dropped;
    if (!current)
      emptyAliases[io.name] = dropped.length ? 'linked' : 'optional';
  }
  if (Object.keys(emptyAliases).length)
    res.emptyAliases = emptyAliases;
  if (Object.keys(linkedTargets).length)
    res.linkedTargets = linkedTargets;
  return res;
}

export function inspectNode(state: StateTree, uuid: string): InspectedNode | undefined {
  const [node, path] = state.nodeTree.find((item) => item.uuid === uuid) ?? [];
  if (!node || !path)
    return undefined;
  const item = node.getItem();
  const res: InspectedNode = {
    ...nodeRef(node, path), friendlyName: item.config.friendlyName, type: item.nodeType, isReadonly: item.isReadonly,
  };
  const descriptions: Record<string, any> = {};
  for (const name of item.nodeDescription.getStateNames()) {
    const val = item.nodeDescription.getState(name);
    if (val === undefined)
      continue;
    descriptions[name] = name === 'tags' ? flattenTags(val) : val;
  }
  if (Object.keys(descriptions).length)
    res.descriptions = descriptions;

  if (isFuncCallNode(item)) {
    Object.assign(res, splitIO(item.config.io ?? []));
    res.callState = item.funcCallState$.value;
    res.validations = item.validationInfo$.value;
    res.consistency = item.consistencyInfo$.value;
    res.meta = Object.fromEntries(Object.entries(item.metaInfo$.value).map(([k, meta$]) => [k, meta$?.value]));
  } else {
    const store = item.getStateStore();
    const names = store.getStateNames();
    if (names.length)
      res.states = Object.fromEntries(names.map((name) => [name, store.getState(name)]));
    res.pipelineValidations = item.pipelineValidations$.value;
  }
  return res;
}

/** Plain JSON for the Inspector: DG objects, rx subjects, collections and handlers become readable values. */
export function toInspectorJSON(value: any): any {
  return value ? JSON.parse(JSON.stringify(value, inspectorReplacer)) : {};
}

/** Config as plain JSON, with each step's io split into inputs and outputs like in {@link inspectNode}. */
export function inspectConfig(config: any): any {
  return splitConfigIO(toInspectorJSON(config));
}

function nodeRef(node: TreeNode<StateTreeNode>, path: Readonly<NodePath>): InspectedNodeRef {
  const item = node.getItem();
  return {uuid: item.uuid, configId: item.config.id, path: formatNodePath(path)};
}

function formatMatched(paths: Record<string, MatchedNodePaths>) {
  return Object.fromEntries(Object.entries(paths).map(([name, matched]) =>
    [name, matched.map((m) => formatNodePath(m.path, m.ioName))]));
}

function splitIO(io: {id: string, type?: string, direction?: string, nullable?: boolean}[]) {
  const inputs: InspectedIO[] = [];
  const outputs: InspectedIO[] = [];
  for (const item of io) {
    const decl: InspectedIO = {name: item.id, type: item.type};
    if (item.nullable)
      decl.nullable = true;
    (item.direction === 'output' ? outputs : inputs).push(decl);
  }
  return {...(inputs.length ? {inputs} : {}), ...(outputs.length ? {outputs} : {})};
}

function splitConfigIO(data: any): any {
  if (data == null || typeof data !== 'object')
    return data;
  if (Array.isArray(data))
    return data.map(splitConfigIO);
  const res: Record<string, any> = {};
  for (const [k, v] of Object.entries(data)) {
    if (k === 'io' && Array.isArray(v) && v[0]?.id && v[0]?.direction)
      Object.assign(res, splitIO(v));
    else
      res[k] = splitConfigIO(v);
  }
  return res;
}

function isLinkIOParsed(v: any): v is LinkIOParsed {
  return v && typeof v.name === 'string' && Array.isArray(v.segments) && !v.type;
}

function isLinkSegment(v: any) {
  return v && ((v.type === 'selector' && Array.isArray(v.ids)) || (v.type === 'tag' && Array.isArray(v.tags)));
}

function inspectorReplacer(_key: string, value: any): any {
  if (isLinkIOParsed(value))
    return formatLinkIO(value);
  if (isLinkSegment(value))
    return formatLinkSegment(value);
  if (value instanceof DG.FuncCall)
    return {'#': 'FuncCall', 'id': value.id, 'func': value.func?.nqName};
  if (value instanceof DG.Func)
    return `#Func(${value.nqName})`;
  if (value instanceof DG.DataFrame)
    return `#DataFrame(${value.rowCount} rows, ${value.columns.length} cols)`;
  if (value instanceof Map)
    return Object.fromEntries(value);
  if (value instanceof Set)
    return [...value];
  if (value instanceof BehaviorSubject)
    return value.value;
  if (typeof value === 'function')
    return '#Handler';
  return value;
}
