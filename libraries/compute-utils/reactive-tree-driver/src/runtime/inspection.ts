import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {NodePath, TreeNode} from '../data/BaseTree';
import {AnnotationLinkKind} from '../data/common-types';
import {getOriginalConfig, PipelineConfigurationProcessed} from '../config/config-processing-utils';
import {formatNodePath} from '../utils';
import {MatchedNodePaths, matchNodeLink} from './link-matching';
import {LinksData} from './LinksState';
import {StateTree} from './StateTree';
import {StateTreeNode} from './StateTreeNodes';

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
  visible?: boolean;
  annotation?: AnnotationLinkKind;
  // the link, rule or check as written in the config
  original?: any;
};

export type InspectedIO = {name: string, type?: string, nullable?: boolean};

export function inspectLinks(state: StateTree): InspectedLink[] {
  return state.linksState.getLinksInfo(true).map((link) => inspectLink(state, link));
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
  if (spec.annotation)
    res.annotation = spec.annotation;
  if (isAction)
    res.visible = state.linksState.actionsVisibility.get(link.uuid) ?? true;
  const original = getOriginalConfig(spec);
  if (original)
    res.original = toInspectorJSON(original);

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

/** Plain JSON for the Inspector: DG objects, rx subjects, collections and handlers become readable values. */
export function toInspectorJSON(value: any): any {
  return value ? JSON.parse(JSON.stringify(value, inspectorReplacer)) : {};
}

/** The config as written, with nested pipelines loaded and each step's io resolved into inputs and outputs. */
export function inspectConfig(config: PipelineConfigurationProcessed): any {
  return toInspectorJSON(originalConfig(config));
}

function originalConfig(processed: any): any {
  const res: Record<string, any> = {...(getOriginalConfig(processed) ?? processed)};
  for (const key of ['steps', 'stepTypes']) {
    if (res[key])
      res[key] = processed[key].map(originalConfig);
  }
  if (processed.io) {
    delete res.io;
    Object.assign(res, splitIO(processed.io));
  }
  return res;
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

function inspectorReplacer(_key: string, value: any): any {
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
