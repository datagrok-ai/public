import {isNonRefSelector, LinkNonRefSelectors, LinkIOParsed, LinkRefSelectors, refSelectorAdjacent, refSelectorAll, refSelectorDirection, refSelectorFindOne, TagRefSelectors} from '../config/LinkSpec';
import {DataActionConfiguraion, FuncCallActionConfiguration, PipelineLinkConfiguration, PipelineMutationConfiguration} from '../config/PipelineConfiguration';
import {BaseTree, NodePath, NodePathSegment, TreeNode} from '../data/BaseTree';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {indexFromEnd, pathToUUID} from '../utils';
import {StateTree} from './StateTree';
import {StateTreeNode} from './StateTreeNodes';

export type LinkSpec = PipelineLinkConfiguration<LinkIOParsed[]>;
export type ActionSpec = DataActionConfiguraion<LinkIOParsed[]> | PipelineMutationConfiguration<LinkIOParsed[]> | FuncCallActionConfiguration<LinkIOParsed[]>;

export type MatchedIO = {
  path: Readonly<NodePath>;
  ioName?: string;
}

export type MatchedNodePaths = Readonly<Array<MatchedIO>>;

export type MatchInfo = {
  spec: LinkSpec | ActionSpec;
  basePath?: Readonly<NodePath>;
  basePathUUID?: string;
  actions: Record<string, MatchedNodePaths>;
  inputs: Record<string, MatchedNodePaths>;
  outputs: Record<string, MatchedNodePaths>;
  // needed to diff order beetween multiple steps of the same kind
  inputsUUID: Map<string, Array<string>>;
  outputsUUID: Map<string, Array<string>>;
  isDefaultValidator?: boolean;
}

type NodeTraverseState = {
  pnode: TreeNode<StateTreeNode>,
  isLastSegment: boolean,
};

export function matchLink(state: StateTree, address: NodePath, spec: LinkSpec): MatchInfo[] | undefined {
  const [rnode] = state.nodeTree.find(((_node, path) => BaseTree.isNodeAddressEq(path, address))) ?? [];
  if (rnode == null)
    return;
  return matchNodeLink(rnode, spec);
}

export function matchNodeLink(rnode: TreeNode<StateTreeNode>, spec: LinkSpec | ActionSpec, basePath?: Readonly<NodePath>) {
  return matchLinkAt(rnode, spec, basePath, false).matched;
}

/** Runs the matching of one link at one node without stopping at the first empty alias. */
export function explainLinkMatch(rnode: TreeNode<StateTreeNode>, spec: LinkSpec | ActionSpec): LinkMatchExplanation {
  return matchLinkAt(rnode, spec, undefined, true).explanation!;
}

function matchLinkAt(
  rnode: TreeNode<StateTreeNode>,
  spec: LinkSpec | ActionSpec,
  basePath: Readonly<NodePath> | undefined,
  explain: boolean,
): {matched?: MatchInfo[], explanation?: LinkMatchExplanation} {
  const basePaths = basePath ? [{path: basePath}] : (spec.base?.length ? expandLinkBase(rnode, spec.base[0]) : undefined);
  const baseName = spec.base?.length ? spec.base[0].name : undefined;
  const explanation: LinkMatchExplanation | undefined = explain ? {aliases: []} : undefined;
  for (const io of spec.not ?? []) {
    const nodes = matchLinkIO(rnode, {}, io, true, false).length;
    if (!explanation) {
      if (nodes)
        return {};
      continue;
    }
    explanation.aliases.push({name: io.name, kind: 'not', optional: false, nodes});
    if (nodes && !explanation.blockedBy)
      explanation.blockedBy = io.name;
  }
  if (explanation && basePaths) {
    explanation.aliases.push({name: baseName!, kind: 'base', optional: false, nodes: basePaths.length});
    if (!basePaths.length) {
      explanation.failedAlias = baseName;
      return {explanation};
    }
  }
  const bases = basePaths == null ? [undefined] : (explanation ? basePaths.slice(0, 1) : basePaths);
  const instances = bases.map((base) => matchLinkInstance(rnode, spec, explain, base, baseName));
  if (explanation) {
    explanation.aliases.push(...instances[0].aliases!);
    explanation.failedAlias = instances[0].failedAlias;
  }
  const matched = instances.map((instance) => instance.matchInfo).filter((x) => !!x);
  return {matched: matched.length ? matched : undefined, explanation};
}

function matchLinkInstance(
  rnode: TreeNode<StateTreeNode>,
  spec: LinkSpec | ActionSpec,
  explain: boolean,
  base?: Readonly<MatchedIO>,
  baseName?: string,
): {matchInfo?: MatchInfo, aliases?: AliasMatch[], failedAlias?: string} {
  const actions: Record<string, MatchedNodePaths> = {};
  const matchInfo: MatchInfo = {
    spec,
    basePath: base?.path,
    actions,
    inputs: {},
    outputs: {},
    inputsUUID: new Map(),
    outputsUUID: new Map(),
  };
  const currentIO: Record<string, MatchedNodePaths> = {};
  if (base && baseName)
    currentIO[baseName] = [base];
  const aliases: AliasMatch[] | undefined = explain ? [] : undefined;
  let failedAlias: string | undefined;

  for (const action of spec.actions ?? []) {
    const parsed = matchLinkIO(rnode, currentIO, action, true, false);
    actions[action.name] = parsed;
    aliases?.push({name: action.name, kind: 'action', optional: false, nodes: parsed.length});
  }

  const ioData = [...spec.from.map((item) => ['inputs', item] as const), ...spec.to.map((item) => ['outputs', item] as const)];
  for (const [kind, io] of ioData) {
    const [skipIO, useDescriptionsStore] = ioMatchFlags(spec, kind === 'outputs', io);
    const dependsOn = explain ? io.segments.find((segment) => segment.ref && !currentIO[segment.ref])?.ref : undefined;
    const paths = dependsOn ? [] : matchLinkIO(rnode, currentIO, io, skipIO, useDescriptionsStore);
    aliases?.push(explainAlias(rnode, currentIO, io, kind, skipIO, paths, dependsOn));
    if (paths.length == 0) {
      if (io.flags?.includes('optional'))
        continue;
      if (!explain)
        return {};
      failedAlias ??= io.name;
      continue;
    }
    if (currentIO[io.name] != null && !explain)
      throw new Error(`Duplicate io name ${io.name} in link ${spec.id}, node: ${rnode.getItem().config.id}`);
    currentIO[io.name] = paths;
    matchInfo[kind][io.name] = paths;
  }
  if (failedAlias)
    return {aliases, failedAlias};
  updateMatchInfoUUIDs(rnode, matchInfo);
  return {matchInfo, aliases};
}

export type AliasMatch = {
  name: string;
  kind: 'not' | 'base' | 'action' | 'input' | 'output';
  optional: boolean;
  // nodes reached by the path, and of those the ones having the io
  nodes: number;
  ios?: number;
  // the alias it refers to matched nothing, so it was not matched
  dependsOn?: string;
};

export type LinkMatchExplanation = {
  aliases: AliasMatch[];
  failedAlias?: string;
  blockedBy?: string;
};

function explainAlias(
  rnode: TreeNode<StateTreeNode>,
  currentIO: Record<string, MatchedNodePaths>,
  io: LinkIOParsed,
  kind: 'inputs' | 'outputs',
  skipIO: boolean,
  paths: MatchedNodePaths,
  dependsOn?: string,
): AliasMatch {
  const entry: AliasMatch = {
    name: io.name, kind: kind === 'inputs' ? 'input' : 'output', optional: !!io.flags?.includes('optional'), nodes: 0,
  };
  if (dependsOn)
    entry.dependsOn = dependsOn;
  else if (skipIO)
    entry.nodes = paths.length;
  else {
    entry.nodes = matchLinkIO(rnode, currentIO, {...io, segments: io.segments.slice(0, -1)}, true, false).length;
    entry.ios = paths.length;
  }
  return entry;
}

function ioMatchFlags(spec: LinkSpec | ActionSpec, isOutput: boolean, io: LinkIOParsed) {
  const skipIO = (isOutput && (spec.type === 'pipeline' || spec.type === 'pipelineValidator')) ||
    !!io.flags?.includes('call');
  const useDescriptionsStore = isOutput && (spec.type === 'nodemeta' || spec.type === 'selector');
  return [skipIO, useDescriptionsStore] as const;
}

function expandLinkBase(
  rnode: TreeNode<StateTreeNode>,
  baseLink: LinkIOParsed,
) {
  const traverse = buildTraverseD([] as Readonly<NodePath>, (state: NodeTraverseState, currentPath, currentSegment?: number) => {
    const segment = baseLink.segments[currentSegment!];
    const {pnode} = state;
    const isLastSegment = baseLink.segments.length - 1 === currentSegment;
    if (segment?.type === 'selector') {
      const {selector, ids} = segment;
      const nextNodes = matchNonRefSegment(pnode, ids, selector as LinkNonRefSelectors);
      return nextNodes.map(([idx, node]) => [{pnode: node.item, isLastSegment}, [...currentPath, {id: node.id, idx}], currentSegment!+1] as const);
    } else if (segment?.type === 'tag') {
      const {selector, tags} = segment;
      const nextNodes = matchNonRefTag(pnode, tags, selector as LinkNonRefSelectors);
      return nextNodes.map(({path, node}) => [{pnode: node, isLastSegment}, [...currentPath, ...path], currentSegment!+1] as const);
    } else
      return [] as const;
  }, 0);

  const initialData = {
    pnode: rnode,
    isLastSegment: baseLink.segments.length === 0,
  };

  const basePaths = traverse(initialData, (acc, {isLastSegment}, path) => {
    if (isLastSegment)
      return [...acc, {path}] as const;
    return acc;
  }, [] as MatchedNodePaths);

  return basePaths;
}

function getRefOrigin(
  rnode: TreeNode<StateTreeNode>,
  currentIO: Record<string, MatchedNodePaths>,
  parsedLink: LinkIOParsed,
  ref: string | undefined,
  selector: LinkRefSelectors,
  path: readonly NodePathSegment[],
) {
  if (ref == null)
    return;
  const io = currentIO[ref];
  if (io == null)
    throw new Error(`Node ${rnode.getItem().config.id} referenced unknown io ${ref} in ${parsedLink.name}`);
  const refOrigin = refSelectorDirection(selector) === 'before' ? io[0] : indexFromEnd(io)!;
  if (!BaseTree.isNodeChildOrEq(path, refOrigin.path))
    return;
  return refOrigin;
}

export function isActionVisible(
  rnode: TreeNode<StateTreeNode>,
  spec: ActionSpec,
): boolean {
  const currentIO: Record<string, MatchedNodePaths> = {};
  for (const io of spec.hideWhen ?? []) {
    const paths = matchLinkIO(rnode, currentIO, io, true, false);
    if (paths.length > 0) return false;
  }
  for (const io of spec.showWhen ?? []) {
    if (io.flags?.includes('optional')) continue;
    const paths = matchLinkIO(rnode, currentIO, io, true, false);
    if (paths.length === 0) return false;
  }
  return true;
}

function matchLinkIO(
  rnode: TreeNode<StateTreeNode>,
  currentIO: Record<string, MatchedNodePaths>,
  parsedLink: LinkIOParsed,
  skipIO: boolean,
  useDescriptionStore: boolean,
): MatchedNodePaths {
  const traverse = buildTraverseD([] as Readonly<NodePath>, (state: NodeTraverseState, currentPath, currentSegment?: number) => {
    const segment = parsedLink.segments[currentSegment!];
    const {pnode} = state;
    const isLastSegment = parsedLink.segments.length - (skipIO ? 1 : 2) === currentSegment;
    if (segment?.type === 'selector') {
      const {ref, selector, ids, stopIds} = segment;
      let nextNodes: ReturnType<typeof matchNonRefSegment> = [];
      if (isNonRefSelector(selector))
        nextNodes = matchNonRefSegment(pnode, ids, selector);
      else {
        let originIdx = undefined;
        const refOrigin = getRefOrigin(rnode, currentIO, parsedLink, ref, selector, currentPath);
        if (refOrigin)
          originIdx = refOrigin.path[currentSegment!].idx;
        nextNodes = matchRefSegment(pnode, ids, selector, originIdx, stopIds);
      }
      return nextNodes.map(([idx, node]) => [{pnode: node.item, isLastSegment}, [...currentPath, {id: node.id, idx}], currentSegment!+1] as const);
    } else if (segment?.type === 'tag') {
      const {ref, selector, tags} = segment;
      let nextNodes: ReturnType<typeof matchNonRefTag> = [];
      if (isNonRefSelector(selector))
        nextNodes = matchNonRefTag(pnode, tags, selector);
      else {
        const refOrigin = getRefOrigin(rnode, currentIO, parsedLink, ref, selector, currentPath);
        if (!refOrigin)
          return [];
        const refNode = rnode.getNode(refOrigin.path);
        if (!refNode)
          return [];
        nextNodes = matchRefTag(pnode, tags, selector, refNode.getItem().uuid);
      }
      return nextNodes.map(({path, node}) => [{pnode: node, isLastSegment}, [...currentPath, ...path], currentSegment!+1] as const);
    }
    return [];
  }, 0);

  const initialData = {
    pnode: rnode,
    isLastSegment: parsedLink.segments.length === (skipIO ? 0 : 1),
  };

  const paths = traverse(initialData, (acc, state, path) => {
    const {pnode: node, isLastSegment} = state;
    if (isLastSegment && !skipIO) {
      const ioSegment = indexFromEnd(parsedLink.segments)!;
      if (ioSegment.type === 'tag')
        throw new Error(`Link ${parsedLink.name}, path ${JSON.stringify(path)} is ending with tag instead of io selector`);
      const ioName = ioSegment.ids[0];
      const item = node.getItem();
      const names = useDescriptionStore ? item.nodeDescription.getStateNames() : item.getStateStore().getStateNames();
      const state = names.find((name) => name === ioName);
      if (state)
        return [...acc, {path, ioName}] as const;
      return acc;
    } else if (isLastSegment && skipIO) {
      const p = {path};
      return [...acc, p] as const;
    }
    return acc;
  }, [] as MatchedNodePaths);

  return paths;
}

function matchNonRefSegment(
  pnode: TreeNode<StateTreeNode>,
  ids: string[],
  selector: LinkNonRefSelectors,
) {
  const idsSet = new Set(ids);
  const matchingNodes = [...pnode.getChildren().entries()].filter(([, c]) => idsSet.has(c.id));
  if (matchingNodes.length === 0)
    return [];
  if (selector === 'all' || selector === 'expand')
    return matchingNodes;
  if (selector === 'first')
    return [matchingNodes[0]];
  if (selector === 'last')
    return [indexFromEnd(matchingNodes)!];
  throw new Error(`Unknown segement mode ${selector}`);
}

function matchNonRefTag(pnode: TreeNode<StateTreeNode>, tags: string[], selector: LinkNonRefSelectors) {
  const matchingNodes = pnode.traverse((acc, node, path) => {
    const item = node.getItem();
    if (includesAll(item.config.tags ?? [], tags))
      acc!.push({path, node});

    return acc;
  }, [] as {path: readonly NodePathSegment[], node: TreeNode<StateTreeNode>}[]);
  if (matchingNodes.length === 0)
    return [];
  if (selector === 'all' || selector === 'expand')
    return matchingNodes;
  if (selector === 'first')
    return [matchingNodes[0]];
  if (selector === 'last')
    return [indexFromEnd(matchingNodes)!];
  throw new Error(`Unknown tag mode ${selector}`);
}

function matchRefSegment(
  pnode: TreeNode<StateTreeNode>,
  ids: string[],
  selector: LinkRefSelectors,
  originIdx?: number,
  stopIds?: string[],
) {
  const idsSet = new Set(ids);
  const stopSet = new Set(stopIds);
  const allNodeEntries = [...pnode.getChildren().entries()];

  if (selector === 'same') {
    if (originIdx == null)
      return [];
    const target = allNodeEntries[originIdx];
    if (ids.length > 0 && !idsSet.has(target[1].id))
      return [];
    return [target];
  }

  const selDirection = refSelectorDirection(selector);
  const partitionedByOrigin = originIdx == null ?
    allNodeEntries :
    allNodeEntries.filter(([idx]) => selDirection === 'before' ? idx < originIdx : idx > originIdx);

  let matchingNodes = partitionedByOrigin;
  if (stopSet?.size) {
    if (selDirection === 'before') {
      const startIdx = partitionedByOrigin.findLastIndex(([, node]) => stopSet.has(node.id));
      if (startIdx >= 0)
        matchingNodes = partitionedByOrigin.slice(startIdx+1);
    } else {
      const stopIdx = partitionedByOrigin.findIndex(([, node]) => stopSet.has(node.id));
      if (stopIdx >= 0)
        matchingNodes = partitionedByOrigin.slice(0, stopIdx);
    }
  }
  matchingNodes = matchingNodes.filter(([, c]) => idsSet.has(c.id));

  if (matchingNodes.length === 0)
    return [];
  if (refSelectorAll(selector))
    return matchingNodes;

  if (refSelectorAdjacent(selector)) {
    if (selDirection === 'before') {
      const adj = indexFromEnd(matchingNodes)!;
      const idx = adj[0];
      return idx + 1 === originIdx ? [adj] : [];
    } else {
      const adj = matchingNodes[0];
      const idx = adj[0];
      return idx - 1 === originIdx ? [adj] : [];
    }
  }
  if (refSelectorFindOne(selector)) {
    const items = selDirection === 'before' ? indexFromEnd(matchingNodes) : matchingNodes[0];
    return items ? [items] : [];
  }
  throw new Error(`Unknown segement mode ${selector}`);
}

function matchRefTag(
  pnode: TreeNode<StateTreeNode>,
  tags: string[],
  selector: TagRefSelectors,
  refUuid: string,
) {
  type MatchData = {path: readonly NodePathSegment[], node: TreeNode<StateTreeNode>};
  type MatchAcc = {
    before: MatchData[],
    same?: MatchData,
    after: MatchData[],
    isRefVisited: boolean,
  };

  const matchingNodesData = pnode.traverse((acc, node, path) => {
    const item = node.getItem();
    const {uuid} = item;
    if (includesAll(item.config.tags ?? [], tags)) {
      if (uuid === refUuid)
        acc!.same = {path, node};
      else if (!acc.isRefVisited)
        acc.before.push({path, node});
      else
        acc.after.push({path, node});
    }
    if (!acc.isRefVisited && uuid === refUuid)
      acc.isRefVisited = true;
    return acc;
  }, {before: [], same: undefined, after: [], isRefVisited: false} as MatchAcc);

  if (selector === 'same')
    return matchingNodesData.same ? [matchingNodesData.same] : [];


  const selDirection = refSelectorDirection(selector);
  if (refSelectorFindOne(selector)) {
    const items = selDirection === 'before' ? indexFromEnd(matchingNodesData.before) : matchingNodesData.after[0];
    return items ? [items] : [];
  }

  if (refSelectorAll(selector)) {
    if (selDirection === 'before')
      return matchingNodesData.before ?? [];
    else
      return matchingNodesData.after ?? [];
  }

  throw new Error(`Unknown tag mode ${selector}`);
}

function includesAll(target: string[], toInclide: string[]) {
  const targetSet = new Set(target);
  for (const item of toInclide) {
    if (!targetSet.has(item))
      return false;
  }
  return true;
}

export function updateMatchInfoUUIDs(rnode: TreeNode<StateTreeNode>, matchInfo: MatchInfo) {
  if (matchInfo.basePath) {
    const uuids = matchedPathsToUUIDs(rnode, [{path: matchInfo.basePath}]);
    matchInfo.basePathUUID = [...uuids.values()][0];
  }
  for (const [name, input] of Object.entries(matchInfo.inputs)) {
    const s = matchedPathsToUUIDs(rnode, input);
    matchInfo.inputsUUID.set(name, s);
  }
  for (const [name, output] of Object.entries(matchInfo.outputs)) {
    const s = matchedPathsToUUIDs(rnode, output);
    matchInfo.outputsUUID.set(name, s);
  }
}

function matchedPathsToUUIDs(rnode: TreeNode<StateTreeNode>, paths: MatchedNodePaths): string[] {
  const res: string[] = [];
  for (const path of paths) {
    const uuids = pathToUUID(rnode, path.path);
    if (path.ioName)
      uuids.push(path.ioName);
    res.push(uuids.join('/'));
  }
  return res;
}
