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
  inputsUUID: Map<string, Set<string>>;
  outputs: Record<string, MatchedNodePaths>;
  outputsUUID: Map<string, Set<string>>;
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
  const basePaths = basePath ? [{path: basePath}] : (spec.base?.length ? expandLinkBase(rnode, spec.base[0]) : undefined);
  const baseName = spec.base?.length ? spec.base[0].name : undefined;
  if (spec.not) {
    const currentIO: Record<string, MatchedNodePaths> = {};
    for (const io of spec.not) {
      const paths = matchLinkIO(rnode, currentIO, io, true, false);
      if (paths?.length)
        return undefined;
    }
  }
  if (basePaths == null) {
    const matchInfo = matchLinkInstance(rnode, spec);
    return matchInfo ? [matchInfo] : undefined;
  } else {
    const instances = basePaths.map((basePath) => matchLinkInstance(rnode, spec, basePath, baseName)).filter((x) => !!x);
    return instances.length ? instances : undefined;
  }
}

function matchLinkInstance(
  rnode: TreeNode<StateTreeNode>,
  spec: LinkSpec | ActionSpec,
  base?: Readonly<MatchedIO>,
  baseName?: string,
): MatchInfo | undefined {
  const actions: Record<string, MatchedNodePaths> = {};
  const matchInfo: MatchInfo = {
    spec,
    basePath: base?.path,
    actions,
    inputs: {},
    inputsUUID: new Map(),
    outputs: {},
    outputsUUID: new Map(),
  };
  const currentIO: Record<string, MatchedNodePaths> = {};
  if (base && baseName)
    currentIO[baseName] = [base];

  for (const action of spec.actions ?? []) {
    const parsed = matchLinkIO(rnode, currentIO, action, true, false);
    actions[action.name] = parsed;
  }

  const ioData = [...spec.from.map((item) => ['inputs', item] as const), ...spec.to.map((item) => ['outputs', item] as const)];
  for (const [kind, io] of ioData) {
    const skipIO = (spec.type === 'pipeline' && kind === 'outputs') || (!!io.flags?.includes('call'));
    const useDescriptionsStore = (spec.type === 'selector' && kind === 'outputs');
    const paths = matchLinkIO(rnode, currentIO, io, skipIO, useDescriptionsStore);
    if (paths.length == 0) {
      if (io.flags?.includes('optional'))
        continue;
      else
        return;
    }
    if (currentIO[io.name] != null)
      throw new Error(`Duplicate io name ${io.name} in link ${rnode.getItem().config.id}`);
    currentIO[io.name] = paths;
    matchInfo[kind][io.name] = paths;
  }
  updateMatchInfoUUIDs(rnode, matchInfo);
  return matchInfo;
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
  ref: string,
  selector: LinkRefSelectors,
  path: readonly NodePathSegment[],
) {
  const io = currentIO[ref];
  if (io == null)
    throw new Error(`Node ${rnode.getItem().config.id} referenced unknown io ${ref} in ${parsedLink.name}`);
  const refOrigin = refSelectorDirection(selector) === 'before' ? io[0] : indexFromEnd(io)!;
  if (!BaseTree.isNodeChildOrEq(path, refOrigin.path))
    throw new Error(`Node ${rnode.getItem().config.id} reference path ${JSON.stringify(refOrigin.path)} is different from current ${JSON.stringify(path)}`);
  return refOrigin;
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
        if (ref) {
          const refOrigin = getRefOrigin(rnode, currentIO, parsedLink, ref, selector, currentPath);
          originIdx = refOrigin.path[currentSegment!].idx;
        }
        nextNodes = matchRefSegment(pnode, ids, selector, originIdx, stopIds);
      }
      return nextNodes.map(([idx, node]) => [{pnode: node.item, isLastSegment}, [...currentPath, {id: node.id, idx}], currentSegment!+1] as const);
    } else if (segment?.type === 'tag') {
      const {ref, selector, tags} = segment;
      let nextNodes: ReturnType<typeof matchNonRefTag> = [];
      if (isNonRefSelector(selector))
        nextNodes = matchNonRefTag(pnode, tags, selector);
      else {
        if (!ref)
          return [];
        const refOrigin = getRefOrigin(rnode, currentIO, parsedLink, ref, selector, currentPath);
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

function matchedPathsToUUIDs(rnode: TreeNode<StateTreeNode>, paths: MatchedNodePaths) {
  const res = new Set<string>;
  for (const path of paths) {
    const uuids = pathToUUID(rnode, path.path);
    if (path.ioName)
      uuids.push(path.ioName);
    res.add(uuids.join('/'));
  }
  return res;
}
