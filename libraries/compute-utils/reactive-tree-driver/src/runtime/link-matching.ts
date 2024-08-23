import {isNonRefSelector, LinkNonRefSelectors, LinkParsed, LinkRefSelectors, refSelectorAdjacent, refSelectorAll, refSelectorDirection, refSelectorFindOne} from '../config/LinkSpec';
import {PipelineLinkConfigurationBase} from '../config/PipelineConfiguration';
import {isFuncCallState} from '../config/PipelineInstance';
import {BaseTree, NodePath, TreeNode} from '../data/BaseTree';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {indexFromEnd} from '../utils';
import {StateTree} from './StateTree';
import {isFuncCallNode, StateTreeNode} from './StateTreeNodes';

export type LinkSpec = PipelineLinkConfigurationBase<LinkParsed[]>;

type matchedIO = {
  path: Readonly<NodePath>;
  ioName?: string;
}

type MatchedNodePaths = Readonly<Array<matchedIO>>;

export type MatchInfo = {
  basePath?: Readonly<NodePath>;
  inputs: Record<string, MatchedNodePaths>;
  outputs: Record<string, MatchedNodePaths>;
}

export function matchLink(state: StateTree, address: NodePath, spec: LinkSpec): MatchInfo[] | undefined {
  const [rnode] = state.find(((_node, path) => BaseTree.isNodeAddressEq(path, address))) ?? [];
  if (rnode == null)
    return;
  const basePaths = spec.base?.length ? expandLinkBase(rnode, spec.base[0]) : undefined;
  const baseName = spec.base?.length ? spec.base[0].name : undefined;
  if (basePaths == null) {
    const matchInfo = matchLinkInstance(rnode, spec);
    return matchInfo ? [matchInfo] : undefined;
  } else {
    const instances = basePaths.map((basePath) => matchLinkInstance(rnode, spec, basePath, baseName)).filter((x) => !!x);
    return instances.length ? instances : undefined;
  }
}

function matchLinkInstance(rnode: TreeNode<StateTreeNode>, spec: LinkSpec, base?: Readonly<matchedIO>, baseName?: string): MatchInfo | undefined {
  const matchInfo: MatchInfo = {
    basePath: base?.path,
    inputs: {},
    outputs: {},
  };
  const currentIO: Record<string, MatchedNodePaths> = {};
  if (base && baseName)
    currentIO[baseName] = [base];

  const ioData = [...spec.from.map((item) => ['inputs', item] as const), ...spec.to.map((item) => ['outputs', item] as const)];
  for (const [kind, io] of ioData) {
    const paths = matchLinkIO(rnode, currentIO, io);
    if (paths.length == 0)
      return;
    if (currentIO[io.name] != null)
      throw new Error(`Duplicate io name ${io.name} in link ${rnode.getItem().config.id}`);
    currentIO[io.name] = paths;
    matchInfo[kind][io.name] = paths;
  }
  return matchInfo;
}

function expandLinkBase(rnode: TreeNode<StateTreeNode>, baseLink: LinkParsed) {
  const traverse = buildTraverseD([] as Readonly<NodePath>, (pnode: TreeNode<StateTreeNode>, path, level?: number) => {
    const segment = baseLink.segments[level!];
    if (!segment)
      return [] as const;
    const {selector, ids} = segment;
    const nextNodes = matchNonRefSegment(pnode, ids, selector as LinkNonRefSelectors);
    return nextNodes.map(([idx, node]) => [node.item, [...path, {id: node.id, idx}], level!+1] as const);
  }, 0);

  const basePaths = traverse(rnode, (acc, _node, path) => {
    if (path.length < baseLink.segments.length)
      return acc;
    return [...acc, {path}] as const;
  }, [] as MatchedNodePaths);

  return basePaths;
}


function matchLinkIO(rnode: TreeNode<StateTreeNode>, currentIO: Record<string, MatchedNodePaths>, parsedLink: LinkParsed): MatchedNodePaths {
  const traverse = buildTraverseD([] as Readonly<NodePath>, (pnode: TreeNode<StateTreeNode>, path, level?: number) => {
    const segment = parsedLink.segments[level!];
    if (!segment)
      return [] as const;
    const {ref, selector, ids, stopIds} = segment;
    if (isNonRefSelector(selector)) {
      const nextNodes = matchNonRefSegment(pnode, ids, selector);
      return nextNodes.map(([idx, node]) => [node.item, [...path, {id: node.id, idx}], level!+1] as const);
    } else {
      let originIdx = undefined;
      const selDirection = refSelectorDirection(selector);
      if (ref) {
        const io = currentIO[ref];
        if (io == null)
          throw new Error(`Node ${rnode.getItem().config.id} referenced unknown io ${ref} in ${parsedLink.name}`);
        const refOrigin = selDirection === 'before' ? io[0] : indexFromEnd(io)!;
        if (!BaseTree.isNodeChild(path, refOrigin.path))
          throw new Error(`Node ${rnode.getItem().config.id} reference path ${JSON.stringify(refOrigin.path)} is different from current ${JSON.stringify(path)}`);
        originIdx = refOrigin.path[level!].idx;
      }
      const nextNodes = matchRefSegment(pnode, ids, selector, originIdx, stopIds);
      return nextNodes.map(([idx, node]) => [node.item, [...path, {id: node.id, idx}], level!+1] as const);
    }
  }, 0);

  const paths = traverse(rnode, (acc, node, path) => {
    if (path.length === parsedLink.segments.length - 1) {
      const ioSegment = indexFromEnd(parsedLink.segments)!;
      const ioName = ioSegment.ids[0];
      const item = node.getItem();
      if (isFuncCallNode(item)) {
        const io = item.config.io!;
        const idx = io.findIndex((ioItem) => ioItem.id === ioSegment.ids[0]);
        if (idx >= 0)
          return [...acc, {path, ioName}] as const;
      }
      return acc;
    }
    return acc;
  }, [] as MatchedNodePaths);
  return paths;
}

function matchNonRefSegment(pnode: TreeNode<StateTreeNode>, ids: string[], selector: LinkNonRefSelectors) {
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

function matchRefSegment(pnode: TreeNode<StateTreeNode>, ids: string[], selector: LinkRefSelectors, originIdx?: number, stopIds?: string[]) {
  const idsSet = new Set(ids);
  const stopSet = new Set(stopIds);
  const allNodeEntries = [...pnode.getChildren().entries()];

  if (selector === 'same') {
    if (originIdx == null)
      return [];
    const target = allNodeEntries[originIdx];
    if (idsSet.has(target[1].id))
      return [target];
    return [];
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
    if (selDirection === 'before')
      return [indexFromEnd(matchingNodes)!];
    else
      return [matchingNodes[0]];
  }
  throw new Error(`Unknown segement mode ${selector}`);
}
