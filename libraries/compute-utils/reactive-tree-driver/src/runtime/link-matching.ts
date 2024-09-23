import { isNonRefSelector, LinkNonRefSelectors, LinkParsed, LinkRefSelectors, refSelectorAdjacent, refSelectorAll, refSelectorDirection, refSelectorFindOne, refSelectorIncludesOrigin } from '../config/LinkSpec';
import {PipelineLinkConfigurationBase} from '../config/PipelineConfiguration';
import {BaseTree, NodePath, TreeNode} from '../data/BaseTree';
import {buildTraverseD} from '../data/graph-traverse-utils';
import {indexFromEnd} from '../utils';
import {StateTree} from './StateTree';
import {StateTreeNode} from './StateTreeNodes';

export type LinkSpec = PipelineLinkConfigurationBase<LinkParsed[]>;

type NodePathsRO = Readonly<Array<Readonly<NodePath>>>;

export type MatchInfo = {
  inputs: Record<string, NodePathsRO>;
  outputs: Record<string, NodePathsRO>;
}

export function matchLink(state: StateTree, address: NodePath, spec: LinkSpec): MatchInfo | undefined {
  const [rnode] = state.find(((_node, path) => BaseTree.isNodeAddressEq(path, address))) ?? [];
  if (rnode == null)
    return;
}

function expandBase() {
  // TODO
}

function matchLinkIO(rnode: TreeNode<StateTreeNode>, currentIO: Record<string, NodePathsRO>, parsedLink: LinkParsed): NodePathsRO {
  const traverse = buildTraverseD([] as Readonly<NodePath>, (pnode: TreeNode<StateTreeNode>, path, level?: number) => {
    const segment = parsedLink.segments[level!];
    if (!segment)
      return [] as const;
    const {ref, selector, id, args} = segment;
    if (isNonRefSelector(selector)) {
      const nextNodes = matchNonRefSegment(pnode, id, selector);
      return nextNodes.map(([idx, node]) => [node.item, [...path, {id: node.id, idx}], level!+1] as const);
    } else {
      let originIdx = undefined;
      const selDirection = refSelectorDirection(selector);
      if (ref) {
        const io = currentIO[ref];
        if (io == null)
          throw new Error(`Node ${rnode.getItem().config.id} referenced unknown input ${ref} in output ${parsedLink.name}`);
        const refOrigin = selDirection === 'before' ? io[0] : indexFromEnd(io)!;
        if (!BaseTree.isNodeAddressPrefix(refOrigin, path, level!))
          throw new Error(`Node ${rnode.getItem().config.id} referenced input in different segment from output ${parsedLink.name}`);
        originIdx = refOrigin[level!].idx;
      }
      const stopIds = (args) ? new Set(segment.args) : undefined;
      const nextNodes = matchRefSegment(pnode, segment.id, selector, originIdx, stopIds);
      return nextNodes.map(([idx, node]) => [node.item, [...path, {id: node.id, idx}], level!+1] as const);
    }
  }, parsedLink.segments.length);

  const paths = traverse(rnode, (acc, _node, path) => {
    if (path.length < parsedLink.segments.length)
      return acc;
    return [...acc, path] as const;
  }, [] as NodePathsRO);
  return paths;
}

function matchNonRefSegment(pnode: TreeNode<StateTreeNode>, id: string, mode: LinkNonRefSelectors) {
  const matchingNodes = [...pnode.getChildren().entries()].filter(([, c]) => c.id === id);
  if (matchingNodes.length === 0)
    return [];
  if (mode === 'all')
    return matchingNodes;
  if (mode === 'first')
    return [matchingNodes[0]];
  if (mode === 'last')
    return [indexFromEnd(matchingNodes)!];
  throw new Error(`Unknown segement mode ${mode}`);
}

function matchRefSegment(pnode: TreeNode<StateTreeNode>, id: string, selector: LinkRefSelectors, originIdx?: number, stopIds?: Set<string>) {
  const selDirection = refSelectorDirection(selector);
  const includesOrigin = refSelectorIncludesOrigin(selector);

  const allNodeEntries = [...pnode.getChildren().entries()];

  const partitionedByOrigin = originIdx == null ?
    allNodeEntries :
    allNodeEntries.filter(([idx]) => {
      if (includesOrigin)
        return selDirection === 'before' ? idx <= originIdx : idx >= originIdx;
      else
        return selDirection === 'before' ? idx < originIdx : idx > originIdx;
    });

  let matchingNodes = partitionedByOrigin;
  if (stopIds?.size) {
    const stopIdx = selDirection === 'before' ?
      partitionedByOrigin.findLastIndex(([, node]) => stopIds.has(node.id)) :
      partitionedByOrigin.findIndex(([, node]) => stopIds.has(node.id));
    if (stopIdx > 0)
      matchingNodes = partitionedByOrigin.slice(0, stopIdx+1);
  }
  matchingNodes = matchingNodes.filter(([, c]) => c.id === id);

  if (matchingNodes.length === 0)
    return [];
  if (refSelectorAll(selector)) {
    return matchingNodes;
  }
  if (selector === 'same') {
    const idx = matchingNodes[0][0];
    return idx == originIdx ? [matchingNodes[0]] : [];
  }
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
    if (selDirection === 'before') {
      return [indexFromEnd(matchingNodes)!];
    } else {
      return [matchingNodes[0]];
    }
  }
  throw new Error(`Unknown segement mode ${selector}`);
}
