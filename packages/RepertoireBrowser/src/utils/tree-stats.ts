//@ts-ignore
import {Utils, Newick} from '@phylocanvas/phylocanvas.gl';

interface TreeNodeInfo {
    branchLength: number,
    id: string,
    isLeaf: boolean,
    name: string,
}

interface TreeStats {
    totalLeaves: number;
    leavesIntersected: number;
    totalSubtreeLength: number;
}

type TreeStatsKeys = keyof TreeStats;
// eslint-disable-next-line no-unused-vars
type TreeStatColumns = {[key in TreeStatsKeys]: number[]};
type IdParser = (name: string) => string;

export class TreeAnalyzer {
  protected items: Set<string>;
  protected parseId: IdParser;

  constructor(items: string[] = [], parser: IdParser = (s) => s) {
    this.items = new Set(items);
    this.parseId = parser;
  }

  protected _findIntersection(nodes: {[key: string]: TreeNodeInfo}): string[] {
    const isect = [];

    for (const [id, node] of Object.entries(nodes)) {
      if (!node.isLeaf)
        continue;

      const parsedId = this.parseId(id);

      if (this.items.has(parsedId))
        isect.push(parsedId);
    }
    return isect;
  }

  private _traverseTree(nwk: string): TreeStats {
    const tree = Newick.parse_newick(nwk);
    const traverse = Utils.treeTraversal(tree);
    const root = traverse.rootNode;
    const isect = this._findIntersection(traverse.nodeById);
    const stats = {
      totalLeaves: root.totalLeaves,
      leavesIntersected: isect.length,
      totalSubtreeLength: root.totalSubtreeLength,
    };
    return stats;
  }

  analyze(trees: string[]): TreeStatColumns {
    const regexp = new RegExp(/^(\(*[^:]+:[^,]+\)*)+$/);
    // eslint-disable-next-line no-unused-vars
    const stats: TreeStatColumns = {
      totalLeaves: [],
      leavesIntersected: [],
      totalSubtreeLength: [],
    };
    const statsKeys = Object.keys(stats);

    for (const k of statsKeys)
      stats[k] = [];

    for (const nwk of trees) {
      const s = regexp.test(nwk.trim()) ? this._traverseTree(nwk) : statsKeys;

      for (const k of statsKeys)
        stats[k].push(s[k]);
    }
    return stats;
  }
}
