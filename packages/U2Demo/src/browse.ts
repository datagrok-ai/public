/* The Browse sub-tree under Apps ▸ Dev ▸ U2 Demo: the same 2-level registry the in-app nav shows,
   filled lazily when the platform expands the app node (meta.role: appTreeBrowser). */
import * as DG from 'datagrok-api/dg';
import {DEMO_TREE, DemoLeaf} from './nav';

export function buildBrowseTree(treeNode: DG.TreeViewGroup, open: (leaf: DemoLeaf) => void): void {
  for (const group of DEMO_TREE) {
    const node = treeNode.getOrCreateGroup(group.label, null, false);
    for (const leaf of group.children)
      node.item(leaf.label).onSelected.subscribe(() => open(leaf));
  }
}
