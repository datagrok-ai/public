import {Control, VirtualTree, component, computed, notify, span} from '@datagrok-libraries/u2';
import type {TreeNode} from '@datagrok-libraries/u2';
import './tree.css';

export function filesTree(roots: Record<string, string[]>, list: (path: string) => Promise<TreeNode[]>): Control {
  return component(() => {
    const tree = new VirtualTree({contextActions: (node) => [
      {text: 'Open', icon: 'folder-open', run: () => {
        notify.info(`Opening ${node.label}`);
      }},
      {text: 'Copy path', icon: 'copy', run: () => {
        notify.info(`Copied ${node.id}`);
      }},
    ]});
    tree.root.classList.add('files-tree');
    tree.setRoots(Object.entries(roots).map(([root, folders]) => ({id: root, label: root, children: folders.map(
      (folder) => ({id: `${root}/${folder}`, label: folder, children: () => list(`${root}/${folder}`)}))})));
    return [tree, span(computed(() => tree.selectedNode.value?.id ?? 'Nothing selected'), 'files-selected')];
  }, {name: 'files-tree'});
}

// Live demo
async function listFolder(path: string): Promise<TreeNode[]> {
  await new Promise((resolve) => setTimeout(resolve, 400));
  if (path.endsWith('Restricted'))
    throw new Error('You have no access to this folder');
  return Array.from({length: 5}, (_, i) => ({id: `${path}/plate-${i + 1}.csv`, label: `plate-${i + 1}.csv`}));
}

export function render(host: HTMLElement): Control {
  const demo = filesTree({Screening: ['2026-10', '2026-09'], Assays: ['Kinase', 'Restricted']}, listFolder);
  host.append(demo.root);
  return demo;
}
