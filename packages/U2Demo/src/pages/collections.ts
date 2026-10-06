import {
  computed,
  divV, divH, span, h3, button,
  VirtualList, VirtualTree, TreeNode,
} from '@datagrok-libraries/u2';
import {readout} from './common';

export function listsPage(): HTMLElement {
  const list = new VirtualList<number>({
    keyOf: (i) => String(i),
    render: (item, index, row) => {
      row.title = `row ${index}`;
      return span(`item ${item} ${'▪'.repeat(1 + item % 7)}`);
    },
  });
  list.setItems(Array.from({length: 100000}, (_, i) => i));
  list.root.classList.add('u2demo-list');

  return divV([
    span('The list renders only the visible rows — the one below holds 100,000 items.', 'u2demo-hint'),
    h3('VirtualList'), list,
    readout('selectedIndex', computed(() => String(list.selectedIndex.value))),
  ], 'u2demo-page');
}

export function treesPage(): HTMLElement {
  const projects: TreeNode<string>[] = [{
    id: 'demo', label: 'Demo project', children: [
      {id: 'demo/tables', label: 'Tables', children: [
        {id: 'demo/tables/demog', label: 'demog'},
        {id: 'demo/tables/cars', label: 'cars'},
      ]},
      {id: 'demo/queries', label: 'Queries', children: [
        {id: 'demo/queries/orders', label: 'orders by country'},
      ]},
    ],
  }, {
    id: 'remote', label: 'Server (lazy)', children: async () => {
      await new Promise((resolve) => setTimeout(resolve, 600));
      return Array.from({length: 5}, (_, i) => ({id: `remote/${i}`, label: `dataset-${i}.csv`}));
    },
  }];
  const tree = new VirtualTree<string>();
  tree.setRoots(projects);
  tree.root.classList.add('u2demo-list');

  return divV([
    h3('VirtualTree'), tree,
    divH([
      button('Reveal demog', () => void tree.expandPath(['demo', 'demo/tables', 'demo/tables/demog'])),
      span('— expands ancestors (awaiting lazy loads), then selects and reveals.'),
    ], 'u2demo-row'),
    readout('selectedNode', computed(() => tree.selectedNode.value?.label ?? '(none)')),
  ], 'u2demo-page');
}
