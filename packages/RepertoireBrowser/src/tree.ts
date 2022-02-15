import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class TreeBrowser {// extends DG.JsViewer {
  title: string;
  phyloTreeViewer: DG.Viewer;
  networkViewer: DG.Viewer;

  // constructor() {
  //   super();
  //   this.title = this.string('title', 'Phylogenetic tree');
  // }

  async init(df: DG.DataFrame, mlbView: DG.TableView) {
    const treeCol = df.col('TREE');
    const cloneId = df.col('CLONE');
    let processed = null;

    for (let i = 1; i < treeCol.length; i++) {
      if (DG.Func.find({name: '_newickToDf'}).length == 0)
        return grok.shell.warning('Newick parser is unavailable');

      // TODO: switch call from system to local import.
      const t = await grok.functions.call('PhyloTreeViewer:_newickToDf', {newick: treeCol.get(i), filename: 'nwk'});
      const p = t.col('parent');
      const c = t.col('node');

      const id = cloneId.get(i);
      t.rows.removeAt(0);
      t.columns.addNewString('clone').init((_) => id);
      t.columns.addNewInt('edgeColor').init((_) => i);

      for (let k = 0; k < t.rowCount; k++) {
        const n1 = p.get(k);
        if (n1 == 'root') p.set(k, `root-${id}`);
        else if (n1.startsWith('node-')) p.set(k, `${n1}-${id}`);
        const n2 = c.get(k);
        if (n2.startsWith('node-')) c.set(k, `${n2}-${id}`);
      }
      if (processed == null) processed = t;
      else processed.append(t, true);
    }

    grok.data.linkTables(df, processed, ['clone'], ['clone'], [DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    treeCol.semType = 'newick';
    df.currentRowIdx = 1;

    const tree = DG.Viewer.fromType('PhyloTree', df);
    const network = DG.Viewer.fromType(DG.VIEWER.NETWORK_DIAGRAM, processed, {
      node1: 'node',
      node2: 'parent',
      edgeColorColumnName: 'edgeColor',
    });

    const treeNode = mlbView.dockManager.dock(tree, DG.DOCK_TYPE.DOWN);
    mlbView.dockManager.dock(network, DG.DOCK_TYPE.RIGHT, treeNode);

    network.onEvent('d4-network-diagram-node-click').subscribe((args) => {
      console.warn(['d4-network-diagram-node-click', args, typeof args]);
    });
    network.onEvent('d4-network-diagram-edge-click').subscribe((args) => {
      console.warn(['d4-network-diagram-edge-click', args, typeof args]);
    });

    // this.phyloTreeViewer = tree;
    // this.networkViewer = network;
  }

  get root(): HTMLElement {
    const title = ui.h1(this.title, {style: {'align-self': 'center', 'alignContent': 'center'}});
    if (this.phyloTreeViewer && this.networkViewer) {
      [this.phyloTreeViewer.root, this.networkViewer.root].forEach((v) => v.style.width = 'auto');
      return ui.divV([title, ui.divH([this.phyloTreeViewer.root, this.networkViewer.root])]);
    }
    return title;
  }
}
