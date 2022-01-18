import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";


export class TreeBrowser {
  async init(df) {
    let treeCol = df.col('TREE');
    let cloneId = df.col('CLONE');
    let processed = null;

    for (let i = 1; i < treeCol.length; i++) {
      if (DG.Func.find({name: '_newickToDf'}).length == 0) return grok.shell.warning('Newick parser is unavailable');
      let t = await grok.functions.call('PhyloTreeViewer:_newickToDf', { newick: treeCol.get(i), filename: 'nwk' });
      let p = t.col('parent');
      let c = t.col('node');

      let id = cloneId.get(i);
      t.rows.removeAt(0);
      t.columns.addNewString('clone').init((_) => id);
      t.columns.addNewInt('edgeColor').init((_) => i);

      for (let k = 0; k < t.rowCount; k++) {
        let n1 = p.get(k);
        if (n1 == 'root') p.set(k, `root-${id}`)
        else if (n1.startsWith('node-')) p.set(k, `${n1}-${id}`);
        let n2 = c.get(k);
        if (n2.startsWith('node-')) c.set(k, `${n2}-${id}`);
      }
      if (processed == null) processed = t;
      else processed.append(t, true);
    }

    grok.data.linkTables(df, processed, ['clone'], ['clone'], [DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    treeCol.semType = 'newick';
    df.currentRow = 1;

    let view = grok.shell.addTableView(df);
    let grid = DG.Viewer.grid(processed);
    let tree = DG.Viewer.fromType('PhyloTree', df);
    let trees = DG.Viewer.fromType(DG.VIEWER.NETWORK_DIAGRAM, processed, {
      node1: 'node',
      node2: 'parent',
      edgeColorColumnName: 'edgeColor'
    });  

    view.addViewer(tree);
    view.dockManager.dock(tree, DG.DOCK_TYPE.RIGHT);

    view.addViewer(grid);
    let gridNode = view.dockManager.dock(grid, DG.DOCK_TYPE.DOWN);

    view.dockManager.dock(trees, DG.DOCK_TYPE.RIGHT, gridNode);
  }
}
