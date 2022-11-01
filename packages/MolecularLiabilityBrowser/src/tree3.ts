import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class Tree3Browser {
  async init(df: DG.DataFrame, mlbView: DG.TableView) {
    const treeCol = df.getCol('TREE');
    const cloneId = df.getCol('CLONE');
    let processed = null;
    for (let i = 1; i < treeCol.length; i++) {
      if (DG.Func.find({name: '_newickToDf'}).length == 0) return grok.shell.warning('Newick parser is unavailable');
      const t = await grok.functions.call('PhyloTreeViewer:_newickToDf',
        {newick: treeCol.get(i), name: cloneId.get(i)});
      const p = t.col('parent');
      const c = t.col('node');
      const id = cloneId.get(i);
      t.rows.removeAt(0);
      t.columns.addNewString('clone').init((_: number) => id);
      t.columns.addNewInt('edgeColor').init((_: number) => i);
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
    //let view = grok.shell.addTableView(df);
    const grid = DG.Viewer.grid(processed);
    const tree = DG.Viewer.fromType('PhyloTree', df);
    // let trees = DG.Viewer.fromType(DG.VIEWER.NETWORK_DIAGRAM, processed, {
    //   node1: 'node',
    //   node2: 'parent',
    //   edgeColorColumnName: 'edgeColor'
    // });
    //view.addViewer(tree);
    const treeNode = mlbView.dockManager.dock(tree, DG.DOCK_TYPE.RIGHT);
    //mlbView.addViewer(grid);
    //let gridNode = mlbView.dockManager.dock(grid, DG.DOCK_TYPE.DOWN);
    //mlbView.dockManager.dock(trees, DG.DOCK_TYPE.RIGHT, treeNode);
  }
}
