import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {setDifference} from '@datagrok-libraries/utils/src/operations';

export class TreeBrowser {// extends DG.JsViewer {
  title: string;
  phyloTreeViewer: DG.Viewer;
  networkViewer: DG.Viewer;
  dataFrame: DG.DataFrame;
  network: DG.DataFrame;
  leafs: TreeLeafs = {};

  // constructor() {
  //   super();
  //   this.title = this.string('title', 'Phylogenetic tree');
  // }

  /**
   * Converts trees column into network diagram data frame.
   * @param {DG.Column} treeCol Column with trees in newick format.
   * @param {DG.Column} cloneId Clone ids.
   * @return {(Promise<DG.DataFrame | undefined>)} Data frame or undefined if the newick converter was failed.
   * node parent  distance  clone edgeColor
   * node-0-349 root-349  0.00  349 1
   * node-1-349 node-0-349  0.02  349 1
   * node-2-349 node-1-349  0.01  349 1
   * ...
   */
  private async _unpackNewickTrees(treeCol: DG.Column, cloneId: DG.Column): Promise<DG.DataFrame | undefined> {
    let processed: DG.DataFrame | undefined;

    if (DG.Func.find({name: '_newickToDf'}).length == 0)
      grok.shell.warning('Newick parser is unavailable');
    else {
      for (let i = 1; i < treeCol.length; i++) {
      // TODO: switch call from system to local import.
        const t: DG.DataFrame = await grok.functions.call(
          'PhyloTreeViewer:_newickToDf',
          {newick: treeCol.get(i), filename: 'nwk'},
        );
        const p = t.col('parent');
        const c = t.col('node');
        const id = cloneId.get(i);

        this.leafs[id] = {index: [i], items: setDifference(c.toList(), p.toList())};

        t.rows.removeAt(0);
        t.columns.addNewString('clone').init((_) => id);
        t.columns.addNewInt('edgeColor').init((_) => i);

        for (let k = 0; k < t.rowCount; k++) {
          const n1 = p.get(k);

          if (n1 == 'root')
            p.set(k, `root-${id}`);
          else if (n1.startsWith('node-'))
            p.set(k, `${n1}-${id}`);

          const n2 = c.get(k);

          if (n2.startsWith('node-'))
            c.set(k, `${n2}-${id}`);
          else {
            if (!this.leafs[n2])
              this.leafs[n2] = {index: [], items: []};

            this.leafs[n2]['index'].push(i);
          }
        }
        if (processed === undefined)
          processed = t;
        else
          processed.append(t, true);
      }
    }
    return processed;
  }

  /**
   * Initializes tree browser.
   * @param {DG.DataFrame} df Table to convert.
   * CLONE ... TREE
   * 123 ... ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930...
   * 123 ... (Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386)...
   * 456 ... (Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:...
   * @param {DG.TableView} mlbView
   * @return {Promise<void>}
   */
  async init(df: DG.DataFrame, mlbView: DG.TableView): Promise<void> {
    this.dataFrame = df;

    const treeCol = df.col('TREE');
    const cloneId = df.col('CLONE');

    const processed: DG.DataFrame = await this._unpackNewickTrees(treeCol, cloneId);

    this.network = processed;
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
      this._onNetworkDiagramNodeClick(args);
    });

    network.onEvent('d4-network-diagram-edge-click').subscribe((args) => {
      this._onNetworkDiagramEdgeClick(args);
    });

    this.phyloTreeViewer = tree;
    // this.networkViewer = network;
  }

  selectTree(index: number) {
    const maxIndex = this.dataFrame.rowCount;
    this.dataFrame.currentRowIdx = index >= maxIndex ? maxIndex - 1 : (index < 0 ? 0 : index);
    this.phyloTreeViewer.onFrameAttached(this.dataFrame);
    console.warn(['selectTree', this.dataFrame.currentRow.idx, this.dataFrame.currentRowIdx]);
  }

  private _selectNode(nodeId: string) {
    let cloneId: string | undefined;

    if (nodeId.startsWith('node') || nodeId.startsWith('root')) {
      const s = nodeId.split('-');
      cloneId = s[s.length - 1];
    }

    const leafs = cloneId ? this.leafs[cloneId] : {index: this.leafs[nodeId], items: nodeId};

    console.warn([nodeId, cloneId, leafs]);

    this.selectTree(leafs['index'][0]);
  }

  /**
   * Called when mouse clicked a node on the network diagram.
   * @param {ClickEventData} data Event data containing the clicked node ID.
   */
  private _onNetworkDiagramNodeClick(data: ClickEventData) {
    this._selectNode((data.args as NodeClickArgs).nodeId);
  }

  /**
   * Called when mouse clicked an edge on the network diagram.
   * @param {ClickEventData} data Event data containing the clicked edge ID.
   */
  private _onNetworkDiagramEdgeClick(data: ClickEventData) {
    const edgeId = (data.args as EdgeClickArgs).edgeId;

    console.warn(['_onNetworkDiagramEdgeClick', edgeId]);

    const parent = this.network.col('parent');
    const nodeId = parent.get(edgeId);

    this._selectNode(nodeId);
    console.warn([edgeId, nodeId]);
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

interface NodeClickArgs {
  nodeId: string;
}
interface EdgeClickArgs {
  edgeId: any;
}
interface ClickEventData {
  dart: any;
  args: NodeClickArgs | EdgeClickArgs;
}

interface TreeMap {
  index: number[];
  items: string[]
}

type TreeLeafs = {[key: string]: TreeMap};
