import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';
import * as u from '@datagrok-libraries/utils';

import {PickingInfo} from '@deck.gl/core/typed';
import {TooltipContent} from '@deck.gl/core/typed/lib/tooltip';
import {Rect} from '@deck.gl/core/typed/passes/layers-pass';

// //@ts-ignore
// const oldInit = bio.PhylocanvasGL.prototype.init;
// //@ts-ignore
// bio.PhylocanvasGL.prototype.init = function() {
//   oldInit(arguments);
// };

export class TreeInGridCellApp {
  private phylocanvasGlSvc: bio.PhylocanvasGlServiceBase | null = null;

  private view: DG.TableView | null = null;

  private _df: DG.DataFrame;
  get df(): DG.DataFrame { return this._df; }

  async init(): Promise<void> {
    this.phylocanvasGlSvc = await bio.getPhylocanvasGlService();
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const csv: string = `id,TREE
0,"(a:1.0,b:1.0);","text of tree zero size" 
1,"(a:1.0);","tree with one leaf"
2,"(a:1.0,b:1.4);","two leaves in tree"
3,"((a:1.0,b:1.2)ab:0.2,c:1.0);", "tree of three leaves and one internal node"
4,"((a:1.0,b:1.2)ab:0.2,(c:0.6,d:1.3)cd:0.1);","some text for four leaves"
5,"(((a:1.0,b:1.2)ab:0.2,(c:0.6,d:1.3)cd:0.1)abcd:0.4,e:0.8);","I was too lazy to invent text here"
6,"(((a:1.0,b:1.2)ab:0.2,(c:0.6,d:1.3)cd:0.1)abcd:0.4,(e:0.8,f:0.6)ef:0.4);","short text for large tree"
7,"((a:1.,b:0.6)ab:0.2,(c:1.2,d:0.3)cd:0.1)abcd:0.1;","text7"
8,"(a:1.0,b:0.7)ab:0.2;", "text8"
9,"((a:1.,b:0.6)ab:0.2,c:1.2);","text9"
10,"((a:1.,b:0.6)ab:0.2,(c:0.2,d:0.3)cd:0.1);","text10"
11,"((a:1.,b:0.6)ab:0.2,(c:1.2,d:0.3)cd:0.1)abcd:0.1;","text11"
`;
    const df = DG.DataFrame.fromCsv(csv);

    await this.setData(df);
  }

  async setData(df: DG.DataFrame): Promise<void> {
    await this.destroyView();

    this._df = df;

    await this.buildView();
  }

  // -- View --

  async destroyView(): Promise<void> {
    if (this.view) {
      this.view.close();
      this.view = null;
    }
  }

  async buildView(): Promise<void> {
    if (!this.view) {
      this.view = grok.shell.addTableView(this.df);
      this.view.path = this.view.basePath = '/apps/PhyloTreeViewer/TreeInGridCell';

      this.view.grid.props.rowHeight = 150;
      const treeGCol: DG.GridColumn = this.view.grid.columns.byName('TREE')!;
      treeGCol.width = 300;
      this.view.grid.onCellRender.subscribe(this.gridOnCellRender.bind(this));
    }
  }

  // -- Handle events --

  gridOnCellRender(args: DG.GridCellRenderArgs): void {
    const gCell = args.cell;
    if (
      gCell.isTableCell && gCell.gridColumn.column &&
      gCell.gridColumn.column.name == 'TREE'
      // && gCell.tableRowIndex !== undefined && [0, 2, 4].includes(gCell.tableRowIndex!)
    ) {
      try {
        const name: string = `gridRow = ${gCell.gridRow}`;
        const bd = args.bounds;
        const gCtx: CanvasRenderingContext2D = args.g;
        //console.debug('PTV: TreeInGridCell.gridOnCellRender() start ' + `name: ${name}, bd: ${rectToString(bd)} `);

        const nwkStr: string = gCell.cell.value;
        const nwkRoot: bio.NodeType = bio.Newick.parse_newick(nwkStr);

        const nodeShape: string = bio.Shapes.Circle;
        const treeType: string = bio.TreeTypes.Rectangular;
        this.phylocanvasGlSvc!.render({
          name: name,
          backColor: gCell.grid.props.backColor,
          props: {
            size: {width: bd.width - 2, height: bd.height - 3},
            treeType: treeType,
            nodeSize: 3,
            nodeShape: nodeShape,
            treeToCanvasRatio: 0.95,
            padding: 5,
            source: {type: 'biojs', data: nwkRoot},
          },
          onAfterRender: (canvas: HTMLCanvasElement) => {
            this.phylocanvasGlSvc!.renderOnGridCell(gCtx, bd, gCell, canvas);
          }
        }, gCell.tableRow?.idx);
      } catch (err) {
        console.error(u.errorToConsole(err));
      } finally {
        args.preventDefault();
        //console.debug('PTV: TreeInGridCell.gridOnCellRender() end ' + `name: ${name}`);
      }
    }
  }
}