import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

// //@ts-ignore
// const oldInit = bio.PhylocanvasGL.prototype.init;
// //@ts-ignore
// bio.PhylocanvasGL.prototype.init = function() {
//   oldInit(arguments);
// };

export class TreeInGridCellApp {
  private view: DG.TableView | null = null;

  private _df: DG.DataFrame;
  get df(): DG.DataFrame { return this._df; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const csv: string = `id,TREE
0,"(a:1.0,b:0.7)ab:0.2;"
1,"((a:1.,b:0.6)ab:0.2,c:1.2);"
2,"((a:1.,b:0.6)ab:0.2,c:1.2)abc:0.1;"
3,"((a:1.,b:0.6)ab:0.2,(c:1.2,d:0.3)abc:0.1;"
4,"(a:1.0,b:0.7)ab:0.2;"
5,"((a:1.,b:0.6)ab:0.2,c:1.2);"
6,"((a:1.,b:0.6)ab:0.2,c:1.2)abc:0.1;"
7,"((a:1.,b:0.6)ab:0.2,(c:1.2,d:0.3)abc:0.1;"
8,"(a:1.0,b:0.7)ab:0.2;"
9,"((a:1.,b:0.6)ab:0.2,c:1.2);"
10,"((a:1.,b:0.6)ab:0.2,c:1.2)abc:0.1;"
11,"((a:1.,b:0.6)ab:0.2,(c:1.2,d:0.3)abc:0.1;"
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
      this.view.path = '/apps/PhyloTreeViewer/TreeInGridCell';

      this.view.grid.props.rowHeight = 120;
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
      gCell.gridColumn.column.name == 'TREE' &&
      gCell.tableRowIndex !== undefined && [0, 2, 4].includes(gCell.tableRowIndex!)
    ) {
      const nwk = gCell.cell.value;

      const bd = args.bounds;
      const gCtx: CanvasRenderingContext2D = args.g;
      gCtx.save();
      try {
        gCtx.beginPath();
        gCtx.rect(bd.x, bd.y, bd.width, bd.height);
        gCtx.clip();

        // renderTree(cc, gc.cell.value, bd);
        const cellDiv = ui.div([], {style: {width: `100px`, height: `100px`, backgroundColor: '#F0F0FF'}});
        const pcgl: bio.PhylocanvasGL = new bio.PhylocanvasGL(cellDiv, {
          size: {width: 250, height: 100},
          //size: {width: bd.width, height: bd.height},
          treeType: bio.TreeTypes.Rectangular,
          nodeSize: 5,
          nodeShape: bio.Shapes.Circle,
          source: nwk,
        });
        pcgl.deck.setProps({useDevicePixels: true});
        pcgl.view.style.backgroundImage = 'none';
        try {
          // //@ts-ignore
          // pcgl.resume();
          //@ts-ignore
          pcgl.render();

          //@ts-ignore
          pcgl.deck.setProps({glOptions: {preserveDrawingBuffer: true}});
          //@ts-ignore
          pcgl.deck.redraw('true');

          //@ts-ignore
          const png = pcgl.exportPNG();

          let k = 11;

          // //@ts-ignore
          // pcgl.deck.animationLoop._createWebGLContext();
          // //@ts-ignore
          // pcgl.deck.animationLoop.redraw();


          //@ts-ignore
          const pCtx: CanvasRenderingContext2D = pcgl.deck.canvas.getContext('2d');

          pCtx.fillStyle = '#FFF0F0';
          // pCtx.fillRect(10, 10, 50, 50);
          pCtx.fillRect(0, 0, pCtx.canvas.width, pCtx.canvas.height);
          //pcgl.deck.redraw('true');

          pCtx.strokeStyle = '#008000';
          //@ts-ignore
          const a = pcgl.getDrawingArea();
          pCtx.strokeRect(a.left, a.top, a.width, a.height);

          gCtx.fillStyle = 'blue';
          gCtx.fillRect(200, 30, 100, 70);

          // TODO: Force draw tree
          gCtx.strokeStyle = 'red';
          gCtx.lineWidth = 2;
          gCtx.beginPath();
          // cc.moveTo(bd.left, bd.top);
          // cc.lineTo(bd.right, bd.bottom);
          gCtx.moveTo(4, 4);
          gCtx.lineTo(8, 8);
          gCtx.closePath();

          // let k = 11;

          //@ts-ignore
          gCtx.drawImage(pCtx.canvas, bd.x+5, bd.y+5);
        } finally {
          pcgl.destroy();
          args.preventDefault();
        }
      } finally {
        gCtx.restore();
      }
    }
  }
}