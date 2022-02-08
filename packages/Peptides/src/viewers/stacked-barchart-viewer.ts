import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {scaleBand, scaleLinear, text} from 'd3';
import {ChemPalette} from '../utils/chem-palette';
import * as rxjs from 'rxjs';
const cp = new ChemPalette('grok');

//TODO: the function should not accept promise. Await the parameters where it is used
export function addViewerToHeader(grid: DG.Grid, viewer: Promise<DG.Widget>) {
  viewer.then((viewer) => {
    const barchart = viewer as StackedBarChart; //TODO: accept specifically StackedBarChart object
    // The following event makes the barchart interactive
    rxjs.fromEvent(grid.overlay, 'click').subscribe((mm: any) => {
      mm = mm as MouseEvent;
      const cell = grid.hitTest(mm.offsetX, mm.offsetY);

      if (cell !== null && cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids')
        barchart.highlight(cell, mm.offsetX, mm.offsetY);
      else
        return;

      if (cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids') 
        barchart.beginSelection(mm);
      else 
        barchart.unhighlight();
    });

    rxjs.fromEvent(grid.overlay, 'mouseout').subscribe((_: any) => {
      barchart.unhighlight();
    });

    barchart.tableCanvas = grid.canvas;
    grid.setOptions({'colHeaderHeight': 130});
    grid.onCellTooltip((cell, x, y) => {
      if (cell.tableColumn) {
        if (['aminoAcids', 'alignedSequence'].includes(cell.tableColumn.semType) ) {
          if ( !cell.isColHeader) {
            cp.showTooltip(cell, x, y);
            return true;
          } else {
            if (barchart.highlighted) {
              let elements: HTMLElement[] = [];
              elements = elements.concat([ui.divText(barchart.highlighted.aaName)]);
              ui.tooltip.show(ui.divV(elements), x, y);
            }
            return true;
          }
        }
      }
    });
    
    grid.onCellRender.subscribe((args) => {
      args.g.save();
      args.g.beginPath();
      args.g.rect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
      args.g.clip();

      if (args.cell.isColHeader && barchart.aminoColumnNames.includes(args.cell.gridColumn.name)) {
        barchart.renderBarToCanvas(
          args.g,
          args.cell,
          args.bounds.x,
          args.bounds.y,
          args.bounds.width,
          args.bounds.height,
        );
        args.preventDefault();
      }
      args.g.restore();
    });
  });
}


export class StackedBarChart extends DG.JsViewer {
    public dataEmptyAA: string;
    highlighted: {'colName' : string, 'aaName' : string} | null = null;
    private ord: { [Key: string]: number; } = {};
    private data: {'name': string, 'data': {'name': string, 'count': number, 'selectedCount': number}[]}[] = [];
    private selectionMode: boolean = false;
    public aminoColumnNames: string[] = [];
    private aminoColumnIndices: {[Key: string]: number} = {};
    private aggregatedTables: {[Key: string]: DG.DataFrame} = {};
    private aggregatedTablesUnselected: {[Key: string]: DG.DataFrame} = {};
    private max = 0;
    private barStats: {[Key: string]: {'name': string, 'count': number, 'selectedCount': number}[]} = {};
    private registered: {[Key: string]: DG.GridCell} = {};
    tableCanvas: HTMLCanvasElement | undefined;

    constructor() {
      super();
      this.dataEmptyAA = this.string('dataEmptyAA', '-');
    }

    init(): void {
      const groups: {[key: string]: string[]} = {
        'yellow': ['C', 'U'],
        'red': ['G', 'P'],
        'all_green': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
        'light_blue': ['R', 'H', 'K'],
        'dark_blue': ['D', 'E'],
        'orange': ['S', 'T', 'N', 'Q'],
      };

      let i = 0;
      for (const value of Object.values(groups)) {
        i++;
        for (const obj of value)
          this.ord[obj] = i;
      }
      this.data = [];
      this.aminoColumnNames = [];
    }

    // Stream subscriptions
    onTableAttached(): void {
      this.init();
      if (this.dataFrame) {
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
        this.computeData(this.dataFrame);
      }
    }

    // Cancel subscriptions when the viewer is detached
    detach(): void {
      this.subs.forEach((sub) => sub.unsubscribe());
    }

    computeData(df: DG.DataFrame): void {
      this.data = [];
      this.aminoColumnNames = [];
      this.aminoColumnIndices = {};

      df.columns.names().forEach((name: string) => {
        if (df.getCol(name).semType === 'aminoAcids' &&
            !df.getCol(name).categories.includes('COOH') &&
            !df.getCol(name).categories.includes('NH2')) {
          this.aminoColumnIndices[name] = this.aminoColumnNames.length + 1;
          this.aminoColumnNames.push(name);
        }
      });

      function getSelectedFilteredBuffer(df: DG.DataFrame) {
        const buf1 = df.selection.getBuffer();
        const buf2 = df.filter.getBuffer();
        const resbuf = new Int32Array(buf1.length);

        for (let i = 0; i < buf2.length; i++)
          resbuf[i] = buf1[i] & buf2[i];

        return resbuf.buffer;
      }

      this.aggregatedTables = {};
      this.aggregatedTablesUnselected = {};
      //TODO: optimize it, why store so many tables?
      const mask = DG.BitSet.fromBytes(getSelectedFilteredBuffer(df), df.rowCount);
      this.selectionMode = mask.trueCount !== df.filter.trueCount;
      this.aminoColumnNames.forEach((name) => {
        this.aggregatedTables[name] = df
          .groupBy([name])
          .whereRowMask(df.filter)
          .add('count', name, `${name}_count`)
          .aggregate();
        
        if (mask.trueCount !== df.filter.trueCount) {
          const aggregatedMask = DG.BitSet.fromBytes(getSelectedFilteredBuffer(df), df.rowCount);
          this.aggregatedTablesUnselected[name] = df
            .groupBy([name])
            .whereRowMask(aggregatedMask)
            .add('count', name, `${name}_count`)
            .aggregate();
        } 
      });

      this.data = [];
      this.barStats = {};

      for (const [name, df] of Object.entries(this.aggregatedTables)) {
        const colObj: {
          'name': string,
          'data': { 'name': string, 'count': number, 'selectedCount': number }[],
        } = {'name': name, 'data': []};
        this.barStats[colObj['name']] = colObj['data'];
        this.data.push(colObj);
        
        for (let i = 0; i < df.rowCount; i++) {
          const amino = df.getCol(name).get(i);
          const aminoCount = df.getCol(`${name}_count`).get(i);
          const aminoObj = {'name': amino, 'count': aminoCount, 'selectedCount': 0};

          if (!amino || amino === this.dataEmptyAA)
            continue;

          colObj['data'].push(aminoObj);
          for (let j = 0; j < this.aggregatedTablesUnselected[name].rowCount; j++) {
            const unsAmino = this.aggregatedTablesUnselected[name].getCol(`${name}`).get(j);
            if (unsAmino == amino) {
              aminoObj['selectedCount'] = this.aggregatedTablesUnselected[name]
                .getCol(`${name}_count`)
                .get(j);
              break;
            }
          }
        }

        colObj['data'] = colObj['data'].sort((o1, o2) => {
          if (this.ord[o1['name']] > this.ord[o2['name']])
            return -1;

          if (this.ord[o1['name']] < this.ord[o2['name']])
            return 1;

          return 0;
        });
      }
      this.max = df.filter.trueCount;
    }

    renderBarToCanvas(
      g: CanvasRenderingContext2D,
      cell: DG.GridCell,
      x: number,
      y: number,
      w: number,
      h: number,
    ): void {
      const name = cell.tableColumn!.name;
      const colNameSize = g.measureText(name).width;
      const barData = this.barStats[name];
      const margin = 0.2;
      const innerMargin = 0.02;
      const selectLineRatio = 0.1;
      let sum = 0;

      barData.forEach((obj) => {
        sum += obj['count'];
      });

      x = x + w * margin;
      y = y + h * margin / 4;
      w = w - w * margin * 2 - 10;
      h = h - h * margin;
      g.fillStyle = 'black';
      g.textBaseline = 'top';
      g.font = `${h * margin / 2}px`;
      g.fillText(name, x + (w - colNameSize) / 2, y + h + h * margin / 4);

      barData.forEach((obj, index) => {
        const sBarHeight = h * obj['count'] / this.max;
        const gapSize = sBarHeight * innerMargin;
        const [, aarOuter,,] = cp.getColorAAPivot(obj['name']);
        const textSize = g.measureText(aarOuter);
        const fontSize = 11;
        const leftMargin = (w - (aarOuter.length > 1 ? fontSize : textSize.width - 8)) / 2;
        const subBartHeight = sBarHeight - gapSize;
        const verticalShift = (this.max - sum) / this.max;
        const absX = x + leftMargin;
        const absY = y + h * verticalShift + gapSize / 2 + subBartHeight / 2 + (aarOuter.length == 1 ? + 4 : 0);

        g.strokeStyle = cp.getColor(obj['name']);
        g.fillStyle = cp.getColor(obj['name']);
        if (textSize.width <= subBartHeight) {
          const origTransform = g.getTransform();

          g.strokeRect(x, y + h * verticalShift + gapSize / 2, w, subBartHeight);

          g.font = `${fontSize}px monospace`;
          g.textAlign = 'center';
          g.textBaseline = 'bottom';

          if (aarOuter.length > 1) {
            g.translate(absX, absY);
            g.rotate(Math.PI / 2);
            g.translate(-absX, -absY);
          }

          g.fillText(aarOuter, absX, absY);
          g.setTransform(origTransform);
        } else {
          g.fillRect(x, y + h * verticalShift + gapSize / 2, w, subBartHeight);
        }

        if (this.selectionMode && obj['selectedCount'] > 0) {
          g.fillStyle = 'rgb(255,165,0)';
          g.fillRect(
            x - w * selectLineRatio * 1.5,
            y + h * verticalShift + gapSize / 2,
            w * selectLineRatio,
            h * obj['selectedCount'] / this.max - gapSize
          );
        }

        // @ts-ignore
        if (this.dataFrame.currentRow[name] === obj['name']) {
          g.strokeStyle = 'rgb(0,0,0)';
          g.strokeRect(x, y + h * verticalShift + gapSize / 2, w, subBartHeight);
        }

        sum -= obj['count'];
      });
      return;
    }

    render(computeData = true): void {
      const df = this.dataFrame!;
      if (computeData)
        this.computeData(df);

      if (this.tableCanvas) {
        for (const name of this.aminoColumnNames)
          this.renderBar(name);
      }
      return;
    }

    onPropertyChanged(property: DG.Property): void {
      super.onPropertyChanged(property);
    }

    register(args: DG.GridCellRenderArgs): void {
      this.registered[args.cell.tableColumn!.name] = args.cell;
    }

    unregister(name: string): void {
      if (this.registered[name])
        delete this.registered[name];
    }


    renderBar(name: string): void {
      if (!this.registered[name] || !this.tableCanvas)
        return;

      const cell = this.registered[name];
      const rect = cell.bounds;
      this.renderBarToCanvas(this.tableCanvas.getContext('2d')!, cell, rect.x, rect.y, rect.width, rect.height);
    }

    highlight(cell: DG.GridCell, offsetX:number, offsetY:number): void {
      if (!cell.tableColumn?.name || !this.aminoColumnNames.includes(cell.tableColumn.name))
        return;
        
      const colName = cell.tableColumn?.name;
      const margin = 0.2;
      const bound = cell.bounds;
      const x = bound.x + bound.width * margin;
      const y = 130 * margin / 8 - 3;
      const w = bound.width - bound.width * margin * 2;
      const h = 130 - 130 * margin / 2;
      const barData = this.barStats[colName];
      let sum = 0;

      barData.forEach((obj) => {
        sum += obj['count'];
      });

      barData.forEach((obj, index) => {
        const sBarHeight = h * obj['count'] / this.max;
        const start = h * (this.max - sum) / this.max;

        if (offsetX >= x &&
            offsetY >= y + start &&
            offsetX <= x + w &&
            offsetY <= y + start + sBarHeight) {
          this.highlighted = {'colName': colName, 'aaName': obj['name']};
        }
        
        sum -= obj['count'];
      });
    }

    unhighlight(): void {
      this.highlighted = null;
    }

    beginSelection(event:any): void {
      if (!this.highlighted || !this.dataFrame)
        return;

      this.dataFrame!.selection.handleClick((i) => {
        return this.highlighted!['aaName'] === (this.dataFrame!.getCol(this.highlighted!['colName']).get(i));
      }, event);
      this.render(true);
    }
}
