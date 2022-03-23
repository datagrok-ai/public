import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import * as ui from 'datagrok-api/ui';
import {MonomerLibrary} from '../monomer-library';
import {PeptidesController} from '../peptides';

export function addViewerToHeader(grid: DG.Grid, barchart: StackedBarChart) {
  if (grid.temp['containsBarchart'])
    return;

  function eventAction(mouseMove: MouseEvent) {
    const cell = grid.hitTest(mouseMove.offsetX, mouseMove.offsetY);
    if (cell !== null && cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids')
      barchart.highlight(cell, mouseMove.offsetX, mouseMove.offsetY, mouseMove);
    else
      return;

    if (cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids')
      barchart.beginSelection(mouseMove);
    else
      barchart.unhighlight();
  }

  // The following events makes the barchart interactive
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'mousemove').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'click').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'mouseout').subscribe(() => barchart.unhighlight());

  barchart.tableCanvas = grid.canvas;

  //Setting grid options
  grid.setOptions({'colHeaderHeight': 130});

  grid.onCellTooltip((cell, x, y) => {
    if (cell.tableColumn && ['aminoAcids', 'alignedSequence'].includes(cell.tableColumn.semType) ) {
      if (!cell.isColHeader) {
        const monomerLib = cell.cell.dataFrame.temp[MonomerLibrary.id];
        PeptidesController.chemPalette.showTooltip(cell, x, y, monomerLib);
      } else {
        if (barchart.highlighted) {
          let elements: HTMLElement[] = [];
          elements = elements.concat([ui.divText(barchart.highlighted.aaName)]);
          ui.tooltip.show(ui.divV(elements), x, y);
        }
      }
    }
    return true;
  });

  grid.onCellRender.subscribe((args) => {
    const context = args.g;
    const boundX = args.bounds.x;
    const boundY = args.bounds.y;
    const boundWidth = args.bounds.width;
    const boundHeight = args.bounds.height;
    const cell = args.cell;
    context.save();
    context.beginPath();
    context.rect(boundX, boundY, boundWidth, boundHeight);
    context.clip();

    if (cell.isColHeader && barchart.aminoColumnNames.includes(cell.gridColumn.name)) {
      barchart.renderBarToCanvas(context, cell, boundX, boundY, boundWidth, boundHeight);
      args.preventDefault();
    }
    context.restore();
  });

  grid.temp['containsBarchart'] = true;
  //FIXME: for some reason barchat didn't show when running analysis. This fixes it, but it's bad. Find a way to fix
  // the problem
  barchart.unhighlight();
}

type stackedBarChartDatatype = {
  'name': string,
  'data': {'name': string, 'count': number, 'selectedCount': number, 'fixedSelectedCount': number}[],
}[];

type bartStatsType = {
  [Key: string]: {'name': string, 'count': number, 'selectedCount': number, 'fixedSelectedCount': number}[],
};

export class StackedBarChart extends DG.JsViewer {
  public dataEmptyAA: string;
  public highlighted: {'colName' : string, 'aaName' : string} | null = null;
  public tableCanvas: HTMLCanvasElement | undefined;
  public aminoColumnNames: string[] = [];
  private ord: { [Key: string]: number; } = {};
  private aminoColumnIndices: {[Key: string]: number} = {};
  private aggregatedTables: {[Key: string]: DG.DataFrame} = {};
  private aggregatedHighlightedTables: {[Key: string]: DG.DataFrame} = {};
  private max = 0;
  private barStats:
    {[Key: string]: {'name': string, 'count': number, 'highlightedCount': number, 'selectedCount': number}[]} = {};
  private selected: {'colName' : string, 'aaName' : string}[] = [];
  private highlightedMask: DG.BitSet | null = null;
  private aggregatedSelectedTables: {[Key: string]: DG.DataFrame} = {};

  constructor() {
    super();
    this.dataEmptyAA = this.string('dataEmptyAA', '-');
  }

  init() {
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
      for (const obj of value)
        this.ord[obj] = i++;
    }

    this.aminoColumnNames = [];
  }

  // Stream subscriptions
  onTableAttached() {
    this.init();
    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.computeData()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.computeData()));
      this.highlightedMask = DG.BitSet.create(this.dataFrame.rowCount);
    }
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  computeData() {
    this.aminoColumnNames = [];
    this.aminoColumnIndices = {};

    this.dataFrame!.columns.names().forEach((name: string) => {
      if (this.dataFrame!.getCol(name).semType === 'aminoAcids' &&
          !this.dataFrame!.getCol(name).categories.includes('COOH') &&
          !this.dataFrame!.getCol(name).categories.includes('NH2')) {
        this.aminoColumnIndices[name] = this.aminoColumnNames.length + 1;
        this.aminoColumnNames.push(name);
      }
    });

    this.aggregatedTables = {};
    this.aggregatedHighlightedTables = {};
    this.aggregatedSelectedTables = {};
    //TODO: optimize it, why store so many tables?
    this.aminoColumnNames.forEach((name) => {
      this.aggregatedTables[name] = this.dataFrame!
        .groupBy([name])
        .whereRowMask(this.dataFrame!.filter)
        .add('count', name, `${name}_count`)
        .aggregate();

      this.aggregatedHighlightedTables[name] = this.dataFrame!
        .groupBy([name])
        .whereRowMask(this.highlightedMask!)
        .add('count', name, `${name}_count`)
        .aggregate();

      this.aggregatedSelectedTables[name] = this.dataFrame!
        .groupBy([name])
        .whereRowMask(this.dataFrame!.selection)
        .add('count', name, `${name}_count`)
        .aggregate();
    });

    this.barStats = {};

    for (const [name, df] of Object.entries(this.aggregatedTables)) {
      const colObj: {
        'name': string,
        'data': { 'name': string, 'count': number, 'highlightedCount': number, 'selectedCount': number}[],
      } = {'name': name, 'data': []};
      const aminoCol = df.getCol(name);
      const aminoCountCol = df.getCol(`${name}_count`);
      this.barStats[colObj['name']] = colObj['data'];

      for (let i = 0; i < df.rowCount; i++) {
        const amino = aminoCol.get(i);
        const aminoCount = aminoCountCol.get(i);
        const aminoObj = {'name': amino, 'count': aminoCount, 'highlightedCount': 0, 'selectedCount': 0};
        const aggHighlightedAminoCol = this.aggregatedHighlightedTables[name].getCol(`${name}`);
        const aggHighlightedCountCol = this.aggregatedHighlightedTables[name].getCol(`${name}_count`);
        const aggSelectedAminoCol = this.aggregatedSelectedTables[name].getCol(`${name}`);
        const aggSelectedCountCol = this.aggregatedSelectedTables[name].getCol(`${name}_count`);

        if (!amino || amino === this.dataEmptyAA)
          continue;

        colObj['data'].push(aminoObj);

        for (const col of [aggHighlightedCountCol, aggSelectedCountCol]) {
          for (let j = 0; j < col.length; j++) {
            const highlightedAmino = aggHighlightedAminoCol.get(j);
            const selectedAmino = aggSelectedAminoCol.get(j);
            const curAmino = (col == aggHighlightedCountCol ? highlightedAmino : selectedAmino);
            if (curAmino == amino) {
              aminoObj[col == aggHighlightedCountCol ? 'highlightedCount' : 'selectedCount'] = col.get(j);
              break;
            }
          }
        }
      }

      colObj['data'].sort((o1, o2) => this.ord[o2['name']] - this.ord[o1['name']]);
    }

    this.max = this.dataFrame!.filter.trueCount;
  }

  renderBarToCanvas(g: CanvasRenderingContext2D, cell: DG.GridCell, x: number, y: number, w: number, h: number) {
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
    w = w - w * margin * 2;
    h = h - h * margin;
    const barWidth = w - 10;
    g.fillStyle = 'black';
    g.textBaseline = 'top';
    g.font = `${h * margin / 2}px`;
    g.fillText(name, x + (w - colNameSize) / 2, y + h + h * margin / 4);

    barData.forEach((obj) => {
      const sBarHeight = h * obj['count'] / this.max;
      const gapSize = sBarHeight * innerMargin;
      const verticalShift = (this.max - sum) / this.max;
      const [color, aarOuter] = PeptidesController.chemPalette.getColorAAPivot(obj['name']);
      const textSize = g.measureText(aarOuter);
      const fontSize = 11;
      const leftMargin = (w - (aarOuter.length > 1 ? fontSize : textSize.width - 8)) / 2;
      const subBartHeight = sBarHeight - gapSize;
      const yStart = h * verticalShift + gapSize / 2;
      const xStart = (w - barWidth) / 2;
      const absX = x + leftMargin;
      const absY = y + yStart + subBartHeight / 2 + (aarOuter.length == 1 ? + 4 : 0);
      const eps = 0.1;

      g.strokeStyle = color;
      g.fillStyle = color;
      if (textSize.width <= subBartHeight) {
        const origTransform = g.getTransform();

        if (color != PeptidesController.chemPalette.undefinedColor) {
          g.fillRect(x + xStart, y + yStart, barWidth, subBartHeight);
          g.fillStyle = 'black';
        } else
          g.strokeRect(x + xStart + 0.5, y + yStart, barWidth - 1, subBartHeight);

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
      } else
        g.fillRect(x + xStart, y + yStart, barWidth, subBartHeight);

      if (obj['selectedCount'] > eps) {
        g.fillStyle = 'rgb(255,165,0)';
        g.fillRect(
          x + xStart - w * selectLineRatio * 2,
          y + yStart,
          barWidth * selectLineRatio,
          h * obj['selectedCount'] / this.max - gapSize,
        );
      }

      if (obj['highlightedCount'] > eps && obj['highlightedCount'] > obj['selectedCount']) {
        g.fillStyle = 'rgb(209,242,251)';
        g.fillRect(
          x + xStart - w * selectLineRatio * 2,
          y + yStart + h * obj['selectedCount'] / this.max - gapSize,
          barWidth * selectLineRatio,
          h * (obj['highlightedCount'] - obj['selectedCount']) / this.max - gapSize,
        );
      }

      sum -= obj['count'];
    });
  }

  highlight(cell: DG.GridCell, offsetX:number, offsetY:number, mouseEvent: MouseEvent) {
    if (!cell.tableColumn?.name || !this.aminoColumnNames.includes(cell.tableColumn.name))
      return;

    const colName = cell.tableColumn?.name;
    const innerMargin = 0.02;
    const margin = 0.2;
    const bound = cell.bounds;
    const height = 130;
    const x = bound.x + bound.width * margin;
    const y = height * margin / 4;
    const w = bound.width - bound.width * margin * 2;
    const h = height - height * margin;
    const barData = this.barStats[colName];
    const barWidth = w - 10;
    let sum = 0;

    barData.forEach((obj) => {
      sum += obj['count'];
    });

    this.highlighted = null;
    barData.forEach((obj) => {
      const sBarHeight = h * obj['count'] / this.max;
      const gapSize = sBarHeight * innerMargin;
      const verticalShift = (this.max - sum) / this.max;
      const subBartHeight = sBarHeight - gapSize;
      const yStart = h * verticalShift + gapSize / 2;
      const xStart = (w - barWidth) / 2;

      if (offsetX >= x + xStart &&
          offsetY >= y + yStart &&
          offsetX <= x + xStart + barWidth &&
          offsetY <= y + yStart + subBartHeight)
        this.highlighted = {'colName': colName, 'aaName': obj['name']};


      sum -= obj['count'];
    });

    if (!this.highlighted)
      return;

    if (mouseEvent.type == 'click') {
      let idx = -1;

      for (let i = 0; i < this.selected.length; ++i) {
        if (JSON.stringify(this.selected[i]) == JSON.stringify(this.highlighted))
          idx = i;
      }

      if (mouseEvent.shiftKey && idx == -1)
        this.selected.push(this.highlighted);

      if (mouseEvent.shiftKey && (mouseEvent.ctrlKey || mouseEvent.metaKey) && idx != -1)
        this.selected.splice(idx, 1);
    }
  }

  unhighlight() {
    this.highlighted = null;
    this.highlightedMask!.setAll(false);
    this.computeData();
  }

  beginSelection(event: MouseEvent) {
    if (!this.dataFrame)
      return;

    this.highlightedMask!.setAll(false);

    this.dataFrame.selection.handleClick((i: number) => {
      for (const high of this.selected) {
        if (high['aaName'] === (this.dataFrame!.getCol(high['colName']).get(i)))
          return true;
      }
      return false;
    }, event);

    if (this.highlighted) {
      this.dataFrame.rows.match({[this.highlighted['colName']]: this.highlighted['aaName']}).highlight();
      this.highlightedMask!.handleClick((i: number) => {
        if (this.highlighted!['aaName'] === (this.dataFrame!.getCol(this.highlighted!['colName']).get(i)))
          return true;
        return false;
      }, event);
    }

    this.computeData();
  }
}
