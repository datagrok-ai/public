import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import * as ui from 'datagrok-api/ui';
import {MonomerLibrary} from '../monomer-library';

import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';

export function addViewerToHeader(grid: DG.Grid, barchart: StackedBarChart): void {
  if (grid.temp['containsBarchart'])
    return;

  const eventAction = (ev: MouseEvent): void => {
    const cell = grid.hitTest(ev.offsetX, ev.offsetY);
    if (cell?.isColHeader && cell.tableColumn?.semType == C.SEM_TYPES.AMINO_ACIDS) {
      const newBarPart = barchart.findAARandPosition(cell, ev);
      barchart._currentBarPart = newBarPart;
      barchart.requestAction(ev, newBarPart);
      barchart.computeData();
    }
  };

  // The following events makes the barchart interactive
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'mousemove').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'click').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'mouseout').subscribe(() => barchart.unhighlight());

  barchart.tableCanvas = grid.canvas;

  //Setting grid options
  grid.setOptions({'colHeaderHeight': 130});

  grid.onCellTooltip((cell, x, y) => {
    if (
      cell.tableColumn &&
      [C.SEM_TYPES.AMINO_ACIDS, C.SEM_TYPES.ALIGNED_SEQUENCE].includes(cell.tableColumn.semType as C.SEM_TYPES)
    ) {
      if (!cell.isColHeader) {
        const monomerLib = cell.cell.dataFrame.temp[MonomerLibrary.id];
        PeptidesModel.chemPalette.showTooltip(cell, x, y, monomerLib);
      } else if (barchart._currentBarPart) {
        let elements: HTMLElement[] = [];
        elements = elements.concat([ui.divText(barchart._currentBarPart.aaName)]);
        ui.tooltip.show(ui.divV(elements), x, y);
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

export class StackedBarChart extends DG.JsViewer {
  dataEmptyAA: string;
  _currentBarPart: type.BarChart.BarPart | null = null;
  tableCanvas: HTMLCanvasElement | undefined;
  aminoColumnNames: string[] = [];
  ord: { [Key: string]: number; } = {};
  aminoColumnIndices: {[Key: string]: number} = {};
  aggregatedFilterTables: type.DataFrameDict = {};
  max = 0;
  barStats: {[Key: string]: type.BarChart.BarStatsObject[]} = {};
  selected: type.BarChart.BarPart[] = [];
  aggregatedSelectedTables: type.DataFrameDict = {};
  model!: PeptidesModel;

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
      for (const obj of value)
        this.ord[obj] = i++;
    }

    this.aminoColumnNames = [];
  }

  // Stream subscriptions
  async onTableAttached(): Promise<void> {
    this.init();
    this.model = await PeptidesModel.getInstance(this.dataFrame);
    // this.controller.init(this.dataFrame);
    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.computeData()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.computeData()));
      this.subs.push(DG.debounce(this.dataFrame.onValuesChanged, 50).subscribe(() => this.computeData()));
    }
  }

  // Cancel subscriptions when the viewer is detached
  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  computeData(): void {
    this.aminoColumnNames = [];
    this.aminoColumnIndices = {};

    this.dataFrame.columns.names().forEach((name: string) => {
      if (this.dataFrame.getCol(name).semType === C.SEM_TYPES.AMINO_ACIDS &&
          !this.dataFrame.getCol(name).categories.includes('COOH') &&
          !this.dataFrame.getCol(name).categories.includes('NH2')) {
        this.aminoColumnIndices[name] = this.aminoColumnNames.length + 1;
        this.aminoColumnNames.push(name);
      }
    });

    this.aggregatedFilterTables = {};
    this.aggregatedSelectedTables = {};
    //TODO: optimize it, why store so many tables?
    this.aminoColumnNames.forEach((name) => {
      this.aggregatedFilterTables[name] = this.dataFrame
        .groupBy([name])
        .whereRowMask(this.dataFrame.filter)
        .add('count', name, `${name}_count`)
        .aggregate();

      this.aggregatedSelectedTables[name] = this.dataFrame
        .groupBy([name])
        .whereRowMask(this.dataFrame.selection)
        .add('count', name, `${name}_count`)
        .aggregate();
    });

    this.barStats = {};

    for (const [name, df] of Object.entries(this.aggregatedFilterTables)) {
      const colData: {'name': string, 'count': number, 'selectedCount': number}[] = [];
      const aminoCol = df.getCol(name);
      const aminoCountCol = df.getCol(`${name}_count`);
      this.barStats[name] = colData;

      for (let i = 0; i < df.rowCount; i++) {
        const amino = aminoCol.get(i);
        const aminoCount = aminoCountCol.get(i);
        const aminoObj = {'name': amino, 'count': aminoCount, 'selectedCount': 0};
        const aggSelectedAminoCol = this.aggregatedSelectedTables[name].getCol(`${name}`);
        const aggSelectedCountCol = this.aggregatedSelectedTables[name].getCol(`${name}_count`);

        if (!amino || amino === this.dataEmptyAA)
          continue;

        colData.push(aminoObj);

        for (let j = 0; j < aggSelectedCountCol.length; j++) {
          const selectedAmino = aggSelectedAminoCol.get(j);
          const curAmino = (selectedAmino);
          if (curAmino == amino) {
            aminoObj['selectedCount'] = aggSelectedCountCol.get(j);
            break;
          }
        }
      }

      colData.sort((o1, o2) => this.ord[o2['name']] - this.ord[o1['name']]);
    }

    this.max = this.dataFrame.filter.trueCount;
  }

  renderBarToCanvas(g: CanvasRenderingContext2D, cell: DG.GridCell, x: number, y: number, w: number, h: number): void {
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
      const [color, aarOuter] = PeptidesModel.chemPalette.getColorAAPivot(obj['name']);
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

        if (color != PeptidesModel.chemPalette.undefinedColor) {
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

      sum -= obj['count'];
    });
  }

  findAARandPosition(cell: DG.GridCell, mouseEvent: MouseEvent): {colName: string, aaName: string} | null {
    if (!cell.tableColumn?.name || !this.aminoColumnNames.includes(cell.tableColumn.name))
      return null;

    const offsetX = mouseEvent.offsetX;
    const offsetY = mouseEvent.offsetY;
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

    const xStart = x + (w - barWidth) / 2;
    for (const obj of barData) {
      const sBarHeight = h * obj['count'] / this.max;
      const gapSize = sBarHeight * innerMargin;
      const verticalShift = (this.max - sum) / this.max;
      const subBartHeight = sBarHeight - gapSize;
      const yStart = y + h * verticalShift + gapSize / 2;

      const isIntersectingX = offsetX >= xStart && offsetX <= xStart + barWidth;
      const isIntersectingY = offsetY >= yStart && offsetY <= yStart + subBartHeight;

      if (isIntersectingX && isIntersectingY)
        return {'colName': colName, 'aaName': obj['name']};

      sum -= obj['count'];
    }

    return null;
  }

  unhighlight(): void {
    ui.tooltip.hide();
    this.computeData();
  }

  /** Requests highlight/select/filter action based on currentBarPart */
  requestAction(event: MouseEvent, barPart: {colName: string, aaName: string} | null): void {
    if (!barPart)
      return;
    let aar = barPart['aaName'];
    let position = barPart['colName'];
    if (event.type === 'click') {
      event.shiftKey ? this.model.modifyCurrentSelection(aar, position) :
        this.model.initCurrentSelection(aar, position);
    } else {
      ui.tooltip.showRowGroup(this.dataFrame, (i) => {
        const currentAAR = this.dataFrame.get(position, i);
        return currentAAR === aar;
      }, event.offsetX, event.offsetY);
    }
  }
}
