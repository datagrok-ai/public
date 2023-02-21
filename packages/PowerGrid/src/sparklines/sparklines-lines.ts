import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {
  getSettingsBase,
  names,
  SparklineType,
  SummarySettingsBase,
  Hit,
  distance,
  createTooltip
} from './shared';

import {fit, sigmoid, FitErrorModel} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';

const minDistance = 5;

// interface for getPos function data
interface getPosConstants {
  b: DG.Rect;
  settings: SparklineSettings;
  cols: DG.Column[];
}

function scaleCoordinates(coordinates: {x: number[], y: number[]}, canvasWidth: number, canvasHeight: number) {
  const xMin = Math.min(...coordinates.x);
  const xMax = Math.max(...coordinates.x);
  const yMin = Math.min(...coordinates.y);
  const yMax = Math.max(...coordinates.y);

  const xDiff = xMax - xMin;
  const yDiff = yMax - yMin;
  const xMid = (xMax + xMin) / 2;
  const yMid = (yMax + yMin) / 2;

  let coeff: number;
  if (xDiff * canvasHeight >= yDiff * canvasWidth)
    coeff = canvasWidth / xDiff;
  else
    coeff = canvasHeight / yDiff;

  const canvasXCenter = canvasWidth / 2;
  const canvasYCenter = canvasHeight / 2;

  const newCoordinates: {x: number[], y: number[]} = {x: [], y: []};
  for (let i = 0; i < coordinates.x.length; i++) {
    newCoordinates.x[i] = canvasXCenter + coeff * (coordinates.x[i] - xMid);
    newCoordinates.y[i] = canvasYCenter + coeff * (coordinates.y[i] - yMid);
  }

  return newCoordinates;
}


function getPos(col: number, row: number, constants: getPosConstants): DG.Point {
  const b = constants.b;
  const settings = constants.settings;
  const cols = constants.cols;
  const gmin = settings.globalScale ? Math.min(...cols.map((c: DG.Column) => c.min)) : 0;
  const gmax = settings.globalScale ? Math.max(...cols.map((c: DG.Column) => c.max)) : 0;
  const r: number = settings.globalScale ? (cols[col].get(row) - gmin) / (gmax - gmin) : cols[col].scale(row);
  return new DG.Point(
    b.left + b.width * (cols.length == 1 ? 0 : col / (cols.length - 1)),
    (b.top + b.height) - b.height * r);
}

interface SparklineSettings extends SummarySettingsBase {
  globalScale: boolean;
  colorCode: boolean;
}

function getSettings(gc: DG.GridColumn): SparklineSettings {
  gc.settings ??= getSettingsBase(gc);
  gc.settings.globalScale ??= false;
  gc.settings.colorCode ??= true;
  return gc.settings;
}

function onHit(gridCell: DG.GridCell, e: MouseEvent): Hit {
  const df = gridCell.grid.dataFrame;

  if (gridCell.bounds.width < 20 || gridCell.bounds.height < 10 || df === void 0) {
    return {
      isHit: false,
      row: -1,
      cols: [],
      activeColumn: -1
    };
  }

  const row = gridCell.cell.row.idx;
  const settings = getSettings(gridCell.gridColumn);
  const b = new DG.Rect(gridCell.bounds.x, gridCell.bounds.y, gridCell.bounds.width,
    gridCell.bounds.height).inflate(-3, -2);

  const cols = df.columns.byNames(settings.columnNames);
  const getPosConstants: getPosConstants = {
    b: b,
    settings: settings,
    cols: cols
  };


  const MousePoint = new DG.Point(e.offsetX, e.offsetY);
  const activeColumn = Math.floor((MousePoint.x - b.left + Math.sqrt(minDistance)) /
    b.width * (cols.length - 1 > 0 ? cols.length - 1 : 1));

  const activePoint = getPos(activeColumn, row, getPosConstants);
  return {
    isHit: distance(activePoint, MousePoint) < minDistance,
    activeColumn: activeColumn,
    row: row,
    cols: cols,
  };
}

export class SparklineCellRenderer extends DG.GridCellRenderer {
  get name() { return SparklineType.Sparkline; }

  get cellType() { return SparklineType.Sparkline; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const hitData = onHit(gridCell, e);
    if (hitData.isHit)
      ui.tooltip.show(ui.divV(createTooltip(hitData.cols, hitData.activeColumn, hitData.row)), e.x + 16, e.y + 16);
    else
      ui.tooltip.hide();
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 20 || h < 10 || df === void 0) return;

    const settings = getSettings(gridCell.gridColumn);
    const b = new DG.Rect(x, y, w, h).inflate(-3, -2);
    g.strokeStyle = 'lightgrey';
    g.lineWidth = 1;


    const row = gridCell.cell.row.idx;
    const cols = df.columns.byNames(settings.columnNames);

    const getPosConstants: getPosConstants = {
      b: b,
      settings: settings,
      cols: cols
    };

    g.beginPath();
    let started = false;
    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        const p = getPos(i, row, getPosConstants);

        if (!started) {
          g.moveTo(p.x, p.y);
          started = true;
        } else {
          g.lineTo(p.x, p.y);
        }
      }
    }
    g.stroke();

    for (let i = 0; i < cols.length; i++) {
      if (!cols[i].isNone(row)) {
        const p = getPos(i, row, getPosConstants);
        const color = settings.colorCode ? DG.Color.getCategoricalColor(i) : DG.Color.blue;
        DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, p.x, p.y, color, 3);
      }
    }
  }

  renderSettings(gridColumn: DG.GridColumn): HTMLElement {
    gridColumn.settings ??= {globalScale: true};
    const settings: SparklineSettings = gridColumn.settings;

    const globalScaleProp = DG.Property.js('globalScale', DG.TYPE.BOOL, {
      description: 'Determines the way a value is mapped to the vertical scale.\n' +
        '- Global Scale OFF: bottom is column minimum, top is column maximum. Use when columns ' +
        'contain values in different units.\n' +
        '- Global Scale ON: uses the same scale. This lets you compare values ' +
        'across columns, if units are the same (for instance, use it for tracking change over time).'
    });

    const normalizeInput = DG.InputBase.forProperty(globalScaleProp, settings);
    normalizeInput.onChanged(() => gridColumn.grid.invalidate());

    const colorCodeScaleProp = DG.Property.js('colorCode', DG.TYPE.BOOL, {
      description: 'Activates color rendering'
    });

    const colorCodeNormalizeInput = DG.InputBase.forProperty(colorCodeScaleProp, settings);
    colorCodeNormalizeInput.onChanged(() => { gridColumn.grid.invalidate(); });

    return ui.inputs([
      normalizeInput,
      ui.columnsInput('Сolumns', gridColumn.grid.dataFrame, (columns) => {
        settings.columnNames = names(columns);
        gridColumn.grid.invalidate();
      }, {
        available: names(gridColumn.grid.dataFrame.columns.numerical),
        checked: settings?.columnNames ?? names(gridColumn.grid.dataFrame.columns.numerical),
      }),
      colorCodeNormalizeInput,
    ]);
  }
}

export class FitCellRenderer extends DG.GridCellRenderer {
  get name() { return 'fit'; }

  get cellType() { return 'fit'; }

  // onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
  //   const hitData = onHit(gridCell, e);
  //   if (hitData.isHit)
  //     ui.tooltip.show(ui.divV(createTooltip(hitData.cols, hitData.activeColumn, hitData.row)), e.x + 16, e.y + 16);
  //   else
  //     ui.tooltip.hide();
  // }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    console.log(`click coordinates: x = ${e.offsetX}, y = ${e.offsetY}`);

    const df = gridCell.grid.dataFrame;

    const column = df.columns.toList().filter((c) => c.type === DG.COLUMN_TYPE.DATA_FRAME)![0];
    const currentDfValue = DG.toJs(column.get(gridCell.cell.row.idx)) as DG.DataFrame;
    if (!currentDfValue) return;

    const columns = currentDfValue.columns.byNames(['x', 'y']);
    const coordinates: {x: number[], y: number[]} = {x: [], y: []};
    coordinates.x = columns[0].toList();
    coordinates.y = columns[1].toList();

    // this to call rerendering??
    const updatedWidth = 160;
    gridCell.gridColumn.width = updatedWidth;

    const canvasCoordinates = scaleCoordinates(coordinates, gridCell.bounds.width, gridCell.bounds.height);

    const markerPxSize = 6;
    for (let i = 0; i < canvasCoordinates.x.length; i++) {
      const markerCenterX = gridCell.bounds.x + canvasCoordinates.x[i];
      const markerCenterY = (gridCell.bounds.y + gridCell.bounds.height) - canvasCoordinates.y[i];
      if (e.offsetX > markerCenterX - (markerPxSize / 2) && e.offsetX < markerCenterX + (markerPxSize / 2) &&
        e.offsetY > markerCenterY - (markerPxSize / 2) && e.offsetY < markerCenterY + (markerPxSize / 2))
        currentDfValue.filter.set(i, !(currentDfValue.filter.get(i)));
    }
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const df = gridCell.grid.dataFrame;

    if (w < 20 || h < 10 || df === void 0) return;

    const column = df.columns.toList().filter((c) => c.type === DG.COLUMN_TYPE.DATA_FRAME)![0];
    const currentDfValue = DG.toJs(column.get(gridCell.cell.row.idx)) as DG.DataFrame;
    const filteredIndexes = currentDfValue.filter.getSelectedIndexes();
    if (!currentDfValue) return;

    const columns = currentDfValue.columns.byNames(['x', 'y']);
    const coordinates: {x: number[], y: number[]} = {x: [], y: []};
    coordinates.x = columns[0].toList();
    coordinates.y = columns[1].toList();

    const updatedWidth = 160;
    gridCell.gridColumn.width = updatedWidth;

    const canvasPointCoordinates = scaleCoordinates(coordinates, updatedWidth, h);

    const filteredDf = currentDfValue
      .groupBy(['x', 'y'])
      .whereRowMask(currentDfValue.filter)
      .aggregate();

    const filteredColumns = filteredDf.columns.byNames(['x', 'y']);
    const filteredCoordinates: {x: number[], y: number[]} = {x: [], y: []};
    filteredCoordinates.x = filteredColumns[0].toList();
    filteredCoordinates.y = filteredColumns[1].toList();

    const minY = filteredColumns[1].min;
    const maxY = filteredColumns[1].max;
    const medY = filteredColumns[1].stats.med;
    let xAtMedY = -1;
    for (let i = 0; i < filteredCoordinates.y.length; i++) {
      if (filteredCoordinates.y[i] === medY) {
        xAtMedY = filteredCoordinates.x[i];
        break;
      }
    }


    // randomize dataframe columns with x and y (grok.data.demo.fit()?)
    // set x like this higher (fixed)
    // set y randomized between [-1; 1]
    // make all markers the same color (like in scatterplot)
    // fix max-min range, make width for renderer cell ??
    // set bigger size for renderer (width height of cell ~ 200x100 and a little demo) !!!!


    // make outliers switch (below)!
    // TODO: add filtering bitset for points (if clicked for example)
    // on click turn off the point (change bitset on false) - change it to cross (DG.MARKER_TYPE??)
    // second click returns it

    const params = [maxY, 1.2, xAtMedY, minY];
    const fitRes = fit(filteredCoordinates, params, sigmoid, FitErrorModel.Constant);

    const fitCoordinates: {x: number[], y: number[]} = {x: [], y: []};
    for (let i = 0; i < filteredCoordinates.x.length; i++) {
      fitCoordinates.x[i] = filteredCoordinates.x[i];
      fitCoordinates.y[i] = fitRes.fittedCurve(fitCoordinates.x[i]);
    }

    const canvasFitCoordinates = scaleCoordinates(fitCoordinates, updatedWidth, h);

    g.strokeStyle = 'black';

    // let index = 0;
    g.beginPath();
    g.moveTo(x + canvasFitCoordinates.x[0], (y + h) - canvasFitCoordinates.y[0]);
    for (let i = 1; i < canvasFitCoordinates.x.length; i++)
      g.lineTo(x + canvasFitCoordinates.x[i], (y + h) - canvasFitCoordinates.y[i]);
    g.stroke();

    for (let i = 0; i < canvasPointCoordinates.x.length; i++) {
      const color = DG.Color.scatterPlotMarker;
      if (filteredIndexes.includes(i)) {
        DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, x + canvasPointCoordinates.x[i],
          (y + h) - canvasPointCoordinates.y[i], color, 4);
      } else {
        DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER, x + canvasPointCoordinates.x[i],
          (y + h) - canvasPointCoordinates.y[i], color, 6);
      }
    }
  }

  // renderSettings(gridColumn: DG.GridColumn): HTMLElement {
  //   gridColumn.settings ??= {globalScale: true};
  //   const settings: SparklineSettings = gridColumn.settings;

  //   const globalScaleProp = DG.Property.js('globalScale', DG.TYPE.BOOL, {
  //     description: 'Determines the way a value is mapped to the vertical scale.\n' +
  //       '- Global Scale OFF: bottom is column minimum, top is column maximum. Use when columns ' +
  //       'contain values in different units.\n' +
  //       '- Global Scale ON: uses the same scale. This lets you compare values ' +
  //       'across columns, if units are the same (for instance, use it for tracking change over time).'
  //   });

  //   const normalizeInput = DG.InputBase.forProperty(globalScaleProp, settings);
  //   normalizeInput.onChanged(() => gridColumn.grid.invalidate());

  //   const colorCodeScaleProp = DG.Property.js('colorCode', DG.TYPE.BOOL, {
  //     description: 'Activates color rendering'
  //   });

  //   const colorCodeNormalizeInput = DG.InputBase.forProperty(colorCodeScaleProp, settings);
  //   colorCodeNormalizeInput.onChanged(() => { gridColumn.grid.invalidate(); });

  //   return ui.inputs([
  //     normalizeInput,
  //     ui.columnsInput('Сolumns', gridColumn.grid.dataFrame, (columns) => {
  //       settings.columnNames = names(columns);
  //       gridColumn.grid.invalidate();
  //     }, {
  //       available: names(gridColumn.grid.dataFrame.columns.numerical),
  //       checked: settings?.columnNames ?? names(gridColumn.grid.dataFrame.columns.numerical),
  //     }),
  //     colorCodeNormalizeInput,
  //   ]);
  // }
}
