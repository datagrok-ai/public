import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {
  fit,
  sigmoid,
  FitErrorModel,
  FitResult
} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';
import {GridColumn} from "datagrok-api/dg";
import {getChartData} from "./fit-data";


export function scaleCoordinates(coordinates: {x: number[], y: number[]}, canvasWidth: number, canvasHeight: number):
  {x: number[], y: number[]} {
  const xMin = Math.min(...coordinates.x);
  const xMax = Math.max(...coordinates.x);
  const yMin = Math.min(...coordinates.y);
  const yMax = Math.max(...coordinates.y);

  const xDiff = xMax - xMin;
  const yDiff = yMax - yMin;
  const xMid = (xMax + xMin) / 2;
  const yMid = (yMax + yMin) / 2;

  const coeff = xDiff * canvasHeight >= yDiff * canvasWidth ? canvasWidth / xDiff : canvasHeight / yDiff;
  const canvasXCenter = canvasWidth / 2;
  const canvasYCenter = canvasHeight / 2;

  const scaledCoordinates: {x: number[], y: number[]} = {x: [], y: []};
  for (let i = 0; i < coordinates.x.length; i++) {
    scaledCoordinates.x[i] = canvasXCenter + coeff * (coordinates.x[i] - xMid);
    scaledCoordinates.y[i] = canvasYCenter + coeff * (coordinates.y[i] - yMid);
  }

  return scaledCoordinates;
}

export function getParams(columns: DG.Column[]): number[] {
  const coordinates: {x: number[], y: number[]} = {x: columns[0].toList(), y: columns[1].toList()};

  const minY = columns[1].min;
  const maxY = columns[1].max;
  const medY = columns[1].stats.med;
  let xAtMedY = -1;
  for (let i = 0; i < coordinates.y.length; i++) {
    if (coordinates.y[i] === medY) {
      xAtMedY = coordinates.x[i];
      break;
    }
  }

  return [maxY, 1.2, xAtMedY, minY];
}


export class FitCellRenderer extends DG.GridCellRenderer {
  get name() { return 'fit-old'; }

  get cellType() { return 'fit-old'; }

  get defaultHeight(): number { return 100; }

  get defaultWidth(): number { return 160; }

  getDefaultSize(gridColumn: GridColumn): {width?: number | null, height?: number | null} {
    return { width: 160, height: 100};
  }

  // onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
  //   const hitData = onHit(gridCell, e);
  //   if (hitData.isHit)
  //     ui.tooltip.show(ui.divV(createTooltip(hitData.cols, hitData.activeColumn, hitData.row)), e.x + 16, e.y + 16);
  //   else
  //     ui.tooltip.hide();
  // }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    const df = gridCell.grid.dataFrame;

    const dfColumn = df.columns.toList().filter((column) => column.type === DG.COLUMN_TYPE.DATA_FRAME)![0];
    const coordinateDf = DG.toJs(dfColumn.get(gridCell.cell.row.idx)) as DG.DataFrame;
    if (!coordinateDf) return;

    const coordinateColumns = coordinateDf.columns.byNames(['x', 'y']);
    const coordinates: {x: number[], y: number[]} = {x: coordinateColumns[0].toList(),
      y: coordinateColumns[1].toList()};

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
        coordinateDf.filter.set(i, !(coordinateDf.filter.get(i)));
    }

    grok.shell.o = gridCell;
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
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
    // fix demo dataset (parameteres?)
    // make viewer (on selected rows with the name??)


    // make normal resizing (width height - make with bounds) ???
    // that's easy - just redefine the defaultHeight and etc. (we can do this)

    // add confidence intervals (top and bottom)
    // ask Leonid about why confidenceTop, confidenceBottom and fittedCurve give the same value

    // make viewer using echarts?? (interactive, synchronizing with the renderer)
    // make markers and lines one color on viewer (blue - blue, green - green, etc.)

    // curves have to be like dose response curve

    const data = getChartData(gridCell);

    const df = gridCell.grid.dataFrame;
    if (w < 20 || h < 10 || df === void 0) return;

    const dfColumn = df.columns.toList().filter((column) => column.type === DG.COLUMN_TYPE.DATA_FRAME)![0];
    const coordinateDf = DG.toJs(dfColumn.get(gridCell.cell.row.idx)) as DG.DataFrame;
    if (!coordinateDf) return;

    const coordinateColumns = coordinateDf.columns.byNames(['x', 'y']);
    const coordinates: {x: number[], y: number[]} = {x: coordinateColumns[0].toList(),
      y: coordinateColumns[1].toList()};

    const updatedWidth = 160;
    gridCell.gridColumn.width = updatedWidth;

    const canvasPointCoordinates = scaleCoordinates(coordinates, gridCell.gridColumn.width, h);

    const filteredCoordinateDf = coordinateDf
      .groupBy(['x', 'y'])
      .whereRowMask(coordinateDf.filter)
      .aggregate();

    const filteredColumns = filteredCoordinateDf.columns.byNames(['x', 'y']);
    const filteredCoordinates: {x: number[], y: number[]} = {x: filteredColumns[0].toList(),
      y: filteredColumns[1].toList()};

    if (data.seriesOptions?.showFitLine) {
      const params = getParams(filteredColumns);

      // caching fit results - not sure if needed
      const fitResultsColumn = gridCell.cell.dataFrame.columns.getOrCreate(`~fit:${gridCell.gridColumn.name}`, DG.TYPE.OBJECT, dfColumn.length);
      if (fitResultsColumn.isNone(gridCell.cell.row.idx))
        fitResultsColumn.set(gridCell.cell.row.idx, fit(filteredCoordinates, params, sigmoid, FitErrorModel.Constant));
      const fitResult: FitResult = fitResultsColumn.get(gridCell.cell.row.idx);

      const fitCurveCoordinates: {x: number[], y: number[]} = {x: [], y: []};
      for (let i = 0; i < filteredCoordinates.x.length; i++) {
        fitCurveCoordinates.x[i] = filteredCoordinates.x[i];
        fitCurveCoordinates.y[i] = fitResult.fittedCurve(fitCurveCoordinates.x[i]);
        // console.group('coordinates');
        // console.log(`fitted curve: x = ${fitCurveCoordinates.x[i]}, y = ${fitCurveCoordinates.y[i]}`);
        // console.log(`confidence top: x = ${fitCurveCoordinates.x[i]},
        //   y = ${fitResult.confidenceTop(fitCurveCoordinates.x[i])}`);
        // console.log(`confidence bottom: x = ${fitCurveCoordinates.x[i]},
        //   y = ${fitResult.confidenceBottom(fitCurveCoordinates.x[i])}`);
        // console.groupEnd();
      }

      const canvasFitCurveCoordinates = scaleCoordinates(fitCurveCoordinates, gridCell.gridColumn.width, h);

      g.strokeStyle = data.seriesOptions.fitLineColor ?? 'black';

      g.beginPath();
      g.moveTo(x + canvasFitCurveCoordinates.x[0], (y + h) - canvasFitCurveCoordinates.y[0]);
      for (let i = 1; i < canvasFitCurveCoordinates.x.length; i++)
        g.lineTo(x + canvasFitCurveCoordinates.x[i], (y + h) - canvasFitCurveCoordinates.y[i]);
      g.stroke();
    }

    if (data.seriesOptions?.showPoints) {
      const filteredIndexes = coordinateDf.filter.getSelectedIndexes();
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

    // const confidenceTopCoordinates: {x: number[], y: number[]} = {x: [], y: []};
    // for (let i = 0; i < filteredCoordinates.x.length; i++) {
    //   confidenceTopCoordinates.x[i] = filteredCoordinates.x[i];
    //   confidenceTopCoordinates.y[i] = fitResult.confidenceTop(confidenceTopCoordinates.x[i]);
    // }

    // const canvasConfidenceTopCoordinates = scaleCoordinates(confidenceTopCoordinates, gridCell.gridColumn.width, h);

    // g.beginPath();
    // g.moveTo(x + canvasConfidenceTopCoordinates.x[0], (y + h) - canvasConfidenceTopCoordinates.y[0]);
    // for (let i = 1; i < canvasConfidenceTopCoordinates.x.length; i++)
    //   g.lineTo(x + canvasConfidenceTopCoordinates.x[i], (y + h) - canvasConfidenceTopCoordinates.y[i]);
    // g.stroke();


    // const confidenceBottomCoordinates: {x: number[], y: number[]} = {x: [], y: []};
    // for (let i = 0; i < filteredCoordinates.x.length; i++) {
    //   confidenceBottomCoordinates.x[i] = filteredCoordinates.x[i];
    //   confidenceBottomCoordinates.y[i] = fitResult.confidenceBottom(confidenceBottomCoordinates.x[i]);
    // }

    // const canvasConfidenceBottomCoordinates = scaleCoordinates(confidenceBottomCoordinates,
    //   gridCell.gridColumn.width, h);

    // g.beginPath();
    // g.moveTo(x + canvasConfidenceBottomCoordinates.x[0], (y + h) - canvasConfidenceBottomCoordinates.y[0]);
    // for (let i = 1; i < canvasConfidenceBottomCoordinates.x.length; i++)
    //   g.lineTo(x + canvasConfidenceBottomCoordinates.x[i], (y + h) - canvasConfidenceBottomCoordinates.y[i]);
    // g.stroke();
  }
}
