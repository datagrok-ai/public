import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {fit, sigmoid, FitErrorModel} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';


function scaleCoordinates(coordinates: {x: number[], y: number[]}, canvasWidth: number, canvasHeight: number):
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

function getParams(columns: DG.Column[]): number[] {
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


export class FitViewer extends DG.JsViewer {
  initialized: boolean = false;
  canvas: HTMLCanvasElement | undefined = undefined;

  constructor() {
    super();
  }

  init() {
    this.initialized = true;
  }

  initChartEventListeners() {
    this.dataFrame.onSelectionChanged.subscribe((_) => {
      const ctx = this.canvas?.getContext('2d')!;
      ctx.clearRect(0, 0, this.canvas?.width!, this.canvas?.height!);

      const selectedIndexes = this.dataFrame.selection.getSelectedIndexes();
      const dfColumn = this.dataFrame.columns.toList().filter((col) => col.type === DG.COLUMN_TYPE.DATA_FRAME)![0];

      const valueList = new Array(selectedIndexes.length) as DG.DataFrame[];
      for (let i = 0, j = 0; i < dfColumn.length; i++) {
        if (selectedIndexes.includes(i))
          valueList[j++] = DG.toJs(dfColumn.get(i));
      }

      for (let i = 0; i < valueList.length; i++) {
        const currentDf = valueList[i];

        const coordinateColumns = currentDf.columns.byNames(['x', 'y']);
        const coordinates: {x: number[], y: number[]} = {x: coordinateColumns[0].toList(),
          y: coordinateColumns[1].toList()};

        const canvasPointCoordinates = scaleCoordinates(coordinates, this.canvas?.width!, this.canvas?.height!);

        const filteredCoordinateDf = currentDf
          .groupBy(['x', 'y'])
          .whereRowMask(currentDf.filter)
          .aggregate();

        const filteredColumns = filteredCoordinateDf.columns.byNames(['x', 'y']);
        const filteredCoordinates: {x: number[], y: number[]} = {x: filteredColumns[0].toList(),
          y: filteredColumns[1].toList()};

        const params = getParams(filteredColumns);
        const fitResult = fit(filteredCoordinates, params, sigmoid, FitErrorModel.Constant);

        const fitCurveCoordinates: {x: number[], y: number[]} = {x: [], y: []};
        for (let j = 0; j < filteredCoordinates.x.length; j++) {
          fitCurveCoordinates.x[j] = filteredCoordinates.x[j];
          fitCurveCoordinates.y[j] = fitResult.fittedCurve(fitCurveCoordinates.x[j]);
        }

        const canvasFitCurveCoordinates = scaleCoordinates(fitCurveCoordinates, this.canvas?.width!,
          this.canvas?.height!);

        ctx.strokeStyle = 'black';

        ctx.beginPath();
        ctx.moveTo(canvasFitCurveCoordinates.x[0], this.canvas?.height! - canvasFitCurveCoordinates.y[0]);
        for (let j = 1; j < canvasFitCurveCoordinates.x.length; j++)
          ctx.lineTo(canvasFitCurveCoordinates.x[j], this.canvas?.height! - canvasFitCurveCoordinates.y[j]);
        ctx.stroke();

        const filteredIndexes = currentDf.filter.getSelectedIndexes();
        for (let j = 0; j < canvasPointCoordinates.x.length; j++) {
          const color = DG.Color.scatterPlotMarker;
          if (filteredIndexes.includes(j)) {
            DG.Paint.marker(ctx!, DG.MARKER_TYPE.CIRCLE, canvasPointCoordinates.x[j],
              this.canvas?.height! - canvasPointCoordinates.y[j], color, 4);
          } else {
            DG.Paint.marker(ctx, DG.MARKER_TYPE.OUTLIER, canvasPointCoordinates.x[j],
              this.canvas?.height! - canvasPointCoordinates.y[j], color, 6);
          }
        }
      }
    });
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onTableAttached() {
    this.init();
    this.initChartEventListeners();

    this.render();
  }

  render() {
    this.canvas = ui.canvas(this.root.clientWidth, this.root.clientHeight);

    this.root.append(this.canvas);
  }
}
