import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {fit, sigmoid, FitErrorModel} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';

import {scaleCoordinates, getParams} from './fit-curve';

// import {FitCellRenderer} from './sparklines/fit-curve';


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

      // const df = this.dataFrame
      //   .groupBy(['dataframes'])
      //   .whereRowMask(this.dataFrame.selection)
      //   .aggregate();

      const selectedIndexes = this.dataFrame.selection.getSelectedIndexes();
      const dfColumn = this.dataFrame.columns.toList().filter((col) => col.type === DG.COLUMN_TYPE.DATA_FRAME)![0];

      const valueList = new Array(selectedIndexes.length) as DG.DataFrame[];
      for (let i = 0, j = 0; i < dfColumn.length; i++) {
        if (selectedIndexes.includes(i))
          valueList[j++] = DG.toJs(dfColumn.get(i));
      }

      for (let i = 0; i < valueList.length; i++) {
        const currentDf = valueList[i];
        // if (!coordinateDf) return;
        // console.log(coordinateDf.toCsv());


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
