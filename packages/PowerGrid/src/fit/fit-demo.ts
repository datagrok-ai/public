import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as fit from './fit-data';
import * as fitMath from "@datagrok-libraries/statistics/src/parameter-estimation/fit-curve";
import {FIT_CELL_TYPE, IFitChartData} from "./fit-data";
import wu from "wu";

function rnd(min: number, max: number) {
  return Math.random() * (max - min) + min;
}

function sigmoid(x: number) {
  return 1 / (1 + Math.exp(-x));
}

function createSigmoidPoints(length: number, step: number):
    { x: Float32Array, y: Float32Array, params: number[] }
{
  const x = new Float32Array(length);
  const y = new Float32Array(length);

  const start = 0 - Math.floor(length / 2) * step;
  const end = start + (step * (length - 1));
  // random curve parameters; fitting will reconstruct them
  const params = [rnd(1, 2), rnd(1, 2), rnd(1, 2), rnd(-5, 5)];

  for (let num = start, i = 0; num <= end; num += step, i++) {
    x[i] = num;
    //y[i] = sigmoid(x[i]);
    y[i] = fitMath.sigmoid(params, x[i]);
  }

  // adding 20% noise
  const range = 0.2 * (Math.max(...y) - Math.min(...y));
  const dx = rnd(0, 4) - 2;
  for (let i = 0; i < length; i++) {
    x[i] += dx;
    y[i] = y[i] + rnd(0, range) - range / 2;
  }

  return {x: x, y: y, params: params};
}

function createDemoDataFrame(rowCount: number, chartsCount: number, chartsPerCell: number) {
  let df = DG.DataFrame.create(rowCount);
  const seriesLength = 15;
  const step = 0.5;

  const dataFrameColumn = df.columns.addNew('dataframes', DG.COLUMN_TYPE.DATA_FRAME); // charts as tables
  for (let i = 0; i < rowCount; i++) {
    const points = createSigmoidPoints(seriesLength, step);

    // this column encodes data as a DataFrame
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('x', points.x),
      DG.Column.fromFloat32Array('y', points.y),
    ]);
    dataFrameColumn.set(i, df);
  }

  for (let colIdx = 0; colIdx < chartsCount; colIdx++) {
    const jsonColumn = df.columns.addNewString('json chart ' + colIdx);          // charts as json
    jsonColumn.semType = fit.FIT_SEM_TYPE;

    for (let i = 0; i < rowCount; i++) {
      const chartData: IFitChartData = {
        //chartOptions: { minX: -10, minY: -2, maxX: 10, maxY: 2},
        series: []
      }

      for (let j = 0; j < chartsPerCell; j++) {
        const points = createSigmoidPoints(seriesLength, step);
        chartData.series?.push({
          parameters: j % 2 == 0 ? points.params : undefined,
          fitLineColor: DG.Color.toHtml(DG.Color.getCategoricalColor(colIdx * chartsPerCell + j)),
          pointColor: DG.Color.toHtml(DG.Color.getCategoricalColor(colIdx * chartsPerCell + j)),
          points: wu.count().take(seriesLength)
            .map(function (i) { return { x: points.x[i], y: points.y[i]}})
            .toArray()
        });
      }
      jsonColumn.set(i, JSON.stringify(chartData));
    }
  }

  return df;
}

export async function curveFitDemo() {
  const df = createDemoDataFrame(30, 5, 2);
  const grid = grok.shell.addTableView(df).grid;
  grid.columns.add({gridColumnName: 'tables', cellType: 'fit-old'}).width = 200;
  grid.props.rowHeight = 150;
}
