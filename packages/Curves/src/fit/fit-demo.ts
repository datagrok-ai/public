/* eslint-disable no-multi-spaces */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {IFitChartData, sigmoid} from '@datagrok-libraries/statistics/src/fit/fit-curve';

import wu from 'wu';
import {FitConstants} from "./const";

/** Returns random number from the interval */
function rnd(min: number, max: number): number {
  return Math.random() * (max - min) + min;
}

export function createSigmoidPoints(length: number, step: number, pointsPerX: number = 1):
    { x: Float32Array, y: Float32Array, params: number[] } {
  const x = new Float32Array(length * pointsPerX);
  const y = new Float32Array(length * pointsPerX);

  // random curve parameters; fitting will reconstruct them
  const params = [rnd(1, 2), rnd(1, 2), rnd(1, 2), rnd(-5, 5)];

  const start = 0 - Math.floor(length / 2) * step;
  const end = start + (step * (length - 1));
  for (let num = start, i = 0; num <= end; num += step, i++) {
    for (let j = 0; j < pointsPerX; j++) {
      x[i * pointsPerX + j] = num - start + 0.1;
      const floatParams = new Float32Array(4);
      floatParams.set(params);
      y[i * pointsPerX + j] = sigmoid(floatParams, num);
    }
  }

  // adding 20% noise
  const range = (pointsPerX === 1 ? 0.2 : 0.4) * (Math.max(...y) - Math.min(...y));
  const minY = Math.min(...y);
  for (let i = 0; i < x.length; i++) {
    y[i] = -minY + y[i] + rnd(0, range);
  }

  return {x: x, y: y, params: params};
}

export function createDemoDataFrame(rowCount: number, chartsCount: number, chartsPerCell: number): DG.DataFrame {
  const df = DG.DataFrame.create(rowCount);
  const seriesLength = 15;
  const step = 0.5;

  // const dataFrameColumn = df.columns.addNew('dataframes', DG.COLUMN_TYPE.DATA_FRAME); // charts as tables
  // for (let i = 0; i < rowCount; i++) {
  //   const points = createSigmoidPoints(seriesLength, step);

  //   // this column encodes data as a DataFrame
  //   const df = DG.DataFrame.fromColumns([
  //     DG.Column.fromFloat32Array('x', points.x),
  //     DG.Column.fromFloat32Array('y', points.y),
  //   ]);
  //   dataFrameColumn.set(i, df);
  // }

  for (let colIdx = 0; colIdx < chartsCount; colIdx++) {
    const pointsPerX = colIdx === 3 ? 5 : 1;

    const jsonColumn = df.columns.addNewString(`json chart ${colIdx}`); // charts as json
    jsonColumn.semType = FitConstants.FIT_SEM_TYPE;
    let charts = colIdx % 2 === 0 ? chartsPerCell : 1;

    for (let i = 0; i < rowCount; i++) {
      const chartData: IFitChartData = {
        //chartOptions: { minX: -10, minY: -2, maxX: 10, maxY: 2},
        series: [],
        chartOptions: {
          showStatistics: charts === 1 ? ['auc'] : [],
          xAxisName: 'Conc.',
          yAxisName: 'Activity',
          title: 'Dose-Response curves'
        }
      };

      for (let j = 0; j < charts; j++) {
        const points = createSigmoidPoints(seriesLength, step, pointsPerX);
        let color = DG.Color.toHtml(DG.Color.getCategoricalColor(colIdx * chartsPerCell + j));
        chartData.series?.push({
          parameters: undefined,
          // TODO: make better parameter generating
          // parameters: j % 2 === 0 ? points.params : undefined,
          fitLineColor: color,
          pointColor: color,
          showCurveConfidenceInterval: charts === 1,
          points: wu.count().take(seriesLength * pointsPerX)
            .map(function(i) { return {x: points.x[i], y: points.y[i]}; })
            .toArray()
        });
      }
      jsonColumn.set(i, JSON.stringify(chartData));
    }
  }

  return df;
}

export async function curveDemo(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  // const df = createDemoDataFrame(30, 5, 2);
  const df = await grok.data.loadTable(`${_package.webRoot}files/curves-demo.csv`);
  const tableView = grok.shell.addTableView(df);
  // tableView.addViewer('MultiCurveViewer');
}
