import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


function getRandomFloat(min: number, max: number, decimals: number) {
  return parseFloat((Math.random() * (max - min) + min).toFixed(decimals));
}

function sigmoid(x: number) {
  return 1 / (1 + Math.exp(-x));
}

function createRandomDataFrame(pointNumber: number, step: number, noise: number, decimals: number) {
  const x = new Float32Array(pointNumber);
  const y = new Float32Array(pointNumber);

  const beginNumber = 0 - Math.floor(pointNumber / 2) * step;
  const endNumber = beginNumber + (step * (pointNumber - 1));

  for (let num = beginNumber, i = 0; num <= endNumber; num += step, i++) {
    x[i] = num;
    y[i] = getRandomFloat(sigmoid(x[i]), sigmoid(x[i]) + noise, decimals);
  }

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromFloat32Array('x', x),
    DG.Column.fromFloat32Array('y', y),
  ]);

  return df;
}

function createDemoDataFrame(rowNumber: number) {
  const pointNumber = 15;
  const step = 0.4;
  const noise = 0.25;
  const decimals = 10;

  const dataFrameColumn = DG.Column.fromType(DG.COLUMN_TYPE.DATA_FRAME, 'dataframes', rowNumber);

  for (let i = 0; i < rowNumber; i++) {
    const df = createRandomDataFrame(pointNumber, step, noise, decimals);
    dataFrameColumn.set(i, df);
  }

  const table = DG.DataFrame.fromColumns([dataFrameColumn]);
  return table;
}

export async function curveFitDemo() {
  const df = createDemoDataFrame(30);
  const grid = grok.shell.addTableView(df).grid;
  grid.columns.add({cellType: 'fit'});
}
