import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


function getRandomFloat(min: number, max: number, decimals: number) {
  return parseFloat((Math.random() * (max - min) + min).toFixed(decimals));
}

function sigmoid(x: number) {
  return 1 / (1 + Math.exp(-x));
}

function createRandomDataFrame(colLength: number) {
  const x = new Float32Array(colLength);
  const y = new Float32Array(colLength);

  const step = 0.5;
  const beginNumber = 0 - Math.floor(colLength / 2) * step;
  const endNumber = beginNumber + (step * (colLength - 1));
  const noise = 0.4;
  const decimals = 10;

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

function createDemoDataFrame() {
  const colLength = 15;
  const rowNumber = 30;
  const dataFrameColumn = DG.Column.fromType(DG.COLUMN_TYPE.DATA_FRAME, 'dataframes', rowNumber);

  for (let i = 0; i < rowNumber; i++) {
    const df = createRandomDataFrame(colLength);
    dataFrameColumn.set(i, df);
  }

  const table = DG.DataFrame.fromColumns([dataFrameColumn]);
  return table;
}

export async function curveFitDemo() {
  const df = createDemoDataFrame();
  const grid = grok.shell.addTableView(df).grid;
  grid.columns.add({cellType: 'fit'});
}
