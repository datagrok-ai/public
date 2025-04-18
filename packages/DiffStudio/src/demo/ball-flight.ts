/** Ball flight model for the Parameter Optimization & Sensitivity Analysis tutorials */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solveDefault} from '../solver-tools';

/** Ball flight simulation tool */
function getBallFlightTable(t1: number, velocity: number, angle: number, dB: number, roB: number): DG.DataFrame {
  const t0 = 0;
  const h = 0.01;
  const xInit = 0;
  const yInit = 0;
  const vxInit = velocity * Math.cos(angle);
  const vyInit = velocity * Math.sin(angle);
  const g = 9.81;
  const cD = 0.47;
  const roA = 1.225;

  // the problem definition
  const odes = {
    name: 'Ball flight',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [xInit, yInit, vxInit, vyInit],
    func: (t: number, _y: Float64Array, _output) => {
      // extract function values
      const vx = _y[2];
      const vy = _y[3];

      // evaluate expressions
      const v = Math.PI * dB ** 3 / 6;
      const mB = roB * v;
      const aB = 0.25 * Math.PI * dB ** 2;
      const vxSqr = vx * vx;
      const vySqr = vy * vy;
      const vSqr = vxSqr + vySqr;
      const vTotal = Math.sqrt(vSqr);
      const cosAlpha = (vTotal > 1.0e-12) ? (vx / vTotal) : 1;
      const cosBeta = (vTotal > 1.0e-12) ? (vy / vTotal) : 0;
      const drag = 0.5 * cD * roA * vSqr * aB;
      const dragX = -drag * cosAlpha;
      const dragY = -drag * cosBeta;

      // compute output
      _output[0] = vx;
      _output[1] = vy;
      _output[2] = dragX / mB;
      _output[3] = -g + dragY / mB;
    },
    tolerance: 0.00005,
    solutionColNames: ['x', 'y', 'vx', 'vy'],
  };

  return solveDefault(odes);
} // getBallFlightTable

/** Clip numeric table by min value of the specified column */
function clipNumericTable(df: DG.DataFrame, colName: string, minVal: number): DG.DataFrame {
  const col = df.col(colName);

  if (col === null)
    throw new Error('Incorrect column name');

  if (col.stats.min >= minVal)
    return df;

  const idx = (col.getRawData() as Float64Array).findIndex((val) => val < minVal);

  const cols: DG.Column[] = [];

  for (const curCol of df.columns) {
    const raw = curCol.getRawData() as Float64Array;
    cols.push(DG.Column.fromFloat64Array(
      curCol.name,
      raw.slice(0, idx),
    ));
  }

  return DG.DataFrame.fromColumns(cols);
}

/** Return dataframe with ball flight simulation */
export function getBallFlightSim(v: number, a: number, dB: number, roB: number): DG.DataFrame {
  let final = 50;
  let solution = getBallFlightTable(final, v, a, dB, roB);

  if (solution.col('y').stats.min > 0) {
    final *= 100;
    solution = getBallFlightTable(final, v, a, dB, roB);
  }
  solution = clipNumericTable(solution, 'y', 0);

  const dist = solution.col('x');
  const height = solution.col('y');

  dist.name = 'Distance';
  height.name = 'Height';

  const result = DG.DataFrame.fromColumns([dist, height]);
  result.name = 'Ball flight';

  return result;
}
