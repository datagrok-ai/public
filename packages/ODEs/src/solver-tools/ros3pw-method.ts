/* The ROS3Pw method implementation

   References:
     [1] https://doi.org/10.1016/j.cam.2015.03.010 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {inverseMatrix, memAlloc, memFree} from '../../wasm/matrix-operations-api';
import {ODEs, max, abs, SAFETY, PSHRNK, PSGROW, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, EPS, tDerivative, jacobian, ERROR_MSG} from './solver-defs';

// The method specific constants
const GAMMA = 0.78867513459481287;

const GAMMA_21 = -2.3660254037844388;
const GAMMA_21_SCLAED = GAMMA / GAMMA_21;
const GAMMA_2 = GAMMA_21 + GAMMA;

const GAMMA_31 = -0.86791218280355165;
const GAMMA_31_SCLAED = GAMMA_31 / GAMMA;
const GAMMA_32 = -0.87306695894642317;
const GAMMA_32_SCLAED = GAMMA_32 / GAMMA;
const GAMMA_3 = GAMMA_31 + GAMMA_32 + GAMMA;

const ALPHA_21 = 2.3660254037844388;
const ALPHA_2 = ALPHA_21;
const ALPHA_31 = 0.5;
const ALPHA_32 = 0.76794919243112270;
const ALPHA_3 = ALPHA_31 + ALPHA_32;

const B_1 = 0.50544867840851759;
const B_2 = -0.11571687603637559;
const B_3 = 0.610268197627858;

const B_HAT_1 = 0.28973180237214197;
const B_HAT_2 = 0.10000000000000001;
const B_HAT_3 = 0.610268197627858;

const R_1 = B_1 - B_HAT_1;
const R_2 = B_2 - B_HAT_2;
const R_3 = B_3 - B_HAT_3;


/** Solve initial value problem using the ROS3Pw method [5]. */
export function ros3pw(odes: ODEs, callBack?: any): DG.DataFrame {
  /** right-hand side of the IVP solved */
  const f = odes.func;

  // operating variables
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  let h = odes.arg.step;
  const hDataframe = h;
  const tolerance = odes.tolerance;

  /** number of solution dataframe rows */
  const rowCount = Math.trunc((t1 - t0) / h) + 1;

  /** dimension of the problem */
  const dim = odes.initial.length;
  const dimSquared = dim * dim;

  /** independent variable values */
  const tArr = new Float32Array(rowCount);

  /** arrays of solution values */
  const yArrs = Array<Float32Array>(dim);

  for (let i = 0; i < dim; ++i)
    yArrs[i] = new Float32Array(rowCount);

  // method routine
  let timeDataframe = t0 + hDataframe;
  let t = t0;
  let tPrev = t0;
  let hNext = 0.0;
  let flag = true;
  let index = 1;
  let errmax = 0;
  let hTemp = 0;
  let tNew = 0;

  // 0 BUFFERS & TEMP STRUCTURES

  /** identity matrix */
  const I = new Float64Array(dim * dim);

  // compute identity matrix
  for (let i = 0; i < dim; ++i) {
    for (let j = 0; j < dim; ++j)
      I[j + i * dim] = (i === j) ? 1 : 0;
  }

  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  const wasmMemory = memAlloc(dimSquared);
  const W = new Float64Array(wasmMemory.buf, wasmMemory.off1, dimSquared);
  const invW = new Float64Array(wasmMemory.buf, wasmMemory.off1, dimSquared);

  const f0 = new Float64Array(dim);
  const k1 = new Float64Array(dim);
  const f1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const f2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const yDer = new Float64Array(dim);
  const hdT = new Float64Array(dim);

  const f0Buf = new Float64Array(dim);
  const f1Buf = new Float64Array(dim);
  const hd = 0;
  let sum = 0;
  let hDivNum = 0;

  // 1. SOLUTION AT THE POINT t0
  tArr[0] = t0;
  for (let i = 0; i < dim; ++i)
    yArrs[i][0] = y[i];

  // 2. COMPUTE NUMERICAL SOLUTION FOR THE POINTS FROM THE INTERVAL (t0, t1)
  while (flag) {
    // compute derivative
    f(t, y, dydt);

    // compute scale vector
    for (let i = 0; i < dim; ++i)
      yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;

    // check end point
    if (t + h > t1) {
      h = t1 - t;
      flag = false;
    }

    // call of adaptive step modified Rosenbrok triple method
    // computation of solution (y), time (t) and next step (hNext)
    while (true) {
      // one stage of the modified Rosenbrok triple approach
      // hdT = h * d * T(t, y, EPS);
      tDerivative(t, y, f, EPS, f0Buf, f1Buf, hdT);
      //hd = h * D / scale;
      for (let i = 0; i < dim; ++i)
        hdT[i] *= hd;

      // The main computations

      // f0 = f(t, y);
      f(t, y, f0);

      // W = I - h * d * J(t, y, EPS);
      jacobian(t, y, f, EPS, f0Buf, f1Buf, W);
      for (let i = 0; i < dimSquared; ++i)
        W[i] = I[i] - hd * W[i];

      // invW = W.inverse();
      inverseMatrix(W, dim, invW);

      // k1 = invW * (f0 + hdT);
      for (let i = 0; i < dim; ++i)
        f0Buf[i] = f0[i] + hdT[i];

      for (let i = 0; i < dim; ++i) {
        sum = 0;

        for (let j = 0; j < dim; ++j)
          sum += invW[j + i * dim] * f0Buf[j];

        k1[i] = sum;
      }

      hDivNum = 0.5 * h;

      // yDer = y + 0.5 * h * k1;
      for (let i = 0; i < dim; ++i)
        yDer[i] = y[i] + hDivNum * k1[i];

      // f1 = f(t + 0.5 * h, yDer);
      f(t + hDivNum, yDer, f1);

      // k2 = invW * (f1 - k1) + k1;
      for (let i = 0; i < dim; ++i)
        f1Buf[i] = f1[i] - k1[i];

      for (let i = 0; i < dim; ++i) {
        sum = k1[i];

        for (let j = 0; j < dim; ++j)
          sum += invW[j + i * dim] * f1Buf[j];

        k2[i] = sum;
      }

      // yOut = y + k2 * h; <--> yTemp
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * k2[i];

      // f2 = f(t + h, yOut);
      f(t + h, yTemp, f2);

      // k3 = invW * (f2 - e32 * (k2 - f1) - 2.0 * (k1 - f0) + hdT);
      for (let i = 0; i < dim; ++i)
        f1Buf[i] = f2[i] - 2 * (k2[i] - f1[i]) - 2.0 * (k1[i] - f0[i]) + hdT[i];

      for (let i = 0; i < dim; ++i) {
        sum = 0;

        for (let j = 0; j < dim; ++j)
          sum += invW[j + i * dim] * f1Buf[j];

        k3[i] = sum;
      }

      // yErr = (k1 - 2.0 * k2 + k3) * h / 6;
      hDivNum = h / 6;

      for (let i = 0; i < dim; ++i)
        yErr[i] = (k1[i] - 2.0 * k2[i] + k3[i]) * hDivNum;

      // estimating error
      errmax = 0;
      for (let i = 0; i < dim; ++i)
        errmax = max(errmax, abs(yErr[i] / yScale[i]));
      errmax /= tolerance;

      // processing the error obtained
      if (errmax > 1) {
        hTemp = SAFETY * h * errmax**PSHRNK;
        h = max(hTemp, REDUCE_COEF * h);
        tNew = t + h;
        if (tNew == t) {
          memFree(wasmMemory.off1);
          memFree(wasmMemory.off2);
          throw new Error(ERROR_MSG.ROSENBROCK_METHOD_FAILS);
        }
      } else {
        if (errmax > ERR_CONTR)
          hNext = SAFETY * h * errmax**PSGROW;
        else
          hNext = GROW_COEF * h;
        t = t + h;

        for (let i = 0; i < dim; ++i)
          y[i] = yTemp[i];

        break;
      }
    } // while (true)

    // compute lineraly interpolated results and store them in dataframe
    while (timeDataframe < t) {
      const cLeft = (t - timeDataframe) / (t - tPrev);
      const cRight = 1.0 - cLeft;

      tArr[index] = timeDataframe;

      for (let j = 0; j < dim; ++j)
        yArrs[j][index] = cRight * y[j] + cLeft * yPrev[j];

      timeDataframe += hDataframe;
      ++index;
    }

    h = hNext;
    tPrev = t;

    for (let i = 0; i < dim; ++i)
      yPrev[i] = y[i];
  } // while (flag)

  // 3. solution at the point t1
  tArr[rowCount - 1] = t1;

  for (let i = 0; i < dim; ++i)
    yArrs[i][rowCount - 1] = y[i];

  // 4. Free wasm buffer
  memFree(wasmMemory.off1);
  memFree(wasmMemory.off2);

  const solutionDf = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array(odes.arg.name, tArr)]);
  yArrs.forEach((arr, idx) => solutionDf.columns.add(DG.Column.fromFloat32Array(odes.solutionColNames[idx], arr)));
  solutionDf.name = odes.name;

  return solutionDf;
} // solveODEs
