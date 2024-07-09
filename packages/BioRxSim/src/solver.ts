import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {STAGE, INIT_COL_IDX} from './constants';

const MIN_COLS_COUNT = 3;

export async function solve(df: DG.DataFrame): Promise<DG.DataFrame> {
    const inpColumns = df.columns;
    const inputsCount = 19; // This is specified by the ODEs

    if (inputsCount !== df.rowCount)
      throw new Error('Incorrect inputs dataframe.');

    const operatingInputs = new Float32Array(inputsCount);
    let current: Float32Array;
    let stages: string[];

    // time vals
    let _t0: number;
    let _t1: number;
    const _h = 0.25

    // constants
    const VLinit = 7.2;
    const VTV = 10;
    const speed = 400;
    const diam = 6;
    const power = 2.1;
    const pH = 7.4;
    const k1red = 0.05604;
    const k1ox = 0.0108;
    const k2Fd = 1.35;
    const k2Fa = 110400000;
    const k2Kd = 0.04038;
    const k2Ka = 120000000;
    const k3FKa = 181200000;
    const k3FKd = 0.01188;
    const k4ox = 0.0108;
    const k4red = 0.05604;
    const ktox = 0.005;
    const krcyst = 0;
    const pO2sat = 100;
    const pKa2MEA = 8.18;
    const H = 1.072069378;
    const R = 0.082;

    const oneStage = async (_t0: number, _t1: number, inputs: Float32Array) => {
    let FFox = inputs[0];
    let KKox = inputs[1];
    let FFred = inputs[2];
    let KKred = inputs[3];
    let Ffree = inputs[4];
    let Kfree = inputs[5];
    let FKred = inputs[6];
    let FKox = inputs[7];
    let MEAthiol = inputs[8];
    let CO2 = inputs[9];
    let yO2P = inputs[10];
    let CYST = inputs[11];
    let VL = inputs[12];
    let Fin = inputs[13];
    let Fper = inputs[14];
    let qin = inputs[15];
    let yO2in = inputs[16];
    let T = inputs[17];
    let P = inputs[18];

    // the problem definition
    let odes = {
        name: 'Bioreactor',
        arg: {name: 't', start: _t0, finish: _t1, step: _h},
        initial: [FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, FKox, MEAthiol, CO2, yO2P, CYST, VL],
        func: (t: number, _y: Float32Array, _output: Float32Array) => {
          // extract function values
          const FFox = _y[0];
          const KKox = _y[1];
          const FFred = _y[2];
          const KKred = _y[3];
          const Ffree = _y[4];
          const Kfree = _y[5];
          const FKred = _y[6];
          const FKox = _y[7];
          const MEAthiol = _y[8];
          const CO2 = _y[9];
          const yO2P = _y[10];
          const CYST = _y[11];
          const VL = _y[12];
    
          // evaluate expressions
          const KF = pow(VL, -0.65) * 0.065 * pow(speed**3 * diam**5 * power / 2.16e12, 0.361);
          const MA = MEAthiol * pow(10, pH - pKa2MEA);
          const qout = qin - KF * (yO2P * H - CO2) * VL * R * T / (P * 1000);
          const OTR = KF * (yO2P * H - CO2);
          const E0 = VLinit / VL;
          const E1 = (MA * E0)**2;
          const E11 = k1red * FFox * E0 * E1;
          const E12 = k1ox * FFred * E0;
          const E31 = k2Fd * FFred * E0;
          const E32 = k2Fa * (Ffree * E0)**2 * E1;
          const E21 = k1red * KKox * E0 * E1;
          const E22 = k1ox * KKred * E0;
          const E41 = k2Kd * KKred * E0;
          const E42 = k2Ka * (Kfree * E0)**2 * E1;
          const E51 = k3FKa * Ffree * E0 * Kfree * E0;
          const E52 = k3FKd * FKred * E0;
          const E61 = k4ox * FKred * E0 * (CYST * E0)**2;
          const E62 = k4red * FKox * E0 * E1;
          const E70 = E0 * CO2;
          const E71 = (MEAthiol * E0)**2;
          const E72 = (E70 >= 0) ? sqrt(E70) : 0;
    
          // compute output
          _output[0] = -E11 + E12;
          _output[1] = -E21 + E22;
          _output[2] = E11 - E12 - E31 + E32;
          _output[3] = E21 - E22 - E41 + E42;
          _output[4] = 2 * (E31 - E32) - E51 + E52;
          _output[5] = 2 * (E41 - E42) - E51 + E52;
          _output[6] = E51 - E52 - E61 + E62;
          _output[7] = E61 - E62;
          _output[8] = 2 * (-E11 + E12 - E21 + E22 + E31 + E41 - E32 - E42 - E62 - ktox * E71 * E72)- (MEAthiol + MA) * (Fin + Fper) / VL;
          _output[9] = (Fin * pO2sat * 0.0022 - 2 * Fper * CO2) / VL + OTR - 0.5 * ktox * E71 * E72;
          _output[10] = -OTR * (VL / (VTV - VL)) * R * T * P + yO2in * qin - yO2P * qout;
          _output[11] = ktox * E71 * E72 - krcyst * CYST * E0 - (Fin + Fper) * CYST / VL;
          _output[12] = Fin - Fper;
        },
        tolerance: 0.00005,
        solutionColNames: ['FFox', 'KKox', 'FFred', 'KKred', 'Ffree', 'Kfree', 'FKred', 'FKox', 'MEAthiol', 'CO2', 'yO2P', 'CYST', 'VL']
    };
    
    let opts = {};
    
    // used Math-functions
    const pow = (x: number, y: number) => Math.pow(x, y);
    const tan = (x: number) => Math.tan(x);
    const sqrt = (x: number) => Math.sqrt(x);
    const exp = (x: number) => Math.exp(x);
    
    // used Math-constants
    const E = Math.E;
    
    // solve the problem
    const solver = await grok.functions.eval('DiffStudio:solveEquations');
    let call = solver.prepare({problem: odes, options: opts});
    await call.call();
    return call.getParamValue('df') as DG.DataFrame;
    }

    // INITIAL STAGES
    current = inpColumns.byIndex(INIT_COL_IDX).getRawData() as Float32Array;
    operatingInputs.set(current);

    _t0 = Number(inpColumns.byIndex(INIT_COL_IDX).name);
    _t1 = Number(inpColumns.byIndex(INIT_COL_IDX + 1).name);
    
    const result = await oneStage(_t0, _t1, operatingInputs);
    const solutionCols = result.columns;
    stages = new Array<string>(result.rowCount).fill('Stage 1');

    // ADD MODELING PHASES
    const addPhases = inpColumns.length - MIN_COLS_COUNT;

    for (let k = 1; k <= addPhases; ++k) {
      // Extract prev. stage point
      const lastRowIdx = result.rowCount - 1;

      let idx = -1;
      for (const col of solutionCols) {
        if (idx > -1) // We do not use the argument column
          operatingInputs[idx] = col.get(lastRowIdx);

        ++idx;
      }

      _t0 = Number(inpColumns.byIndex(INIT_COL_IDX + k).name);
      _t1 = Number(inpColumns.byIndex(INIT_COL_IDX + k + 1).name);

      // Update inputs
      current = inpColumns.byIndex(INIT_COL_IDX + k).getRawData() as Float32Array;

      for (let i = 0; i < inputsCount; ++i)
        operatingInputs[i] += current[i];
      
      // Get stage solution
      const stageSol = await oneStage(_t0, _t1, operatingInputs);

      result.append(stageSol, true);

      stages.push(...new Array<string>(stageSol.rowCount).fill(`Stage ${k + 1}`));
    }

    solutionCols.add(DG.Column.fromStrings(STAGE, stages));

    return result;
}
