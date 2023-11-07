/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {initMatrOperApi, inverseMatrix, memAlloc, memFree} from '../wasm/matrix-operations-api';

import {ODEs, solveODEs} from './solver';
import {getIVP, getScriptLines, getScriptParams, DF_NAME} from './scripting-tools';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initMatrOperApi();
}

//name: check upd
//input: int n = 3
export function checkUpd(n: number) {
  const size = n * n;

  const mem = memAlloc(size);

  const A = new Float64Array(mem.buf, mem.off1, size);
  const invA = new Float64Array(mem.buf, mem.off2, size);

  A[0] = 1; A[1] = 2; A[2] = 3;
  A[3] = 0; A[4] = 1; A[5] = 2;
  A[6] = 0; A[7] = 0; A[8] = 1;
  
  console.log('A:');
  console.log(A);

  inverseMatrix(A, n, invA);

  console.log('inv A:');
  console.log(invA);

  memFree(mem.off1);
  memFree(mem.off2);
}

//name: Example 1
//description: Test example for TypeScript ODEs solver
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time}
//input: double t1 = 10.0 {caption: Final; category: Time}
//input: double h = 0.1 {caption: Step; category: Time}
//input: double x = 0.2 {caption: x; category: Initial values; min: -5; max: 5}
//input: double y = 0.2 {caption: y; category: Initial values; min: -5; max: 5}
//input: double param1 = 1 {caption: param1; category: Parameters; min: -15; max: 15}
//input: double param2 = -1 {caption: param2; category: Parameters; min: -15; max: 15}
//output: dataframe df {caption: Example 1; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function solveExample1(t0: number, t1: number, h: number, x: number, y: number, param1: number, param2: number): DG.DataFrame {
  // constants
  const const1 = 1.0;
  const const2 = 3.0;

  return solveODEs({
    name: 'Example 1',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [x, y],    
    func: (t: number, _y: Float64Array, output: Float64Array) => {    
      // function values
      const x = _y[0];
      const y = _y[1];
    
      // expressions
      const coef1 = const1 + param1;
      const coef2 = const2 + param2 + 0.0;
    
      output[0] = coef1 * y;
      output[1] = coef2 * x;
    },
    tolerance: 0.00005,
    solutionColNames: ['x(t)', 'y(t)']
  });
}

//name: Example 2
//description: Test example for TypeScript ODEs solver
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time}
//input: double t1 = 2.0 {caption: Final; category: Time}
//input: double h = 0.01 {caption: Step; category: Time}
//input: double x = 1 {caption: x; category: Initial values; min: -5; max: 5}
//input: double y = 1 {caption: y; category: Initial values; min: -5; max: 5}
//input: double alpha = 2 {caption: alpha; category: Parameters; min: -2; max: 2}
//input: double beta = 3 {caption: beta; category: Parameters; min: -2; max: 2}
//output: dataframe df {caption: Example 2; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function solveExample2(t0: number, t1: number, h: number, x: number, y: number, alpha: number, beta: number): DG.DataFrame {
  return solveODEs({
    name: 'Example 2',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [x, y],
    func: (t: number, _y: Float64Array, output: Float64Array) => {    
      // function values
      const x = _y[0];
      const y = _y[1];  
    
      output[0] = alpha * x;
      output[1] = beta * y;
  },
    tolerance: 0.00005,
    solutionColNames: ['x(t)', 'y(t)']
  });
}

//name: Example 3
//description: Test example for TypeScript ODEs solver
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time}
//input: double t1 = 5.0 {caption: Final; category: Time}
//input: double h = 0.01 {caption: Step; category: Time}
//input: double x = 1 {caption: x; category: Initial values; min: -5; max: 5}
//input: double y = 0 {caption: y; category: Initial values; min: -5; max: 5}
//input: double alpha = 1 {caption: alpha; category: Parameters; min: -2; max: 2}
//input: double beta = 1 {caption: beta; category: Parameters; min: -2; max: 2}
//output: dataframe df {caption: Example 3; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function solveExample3(t0: number, t1: number, h: number, x: number, y: number, alpha: number, beta: number): DG.DataFrame {
  return solveODEs({
    name: 'Example 3',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [x, y],
    func: (t: number, _y: Float64Array, output: Float64Array) => {  
      // function values
      const x = _y[0];
      const y = _y[1];  
    
      output[0] = -alpha * y;
      output[1] = beta * x;
    },
    tolerance: 0.00005,
    solutionColNames: ['x(t)', 'y(t)']
  });
}

//name: Markov model
//description: Test example for TypeScript ODEs solver: the model M|M|1|1 simulation
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time}
//input: double t1 = 5.0 {caption: Final; category: Time}
//input: double h = 0.01 {caption: Step; category: Time}
//input: double p0 = 1 {caption: p0; category: Initial values; min: 0; max: 1}
//input: double p1 = 0 {caption: p1; category: Initial values; min: 0; max: 1}
//input: double lambda = 1 {caption: lambda; category: Parameters; min: 0.0001; max: 1000}
//input: double mu = 1 {caption: mu; category: Parameters; min: 0.0001; max: 1000}
//output: dataframe df {caption: M|M|1|1 simulation; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function solveMM(t0: number, t1: number, h: number, p0: number, p1: number, lambda: number, mu: number): DG.DataFrame {
  return solveODEs({
    name: 'M|M|1|1 simulation',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [p0, p1],
    func: (t: number, _y: Float64Array, output: Float64Array) => {  
      // function values
      const p0 = _y[0];
      const p1 = _y[1];
    
      // expressions
      const p2 = 1 - p0 - p1;
      const lambda00 = -lambda;
      const lambda10 = mu;
      const lambda01 = lambda;
      const lambda11 = -lambda - mu;
      const lambda21 = mu;  
    
      output[0] = lambda00 * p0 + lambda10 * p1;
      output[1] = lambda01 * p0 + lambda11 * p1 + lambda21 * p2;
    },
    tolerance: 0.00005,
    solutionColNames: ['p0(t)', 'p1(t)']
  });
}

//name: VdP
//description: Van der Pol oscillator simulation
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time}
//input: double t1 = 4200000.0 {caption: Final; category: Time}
//input: double h = 10 {caption: Step; category: Time}
//input: double x1 = -1 {caption: x1; category: Initial values}
//input: double x2 = 1 {caption: x2; category: Initial values}
//input: double mu = 1000000.0 {caption: mu; category: Parameters}
//output: dataframe df {caption: Van der Pol oscillator simulation; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function solveVdP(t0: number, t1: number, h: number, x1: number, x2: number, mu: number): DG.DataFrame {
  return solveODEs({
    name: 'Van der Pol oscillator simulation',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [x1, x2],
    func: (t: number, _y: Float64Array, output: Float64Array) => {  
      // function values
      const x1 = _y[0];
      const x2 = _y[1];
    
      output[0] = x2;
      output[1] = -x1 + mu * (1 - x1 * x1) * x2;
    },
    tolerance: 0.00005,
    solutionColNames: ['x1(t)', 'x2(t)']
  });
}

//name: Bioreactor TS
//description: Controlled fab-arm exchange mechanism simulation using TypeScript solver.
//tags: model
//input: double t0 = 0.0 {caption: Initial; category: Time, min}
//input: double t1 = 1000.0 {caption: Final; category: Time, min}
//input: double h = 0.1 {caption: Step; category: Time, min}
//input: double FFox = 0.2 {units: mmol/L; caption: FF oxidized (FFox); category: Initial values; min: 0.15; max: 0.25} 
//input: double KKox = 0.2 {units: mmol/L; caption: KK oxidized (KKox); category: Initial values; min: 0.15; max: 0.25}
//input: double FFred = 0.1 {units: mmol/L; caption: FF reduced (FFred); category: Initial values; min: 0.98; max: 1.2}
//input: double KKred = 0.1 {units: mmol/L; caption: KK reduced (KKred); category: Initial values; min: 0.98; max: 1.2}
//input: double Ffree = 0.0 {units: mmol/L; caption: F free (Ffree); category: Initial values}
//input: double Kfree = 0.0 {units: mmol/L; caption: K free (Kfree); category: Initial values}
//input: double FKred = 0.0 {units: mmol/L; caption: FK reduced (FKred); category: Initial values}
//input: double FKox = 0.0 {units: mmol/L; caption: FK oxidized (FKox); category: Initial values}
//input: double MEAthiol = 15.0 {units: mmol/L; caption: MEAthiol (MEA); category: Initial values}
//input: double CO2 = 0.12 {units: mmol/L; caption: Dissolved oxygen (CO2); category: Initial values}
//input: double yO2P = 0.209 {units: atm; caption: Atm headspace (yO2P); category: Initial values}
//input: double CYST = 0.0 {units: mmol/L; caption: Cystamine (CYST); category: Initial values}
//input: double VL = 7.2 {units: L; caption: Liquid volume (VL); category: Initial values}
//input: double qin = 1.0 {units: L/min; caption: Gas to headspace; category: Parameters}
//input: double yO2in = 0.21 {units: ; caption: Oxygen mole fraction; category: Parameters}
//input: double He = 1.3 {units: mmol/(L atm); caption: Henry's law constant; category: Parameters}
//input: double T = 300.0 {units: K; caption: System temperature; category: Parameters}
//input: double R = 0.082 {units: L atm/(mol K); caption: Gas constant; category: Parameters}
//input: double P = 1.0 {units: atm; caption: Headspace pressure; category: Parameters}
//input: double TimeToSwitch = 135.0 {units: min; caption: Switch mode time; category: Parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(block: 100, sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100)}
//editor: Compute:RichFunctionViewEditor
export function Bioreactor(t0: number, t1: number, h: number,
  FFox: number, KKox: number, FFred: number, KKred: number, Ffree: number, Kfree: number, FKred: number, FKox: number, MEAthiol: number, CO2: number, yO2P: number, CYST: number, VL: number, 
  qin: number, yO2in: number, He: number, T: number, R: number, P: number, TimeToSwitch: number): DG.DataFrame 
{
  // constants
  const VLinitial = 7.2;
  const Vtotalvessel = 10.0;
  const AgitatorSpeed = 400.0;
  const AgitatorDiameter = 6.0;
  const AgitatorPowerNumber = 2.1;
  const pH = 7.4;
  const k1red = 0.05604;
  const k1ox = 0.0108;
  const k2Fd = 1.35;
  const k2Fa = 110400000.0;
  const k2Kd = 0.04038;
  const k2Ka = 120000000.0;
  const k3FKa = 181200000.0;
  const k3FKd = 0.01188;
  const k4ox = 0.0108;
  const k4red = 0.05604;
  const kthiolox = 0.005;
  const krcyst = 0.0;
  const percentO2saturation = 100.0;
  const pKa2MEA = 8.18;

  return solveODEs({
    name: 'FAE simulation',
    arg: {name: 't', start: t0, finish: t1, step: h},
    initial: [FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, FKox, MEAthiol, CO2, yO2P, CYST, VL],
    func: (t: number, _y: Float64Array, output: Float64Array) => {   
      // function values
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
    
      // expressions
      const constForklasurface = (3.932 * Math.pow((Math.pow(AgitatorSpeed, 3.0) * Math.pow(AgitatorDiameter, 5.0) * AgitatorPowerNumber / 2160000000000), 0.361)) / 60.0;

      const klasurface = Math.pow(VL, -0.65) * constForklasurface;

      const MEAthiolate = MEAthiol * Math.pow(10.0,(pH - pKa2MEA));

      const qout = qin - klasurface*(yO2P*He - CO2) * VL * R * T / (P*1000.0);
  
      const OTR = klasurface*(yO2P*He - CO2);

      const Vg = Vtotalvessel - VL;

      const Fin = t < TimeToSwitch ? (0.0) : (0.025);

      const Fpermeate = t < TimeToSwitch ? (0.025) : (Fin);
  
      const CO2in = percentO2saturation * 7.17 / (32.0 * 100.0);

      const Vres = VLinitial / VL;

      const MEAthiolate_t_by_Vres_t_squared = Math.pow(MEAthiolate * Vres, 2.0);
  
      const FFox_to_FFred = k1red * FFox * Vres * MEAthiolate_t_by_Vres_t_squared;

      const FFred_to_FFox = k1ox * FFred * Vres;

      const FFred_to_Ffree = k2Fd * FFred * Vres;

      const Ffree_to_FFred = k2Fa * Math.pow(Ffree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared;

      const KKox_to_KKred = k1red * KKox * Vres * MEAthiolate_t_by_Vres_t_squared;

      const KKred_to_KKox = k1ox * KKred * Vres;

      const KKred_to_Kfree = k2Kd * KKred * Vres;

      const Kfree_to_KKred = k2Ka * Math.pow(Kfree * Vres, 2.0) * MEAthiolate_t_by_Vres_t_squared;

      const free_to_FKred = k3FKa * Ffree * Vres * Kfree * Vres;

      const FKred_to_free = k3FKd * FKred * Vres;

      const FKred_to_FKox = k4ox * FKred * Vres * Math.pow(CYST * Vres, 2.0);

      const FKox_to_FKred = k4red * FKox * Vres * MEAthiolate_t_by_Vres_t_squared;

      const Vres_t_by_CO2 = Vres * CO2;

      const sqrt_of_Vres_t_by_CO2 = (Vres_t_by_CO2 >= 0.0) ? Math.sqrt(Vres_t_by_CO2) : 0.0;		

      const MEAthiol_t_by_Vres_t_squared = Math.pow(MEAthiol * Vres, 2.0);

    
      output[0] = -FFox_to_FFred + FFred_to_FFox;

      output[1] = -KKox_to_KKred + KKred_to_KKox;

      output[2] = FFox_to_FFred - FFred_to_FFox - FFred_to_Ffree + Ffree_to_FFred;

      output[3] = KKox_to_KKred - KKred_to_KKox - KKred_to_Kfree + Kfree_to_KKred;

      output[4] = 2.0 * FFred_to_Ffree - 2.0 * Ffree_to_FFred - free_to_FKred + FKred_to_free;

      output[5] = 2.0 * KKred_to_Kfree - 2.0 * Kfree_to_KKred - free_to_FKred + FKred_to_free;

      output[6] = free_to_FKred - FKred_to_free - FKred_to_FKox + FKox_to_FKred;

      output[7] = FKred_to_FKox - FKox_to_FKred;

      output[8] = 2.0 * (-FFox_to_FFred + FFred_to_FFox - KKox_to_KKred
			+ KKred_to_KKox + FFred_to_Ffree + KKred_to_Kfree - Ffree_to_FFred
			- Kfree_to_KKred - FKox_to_FKred
			- kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2)
			- (MEAthiol + MEAthiolate) * (Fin + Fpermeate) / VL;

      output[9] = (Fin * CO2in - 2.0 * Fpermeate * CO2) / VL + OTR
	    - 0.5 * kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2;
			
      output[10] = -OTR * (VL / Vg) * R * T * P + yO2in * qin - yO2P * qout;

      output[11] = kthiolox * MEAthiol_t_by_Vres_t_squared * sqrt_of_Vres_t_by_CO2
		  - krcyst * CYST * Vres - (Fin + Fpermeate) * CYST / VL;

      output[12] = Fin - Fpermeate;
    },
    tolerance: 0.00005,
    solutionColNames: ['FFox(t)', 'KKox(t)', 'FFred(t)', 'KKred(t)', 'Ffree(t)', 'Kfree(t)', 'FKred(t)', 'FKox(t)', 'MEAthiol(t)', 'CO2(t)', 'yO2P(t)', 'CYST(t)', 'VL(t)']
  });  
}

//name: solve
//input: object problem
//output: dataframe df
export function solve(problem: ODEs): DG.DataFrame {  
  return solveODEs(problem); 
}

//name: EquaSleek X
//tags: app
export async function EquaSleekX() {
  const odeInput = ui.textInput('', 'Enter equations:\n');  

  const exportBtn = ui.button('export', () => {
    const scriptText = getScriptLines(getIVP(odeInput.value)).join('\n');      
    const script = DG.Script.create(scriptText);
    const sView = DG.ScriptView.create(script);
    grok.shell.addView(sView);
  });

  const solveBtn = ui.bigButton('solve', async () => {
    const ivp = getIVP(odeInput.value);
    const scriptText = getScriptLines(ivp).join('\n');    
    const script = DG.Script.create(scriptText);    
    const params = getScriptParams(ivp);    
    const call = script.prepare(params);
    await call.call();
    df = call.outputs[DF_NAME];
    view.dataFrame = call.outputs[DF_NAME];
    view.name = df.name;
    if (!viewer) {
      viewer = DG.Viewer.lineChart(df, {
        autoLayout: false,
        sharex: true, 
        multiAxis: true,
        multiAxisLegendPosition: "RightCenter",
      });
      view.dockManager.dock(viewer, 'right');
    }
  });
 
  let df = DG.DataFrame.create();
  let view = grok.shell.addTableView(df);
  let viewer: DG.Viewer | null = null;
  view.name = 'ODEs';
  let div = ui.divV([
    odeInput,
    ui.divH([exportBtn, solveBtn])
  ]);  

  view.dockManager.dock(div, 'left');  
}
