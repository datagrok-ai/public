/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { initSolvers } from '../wasm/solving-tools';
import { showHelpPanel } from '../wasm/help-panel';
import { simPKPD, TIME, EFFECT, CENTR_CONC } from './pk-pd-tools';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  await initSolvers();
}

//name: PK-PD
//description: Pharmacokinetic-Pharmacodynamic (PK-PD) simulation.
//tags: model
//input: string compartments = 2 compartment PK {category: PK model; choices: ["1 compartment PK", "2 compartment PK"]}
//input: double dose = 10000.0 {units: um; caption: dose; category: Dosing} [Dosage.]
//input: int dosesCount = 10 {caption: count; category: Dosing} [Number of doses.]
//input: double doseInterval = 12 {units: h; caption: interval; category: Dosing} [Dosing interval.]
//input: double _KAVal = 0.3 {units: ; caption: rate constant; category: PK parameters}
//input: double _CLVal = 2.0 {units: ; caption: clearance; category: PK parameters}
//input: double _V2Val = 4.0 {units: ; caption: central volume; category: PK parameters} [Central compartment volume.]
//input: double _QVal = 1.0 {units: ; caption: intercompartmental rate; category: PK parameters}
//input: double _V3Val = 30.0 {units: ; caption: peripheral volume; category: PK parameters} [Peripheral compartment volume.]
//input: double effRate = 0.2 {units: ; caption: effective rate; category: PD parameters} [Effective compartment rate.]
//input: double _EC50Val = 8.0 {units: ; caption: effect; category: PD parameters} [EC50.]
//output: dataframe simResults {caption: PK-PD simulation; viewer: Line chart(xColumnName: "Time [h]") | Grid(block: 50) }
//editor: Compute:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true
export async function simulatePKPD(compartments: string,
  dose: number, dosesCount: number, doseInterval: number,
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, effRate: number, _EC50Val: number): Promise<DG.DataFrame>
{
  return await simPKPD(compartments, dose, dosesCount, doseInterval, _KAVal, _CLVal, _V2Val, _QVal, _V3Val, effRate, effRate, _EC50Val);
}

//name: PKPD Demo
//description: Pharmacokinetic-Pharmacodynamic (PK-PD) simulation.
//input: string compartments = 2 compartment PK {category: PK model; choices: ["1 compartment PK", "2 compartment PK"]}
//input: double dose = 10000.0 {units: um; caption: dose; category: Dosing} [Dosage.]
//input: int dosesCount = 10 {caption: count; category: Dosing} [Number of doses.]
//input: double doseInterval = 12 {units: h; caption: interval; category: Dosing} [Dosing interval.]
//input: double _KAVal = 0.3 {units: ; caption: rate constant; category: PK parameters}
//input: double _CLVal = 2.0 {units: ; caption: clearance; category: PK parameters}
//input: double _V2Val = 4.0 {units: ; caption: central volume; category: PK parameters} [Central compartment volume.]
//input: double _QVal = 1.0 {units: ; caption: intercompartmental rate; category: PK parameters}
//input: double _V3Val = 30.0 {units: ; caption: peripheral volume; category: PK parameters} [Peripheral compartment volume.]
//input: double effRate = 0.2 {units: ; caption: effective rate; category: PD parameters} [Effective compartment rate.]
//input: double _EC50Val = 8.0 {units: ; caption: effect; category: PD parameters} [EC50.]
//output: dataframe simResults {caption: PK-PD simulation; viewer: Line chart(xColumnName: "Time [h]", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter", autoLayout: "false") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true
export async function simulatePkPdDemo(compartments: string, dose: number, dosesCount: number, doseInterval: number,
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, effRate: number, _EC50Val: number): Promise<DG.DataFrame>
{
  const res = await simPKPD(compartments, dose, dosesCount, doseInterval, _KAVal, _CLVal, _V2Val, _QVal, _V3Val, effRate, effRate, _EC50Val);

  return DG.DataFrame.fromColumns([res.col(TIME)!, res.col(CENTR_CONC)!, res.col(EFFECT)!]);
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation.
//meta.demoPath: Compute | PK-PD modeling
//test: demoSimPKPD() //wait: 100
export async function demoSimPKPD(): Promise<any>  {
  const doeSimpleFunc: DG.Func = await grok.functions.eval('SimPKPD:simulatePkPdDemo');
  const doeSimpleFuncCall = doeSimpleFunc.prepare();
    
  const openModelFunc: DG.Func = await grok.functions.eval('Compute:openModelFromFuncall');
  const openModelFuncCall = openModelFunc.prepare({'funccall': doeSimpleFuncCall});
  openModelFuncCall.call();

  showHelpPanel();
}
