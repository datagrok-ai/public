/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Import for call wasm runtime-system
import { initSolvers, simulateOneCompartmentPK, simulateTwoCompartmentPK } from '../wasm/solving-tools';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  await initSolvers();
}

//name: Simulation PKPD
//tags: app
export async function sim() {
  const dosage = ui.floatInput('Dosage, um', 1000, () => recalculate());
  const doseInterval = ui.floatInput('Dose intrval, h', 12, () => recalculate());
  const compartments = ui.choiceInput('Model ', '2 compartment PK', ['2 compartment PK', '1 compartment PK'], () => recalculate());
  const clearance = ui.floatInput('Clearance', 2, () => recalculate());
  const rateConstant = ui.floatInput('Rate', 0.3, () => recalculate());
  const centralV = ui.floatInput('Central volume', 4, () => recalculate());
  const perV = ui.floatInput('Peripheral volume', 30, () => recalculate());
  const interRate = ui.floatInput('Intercompartmental rate', 1, () => recalculate());
  const effRate = ui.floatInput('Effective rate', 0.2, () => recalculate());
  const effect = ui.floatInput('EC50', 8, () => recalculate());

  let t = await simulate(dosage.value!, doseInterval.value!, compartments.value!,
  clearance.value!, rateConstant.value!, centralV.value!,
  perV.value!, interRate.value!, effRate.value!, effect.value!);

  let lc = DG.Viewer.lineChart(t);

  let recalculate = async (): Promise<void> => {
    t = await simulate(
      dosage.value!, doseInterval.value!, compartments.value!,
      clearance.value!, rateConstant.value!, centralV.value!,
      perV.value!, interRate.value!, effRate.value!, effect.value!);

    lc.dataFrame = t;
  }

  const v = grok.shell.newView('SimPKPD', [
    ui.panel([
      ui.div([
        ui.div([
          ui.divH([ui.h1('Inputs')]),
          ui.divV([
            dosage,
            doseInterval,
            compartments,
            clearance,
            rateConstant,
            centralV,
            perV,
            interRate,
            effRate,
            effect
          ], 'ui-form'),
        ], 'ui-form'),
      ], 'ui-form'),
    ])
  ]);
  v.box = true;

  v.append(lc);
}

export async function simulate(
  dosage: number, doseInterval: number, compartments: string,
  clearance: number, rateConstant: number, centralV: number,
  perV: number, interRate: number, effRate: number, effect: number): Promise<DG.DataFrame> {
  return await grok.functions.call(
    "Simpkpd:pkpd", {
    'dosage': dosage,
    'doseInterval': doseInterval,
    'compartments': compartments,
    'clearance': clearance,
    'rateConstant': rateConstant,
    'centralV': centralV,
    'perV': perV,
    'interRate': interRate,
    'effRate': effRate,
    'effect': effect
  });
}

//name: New PKPD
//description: PK/PD simulation via WebAutoSolver
//tags: model
//input: string compartments = '2 compartment PK' {category: PK model; choices: ["1 compartment PK", "2 compartment PK"]} 
//input: double initial = 0.0 {caption: initial; category: time, hours}
//input: double final = 12.0 {caption: final; category: time, hours}
//input: double step = 0.01 {caption: step; category: time, hours}
//input: double _depotInitial = 0.0 {units: ; caption: depot; category: initial values}
//input: double _centrInitial = 0.0 {units: ; caption: centr; category: initial values}
//input: double _periInitial = 0.0 {units: ; caption: peri; category: initial values}
//input: double _effInitial = 1.0 {units: ; caption: eff; category: initial values}
//input: double _KAVal = 0.3 {units: ; caption: KA; category: parameters}
//input: double _CLVal = 2.0 {units: ; caption: CL; category: parameters}
//input: double _V2Val = 4.0 {units: ; caption: V2; category: parameters}
//input: double _QVal = 1.0 {units: ; caption: Q; category: parameters}
//input: double _V3Val = 30.0 {units: ; caption: V3; category: parameters}
//input: double _KinVal = 0.2 {units: ; caption: Kin; category: parameters}
//input: double _KoutVal = 0.2 {units: ; caption: Kout; category: parameters}
//input: double _EC50Val = 8.0 {units: ; caption: EC50; category: parameters}
//output: dataframe dfSolution {caption: Solution; viewer: Line chart(x: "t", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export async function simPKPD(compartments: string,
  initial: number, final: number, step: number,
  _depotInitial: number, _centrInitial: number, _periInitial: number, _effInitial: number, 
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, _KinVal: number, _KoutVal: number, _EC50Val: number): Promise<DG.DataFrame>
{
  const simFn = (compartments === '1 compartment PK') ? simulateOneCompartmentPK : simulateTwoCompartmentPK;

  return simFn(initial, final, step,
    _depotInitial, _centrInitial, _periInitial, _effInitial, 
    _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val);
}
