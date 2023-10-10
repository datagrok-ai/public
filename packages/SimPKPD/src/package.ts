/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { initSolvers } from '../wasm/solving-tools';
import { simPKPD } from './pk-pd-tools';

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

//name: PK-PD
//description: Pharmacokinetic-Pharmacodynamic (PK-PD) simulation
//tags: model
//input: string compartments = "2 compartment PK" {category: PK model; choices: ["1 compartment PK", "2 compartment PK"]} [Pharmacokinetic model.]
//input: double dose = 10000.0 {units: um; caption: dose; category: Dosing} [Dosage.]
//input: int dosesCount = 10 {caption: count; category: Dosing} [Number of doses.]
//input: double doseInterval = 12 {units: h; caption: interval; category: Dosing} [Dosing interval.]
//input: double _KAVal = 0.3 {units: ; caption: rate constant; category: PK parameters} [Rate constant.]
//input: double _CLVal = 2.0 {units: ; caption: clearance; category: PK parameters} [Clearance.]
//input: double _V2Val = 4.0 {units: ; caption: central volume; category: PK parameters} [Central compartment volume.]
//input: double _QVal = 1.0 {units: ; caption: intercompartmental rate; category: PK parameters} [Intercompartmental rate.]
//input: double _V3Val = 30.0 {units: ; caption: peripheral volume; category: PK parameters} [Peripheral compartment volume.]
//input: double effRate = 0.2 {units: ; caption: effective rate; category: PD parameters} [Effective compartment rate.]
//input: double _EC50Val = 8.0 {units: ; caption: effect; category: PD parameters} [EC50.]
//output: dataframe simResults {caption: PK-PD simulation; viewer: Line chart(xColumnName: "Time [h]") | Grid(block: 50) }
//editor: Compute:RichFunctionViewEditor
export async function simulatePKPD(compartments: string,
  dose: number, dosesCount: number, doseInterval: number,
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, effRate: number, _EC50Val: number): Promise<DG.DataFrame>
{
  return simPKPD(compartments, dose, dosesCount, doseInterval, _KAVal, _CLVal, _V2Val, _QVal, _V3Val, effRate, effRate, _EC50Val);
}
