/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ODEview } from './views/ode-view';

export const _package = new DG.Package();

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
