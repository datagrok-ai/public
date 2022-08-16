/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ODEview } from './views/ode-view';

export const _package = new DG.Package();

//name: SimPKPD
//tags: app
export async function sim() {
  const dosage = ui.floatInput('Dosage, um', 1000);
  const doseInterval = ui.floatInput('Dose intrval, h', 12);
  const compartments = ui.choiceInput('Model ', '2 compartment PK', ['2 compartment PK', '1 compartment PK']);
  const clearance = ui.floatInput('Clearance', 2);
  const rateConstant = ui.floatInput('Rate', 0.3);
  const centralV = ui.floatInput('Central volume', 4);
  const perV = ui.floatInput('Peripheral volume', 30);
  const interRate = ui.floatInput('Intercompartmental rate', 1);
  const effRate = ui.floatInput('Effective rate', 0.2);
  const effect = ui.floatInput('EC50', 8);

  const moleculeSvgDiv = ui.block([]);
  
  const v = grok.shell.newView('Sequence Translator', [
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
      moleculeSvgDiv,
    ])
  ]);
  v.box = true;

  grok.functions.call(
    "Simpkpd:pkpd", {
    'dosage': dosage.value, 
    'doseInterval': doseInterval.value,
    'compartments': compartments.value,
    'clearance': clearance.value, 
    'rateConstant': rateConstant.value,
    'centralV': centralV.value, 
    'perV': perV.value,
    'interRate': interRate.value,
    'effRate': effRate.value,
    'effect': effect.value
  }).then(res => {
    moleculeSvgDiv.append(res);
  });
}
