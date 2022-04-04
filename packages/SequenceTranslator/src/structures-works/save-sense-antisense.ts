import * as ui from 'datagrok-api/ui';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkV3000} from '../structures-works/mol-transformations';

export function saveSdf(as: string, ss: string, oneEntity: boolean) {
  const molSS = sequenceToMolV3000(ss);
  const molAS = sequenceToMolV3000(as, true);
  let result: string;
  if (oneEntity)
    result = linkV3000([molSS, molAS], true) + '\n\n$$$$\n';
  else {
    result =
    molSS + '\n' +
    `>  <Sequence>\nSense Strand\n\n$$$$\n` +
    molAS + '\n' +
    `>  <Sequence>\nAnti Sense\n\n$$$$\n`;
  }

  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', ss.replace(/\s/g, '') + '.sdf');
  element.click();
}

export function saveSenseAntiSense() {
  const moleculeSvgDiv = ui.block([]);
  const ssInput = ui.textInput('Sense Strand 5\' ->3\'', '');
  const asInput = ui.textInput('Anti Sense 3\' ->5\'', '');
  const saveOption = ui.switchInput('Save as one entity', true);
  const saveBtn = ui.button('Save SDF', () => saveSdf(asInput.value, ssInput.value, saveOption.value));

  const saveSection = ui.panel([
    ui.div([
      ui.div([
        ui.divH([ui.h1('Inputs')]),
        ui.divV([
          ui.div([ssInput.root]),
          ui.div([asInput.root]),
          saveOption,
          ui.buttonsInput([saveBtn]),
        ], 'ui-form'),
      ], 'ui-form'),
    ], 'ui-form'),
    moleculeSvgDiv,
  ]);

  return saveSection;
}
