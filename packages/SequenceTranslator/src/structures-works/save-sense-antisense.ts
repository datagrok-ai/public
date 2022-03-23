import * as ui from 'datagrok-api/ui';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkV3000} from '../structures-works/mol-transformations';

export function saveSenseAntiSense() {
  const moleculeSvgDiv = ui.block([]);
  const ssInput = ui.textInput('Sense Strand 5\' ->3\'', '');
  const asInput = ui.textInput('Anti Sense 3\' ->5\'', '');
  const saveOption = ui.switchInput('save as one entity', true);
  const save3dx = ui.switchInput('save for 3dx', true);
  const saveBtn = ui.button('Save SDF', () => {
    const molSS = sequenceToMolV3000(ssInput.value);
    const molAS = sequenceToMolV3000(asInput.value, true);
    let result: string;
    if (saveOption.value)
      result = linkV3000([molSS, molAS], save3dx.value, true) + '\n\n$$$$\n';
    else {
      result =
      molSS + '\n' +
      `>  <Sequence>\nSense Strand\n\n$$$$\n` +
      molAS + '\n' +
      `>  <Sequence>\nAnti Sense\n\n$$$$\n`;
    }

    const element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', ssInput.value.replace(/\s/g, '') + '.sdf');
    element.click();
  });

  const saveSection = ui.panel([
    ui.div([
      ui.div([
        ui.divH([ui.h1('Inputs')]),
        ui.divV([
          ui.div([ssInput.root]),
          ui.div([asInput.root]),
          saveOption,
          save3dx,
          ui.buttonsInput([saveBtn]),
        ], 'ui-form'),
      ], 'ui-form'),
    ], 'ui-form'),
    moleculeSvgDiv,
  ]);

  return saveSection;
}
