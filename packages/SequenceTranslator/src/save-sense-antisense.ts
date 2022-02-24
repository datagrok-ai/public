import * as ui from 'datagrok-api/ui';
import {sequenceToSmiles} from './package';
import * as OCL from 'openchemlib/full.js';

export function saveSenseAntiSense() {
  const ssInput = ui.textInput('Sense Strand 5\' ->3\'', '');
  const asInput = ui.textInput('Anti Sense 3\' ->5\'', '');
  const saveOption = ui.switchInput('save as one entity', false);
  const saveBtn = ui.button('Save SDF', () => {
    const smiSS = sequenceToSmiles(ssInput.value);
    const smiAS = sequenceToSmiles(asInput.value, true);
    let result: string;
    if (saveOption.value)
      result = `${OCL.Molecule.fromSmiles(smiSS + '.' + smiAS).toMolfile()}\n\n$$$$\n`;
    else {
      result =
      `${OCL.Molecule.fromSmiles(smiSS).toMolfile()}\n` +
      `>  <Sequence>\nSense Strand\n\n$$$$\n` +
      `${OCL.Molecule.fromSmiles(smiAS).toMolfile()}\n` +
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
          ui.buttonsInput([saveBtn]),
        ], 'ui-form'),
      ], 'ui-form'),
    ], 'ui-form'),
  ]);

  return saveSection;
}
