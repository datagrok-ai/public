import * as ui from 'datagrok-api/ui';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkV3000} from '../structures-works/mol-transformations';
import {getFormat} from '../structures-works/sequence-codes-tools';

export function saveSdf(as: string, ss: string, oneEntity: boolean, invertSS: boolean, invertAS: boolean) {
  const formatAs = getFormat(as);
  const formatSs = getFormat(ss);
  const molSS = sequenceToMolV3000(ss, invertSS, false, formatSs!);
  const molAS = sequenceToMolV3000(as, invertAS, false, formatAs!);
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
  const ssInput = ui.textInput('Sense Strand', '');
  const asInput = ui.textInput('Anti Sense', '');
  const straight = '5\' ->3\'';
  const inverse = '3\' ->5\'';
  let ssInverse = false;
  let asInverse = false;

  const changeSense = ui.choiceInput('SS direction', straight, [straight, inverse]);
  changeSense.onChanged(() =>{ ssInverse = changeSense.value == inverse; });
  const changeAntiSense = ui.choiceInput('AS direction', straight, [straight, inverse]);
  changeAntiSense.onChanged(() =>{ asInverse = changeAntiSense.value == inverse; });

  const saveOption = ui.switchInput('Save as one entity', true);
  const saveBtn = ui.button('Save SDF', () =>
    saveSdf(asInput.value, ssInput.value, saveOption.value, ssInverse, asInverse));

  const saveSection = ui.panel([
    ui.div([
      ui.div([
        ui.divH([ui.h1('Inputs')]),
        ui.divV([
          ssInput,
          asInput,
          ui.div([changeSense], {style: {width: '40'}}),
          changeSense,
          changeAntiSense,
          saveOption,
          ui.buttonsInput([saveBtn]),
        ], 'ui-form'),
      ], 'ui-form'),
    ], 'ui-form'),
    moleculeSvgDiv,
  ]);

  return saveSection;
}
