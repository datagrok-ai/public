import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkV3000} from '../structures-works/mol-transformations';
import {getFormat} from '../structures-works/sequence-codes-tools';

export async function saveSdf(
  as: string, ss: string, oneEntity: boolean, useChirality: boolean, invertSS: boolean, invertAS: boolean,
) {
  const monomersLibAddress = 'System:AppData/SequenceTranslator/helmLib.json';
  const fileExists = await grok.dapi.files.exists(monomersLibAddress);
  if (!fileExists) {
    // todo: improve behaviour in this case
    grok.shell.warning('Please, provide the file with monomers library');
    return;
  }

  const monomersLib = await grok.dapi.files.readAsText(monomersLibAddress);
  const formatAs = getFormat(as);
  const formatSs = getFormat(ss);
  const molSS = sequenceToMolV3000(ss, invertSS, false, formatSs!, monomersLib);
  const molAS = sequenceToMolV3000(as, invertAS, false, formatAs!, monomersLib!);
  let result: string;
  if (oneEntity)
    result = linkV3000([molSS, molAS], true, useChirality) + '\n$$$$\n';
  else {
    result =
    molSS + '\n' +
    `> <Sequence>\nSense Strand\n$$$$\n` +
    molAS + '\n' +
    `> <Sequence>\nAnti Sense\n$$$$\n`;
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
  const straight = "5 prime -> 3 prime";
  const inverse = "3 prime -> 5 prime";
  let ssInverse = false;
  let asInverse = false;

  const changeSense = ui.choiceInput('SS direction', straight, [straight, inverse]);
  changeSense.onChanged(() => {ssInverse = changeSense.value == inverse;});
  const changeAntiSense = ui.choiceInput('AS direction', straight, [straight, inverse]);
  changeAntiSense.onChanged(() => {asInverse = changeAntiSense.value == inverse;});

  const saveOption = ui.switchInput('Save as one entity', true);
  const chirality = ui.switchInput('Use chiral', true);
  const saveBtn = ui.button('Save SDF', () =>
    saveSdf(asInput.value, ssInput.value, saveOption.value, chirality.value, ssInverse, asInverse));

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
          chirality,
          ui.buttonsInput([saveBtn]),
        ], 'ui-form'),
      ], 'ui-form'),
    ], 'ui-form'),
    moleculeSvgDiv,
  ]);

  return saveSection;
}
