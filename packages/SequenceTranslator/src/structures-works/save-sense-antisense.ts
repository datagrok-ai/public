import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {download} from '../helpers';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkStrandsV3000} from '../structures-works/mol-transformations';
import {getFormat} from '../structures-works/sequence-codes-tools';

export function saveSdf(as: string, ss: string,
                        oneEntity: boolean, useChirality: boolean,
                        invertSS: boolean, invertAS: boolean,
                        as2: string | null = null, invertAS2: boolean | null) {
  const formatAs = getFormat(as);
  const formatSs = getFormat(ss);
  let formatAs2: string | null = null;
  let molAS2: string | null = null;

  const molSS = sequenceToMolV3000(ss, invertSS, false, formatSs!);
  const molAS = sequenceToMolV3000(as, invertAS, false, formatAs!);

  if (as2 != null && as2 != '') {
    formatAs2 = getFormat(as2!);
    molAS2 = sequenceToMolV3000(as2, invertAS2!, false, formatAs2!);
  }

  let result: string;
  if (oneEntity) {
    const antiStrands = molAS2 == null ? [molAS] : [molAS, molAS2];
    result = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChirality) + '\n$$$$\n';

  } else {
    result =
    molSS + '\n' +
    `> <Sequence>\nSense Strand\n$$$$\n` +
    molAS + '\n' +
    `> <Sequence>\nAnti Sense\n$$$$\n`;

    if (molAS2)
      result += molAS2+ '\n' +
      `> <Sequence>\nAnti Sense 2\n$$$$\n`;
  }
  download(ss.replace(/\s/g, '') + '.sdf', encodeURIComponent(result));
}

export function saveSenseAntiSense() {
  const moleculeSvgDiv = ui.block([]);
  const ssInput = ui.textInput('Sense Strand', '');
  const asInput = ui.textInput('Anti Sense', '');
  const asInput2 = ui.textInput('Anti Sense 2', '');
  const straight = "5 prime -> 3 prime";
  const inverse = "3 prime -> 5 prime";
  let ssInverse = false;
  let asInverse = false;
  let as2Inverse = false;

  const changeSense = ui.choiceInput('SS direction', straight, [straight, inverse]);
  changeSense.onChanged(() => {ssInverse = changeSense.value == inverse;});
  const changeAntiSense = ui.choiceInput('AS direction', straight, [straight, inverse]);
  changeAntiSense.onChanged(() => {asInverse = changeAntiSense.value == inverse;});
  const changeAntiSense2 = ui.choiceInput('AS 2 direction', straight, [straight, inverse]);
  changeAntiSense2.onChanged(() => {asInverse = changeAntiSense.value == inverse;});

  const saveOption = ui.switchInput('Save as one entity', true);
  const chirality = ui.switchInput('Use chiral', true);
  const saveBtn = ui.button('Save SDF', () =>
    saveSdf(asInput.value, ssInput.value, saveOption.value, chirality.value, ssInverse, asInverse, asInput2.value, as2Inverse));

  const saveSection = ui.panel([
    ui.div([
      ui.div([
        ui.divH([ui.h1('Inputs')]),
        ui.divV([
          ssInput,
          asInput,
          asInput2,
          ui.div([changeSense], {style: {width: '40'}}),
          changeSense,
          changeAntiSense,
          changeAntiSense2,
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
