import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

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

    if (molAS2) {
      result += molAS2+ '\n' +
      `> <Sequence>\nAnti Sense 2\n$$$$\n`;
    }
  }
  download(ss.replace(/\s/g, '') + '.sdf', encodeURIComponent(result));
}

export function saveSenseAntiSense() {
  //// const v = grok.shell.newView('SDF');

  //let ssDirection = ui.choiceInput('SS direction', '5 prime -> 3 prime', ['5 prime -> 3 prime']);
  //let asDirection = ui.choiceInput('AS direction', '5 prime -> 3 prime', ['5 prime -> 3 prime']);
  //let as2Direction = ui.choiceInput('AS 2 direction', '5 prime -> 3 prime', ['5 prime -> 3 prime']);
  //let saveEntity = ui.boolInput('Save as one entity', true);
  //let useChiral = ui.boolInput('Use chiral', true);

  //let form2 = ui.form([ssDirection, asDirection, as2Direction, saveEntity, useChiral]);
  //form2.className = 'ui-form ui-form-wide';

  //let body = ui.divH([ui.block([form1]), form2]);
  //v.append(body);
  //$(form1).find('textarea').css('flex-grow','1');
  //$(form1).find('label').css('max-width','140px');

  // end of comment

  const moleculeSvgDiv = ui.block([]);
  const inputColHeader = ui.h1('Sequences');
  const ssInput = ui.textInput('Sense Strand', '');
  ssInput.root.style.color = 'red';
  const asInput = ui.textInput('Anti Sense', '');
  const asInput2 = ui.textInput('Anti Sense 2', '');
  const saveEntity = ui.boolInput('Save as one entity', true);
  const useChiral = ui.boolInput('Use chiral', true);

  const straight = '5 prime -> 3 prime';
  const inverse = '3 prime -> 5 prime';
  let ssInverse = false;
  let asInverse = false;
  const as2Inverse = false;

  const ssDirection = ui.choiceInput('SS direction', straight, [straight, inverse]);
  ssDirection.onChanged(() => { ssInverse = ssDirection.value == inverse; });
  const asDirection = ui.choiceInput('AS direction', straight, [straight, inverse]);
  asDirection.onChanged(() => { asInverse = asDirection.value == inverse; });
  const as2Direction = ui.choiceInput('AS 2 direction', straight, [straight, inverse]);
  as2Direction.onChanged(() => { asInverse = asDirection.value == inverse; });

  const saveButton = ui.buttonsInput([
    ui.bigButton('Save SDF', () =>
      saveSdf(
        asInput.value, ssInput.value, saveEntity.value!,
        useChiral.value!, ssInverse, asInverse, asInput2.value, as2Inverse)
    )
  ]);
  //@ts-ignore
  // the above line is recommended by Dmitry because saveButton has wrong return
  // type
  const form1 = ui.form([ssInput, asInput, asInput2, saveButton]);
  form1.className = 'ui-form ui-form-wide';
  const form2 = ui.form([ssDirection, asDirection, as2Direction, saveEntity, useChiral]);
  form2.className = 'ui-form ui-form-wide';

  const body = ui.divH([ui.block([form1]), form2]);
  $(form1).find('textarea').css('flex-grow', '1');
  $(form1).find('label').css('max-width','140px');
  // const saveSection = ui.panel([
  //   ui.div([
  //     ui.divH([
  //       // ui.divH([ui.h1('Inputs')]),
  //       ui.divV([
  //         inputColHeader,
  //         ssInput,
  //         asInput,
  //         asInput2,
  //         ui.div([changeSense], {style: {width: '40'}}),
  //         changeSense,
  //         changeAntiSense,
  //         changeAntiSense2,
  //         saveOption,
  //         chirality,
  //         ui.buttonsInput([saveBtn]),
  //       ], 'ui-form'),
  //     ], 'ui-form'),
  //   ], 'ui-form'),
  //   moleculeSvgDiv,
  // ]);

  return body;
}
