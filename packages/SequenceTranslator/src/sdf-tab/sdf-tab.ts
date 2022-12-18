import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import $ from 'cash-dom';

import {download} from '../utils/helpers';

import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {getFormat} from './sequence-codes-tools';

import {drawMolecule} from '../utils/structures-works/draw-molecule';

function getMolfileForImg(
  ss: string, as: string,
  as2: string | null = null,
  invertSS: boolean, invertAS: boolean, invertAS2: boolean,
  useChiral: boolean
): string {
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

  const antiStrands = molAS2 == null ? [molAS] : [molAS, molAS2];
  const resultingMolfile = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChiral);
  return resultingMolfile;
}

export function saveSdf(as: string, ss: string,
  oneEntity: boolean, useChiral: boolean,
  invertSS: boolean, invertAS: boolean,
  as2: string | null = null, invertAS2: boolean | null
): void {
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
    result = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChiral) + '\n$$$$\n';
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

export function getSdfTab(): HTMLDivElement {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  const moleculeImgDiv = ui.block([]);

  // inputs
  const inputColHeader = ui.h1('Sequences');
  const ssInput = ui.textInput('Sense Strand', '', () => { onInput.next(); });
  ssInput.root.style.color = 'red';
  const asInput = ui.textInput('Anti Sense', '', () => { onInput.next(); });
  const as2Input = ui.textInput('Anti Sense 2', '', () => { onInput.next(); });
  const saveEntity = ui.boolInput('Save as one entity', true);
  const useChiralInput = ui.boolInput('Use chiral', true);

  const straight = '5 prime -> 3 prime';
  const inverse = '3 prime -> 5 prime';
  let invertSS = false;
  let invertAS = false;
  const invertAS2 = false;

  DG.debounce<string>(onInput, 300).subscribe(async () => {
    let molfile = '';
    try {
      molfile = getMolfileForImg(
        ssInput.value, asInput.value, as2Input.value, invertSS, invertAS, invertAS2, useChiralInput.value!
      );
    } finally {
      if (molfile !== '') {
        const canvasWidth = 500;
        const canvasHeight = 170;
        await drawMolecule(canvasWidth, canvasHeight, moleculeImgDiv, molfile);
      }
    }
  });

  const ssDirection = ui.choiceInput('SS direction', straight, [straight, inverse]);
  ssDirection.onChanged(() => { invertSS = ssDirection.value == inverse; });

  const asDirection = ui.choiceInput('AS direction', straight, [straight, inverse]);
  asDirection.onChanged(() => { invertAS = asDirection.value == inverse; });

  const as2Direction = ui.choiceInput('AS 2 direction', straight, [straight, inverse]);
  as2Direction.onChanged(() => { invertAS = asDirection.value == inverse; });

  const saveButton = ui.buttonsInput([
    ui.bigButton('Save SDF', () =>
      saveSdf(
        asInput.value, ssInput.value, saveEntity.value!,
        useChiralInput.value!, invertSS, invertAS, as2Input.value, invertAS2)
    )
  ]);
  //@ts-ignore
  // the above line is recommended by Dmitry because saveButton has wrong return
  // type
  const form1 = ui.form([ssInput, asInput, as2Input, saveButton]);
  form1.className = 'ui-form ui-form-wide';
  const form2 = ui.form([ssDirection, asDirection, as2Direction, saveEntity, useChiralInput]);
  form2.className = 'ui-form ui-form-wide';

  const body = ui.divV([ui.divH([ui.block([form1]), form2]), moleculeImgDiv]);
  $(form1).find('textarea').css('flex-grow', '1');
  $(form1).find('label').css('max-width', '140px');

  return body;
}
