import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import $ from 'cash-dom';

import {download} from '../utils/helpers';

import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {extractAtomDataV3000} from '../utils/structures-works/mol-transformations';
import {getFormat} from './sequence-codes-tools';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

function getMolfileForImg(
  as: string, ss: string,
  useChirality: boolean,
  invertSS: boolean, invertAS: boolean,
  as2: string | null = null, invertAS2: boolean | null
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
  const resultingMolfile = linkStrandsV3000({senseStrands: [molSS], antiStrands: antiStrands}, useChirality);
  return resultingMolfile;
}

async function drawMolecule(
  moleculeImgDiv: HTMLDivElement,
  ss: string, as: string, as2: string | null = null,
  invertSS: boolean, invertAS: boolean, invertAS2: boolean | null,
  useChirality: boolean
): Promise<void> {
  moleculeImgDiv.innerHTML = ''; // clearing childs
  const formCanvasWidth = 500;
  const formCanvasHeight = 170;
  const formCanvas = ui.canvas(
    formCanvasWidth * window.devicePixelRatio, formCanvasHeight * window.devicePixelRatio);
  formCanvas.style.width = `${formCanvasWidth}px`;
  formCanvas.style.height = `${formCanvasHeight}px`;

  formCanvas.addEventListener('click', async () => {
    try {
      const mol = getMolfileForImg(as, ss, useChirality, invertSS, invertAS, as2, invertAS2);
      const addDiv = ui.div([], {style: {overflowX: 'scroll'}});

      // addDiv size required, but now available before dialog show()
      const coordinates = extractAtomDataV3000(mol);
      const cw: number = $(window).width() * 0.80; // addDiv.clientWidth
      const ch: number = $(window).height() * 0.70; // addDiv.clientHeight
      const molWidth: number = Math.max(...coordinates.x) - Math.min(...coordinates.x);
      const molHeight: number = Math.max(...coordinates.y) - Math.min(...coordinates.y);

      const wR: number = cw / molWidth;
      const hR: number = ch / molHeight;
      const r: number = hR; // Math.max(wR, hR);
      const dlgCanvasWidth = r * molWidth;
      const dlgCanvasHeight = r * molHeight;

      const dlgCanvas = ui.canvas(dlgCanvasWidth * window.devicePixelRatio, dlgCanvasHeight * window.devicePixelRatio);
      dlgCanvas.style.width = `${dlgCanvasWidth}px`;
      dlgCanvas.style.height = `${dlgCanvasHeight}px`;

      await grok.functions.call('Chem:canvasMol', {
        x: 0, y: 0, w: dlgCanvas.width, h: dlgCanvas.height, canvas: dlgCanvas,
        molString: mol, scaffoldMolString: '',
        options: {normalizeDepiction: false, straightenDepiction: false}
      });

      addDiv.appendChild(dlgCanvas);
      ui.dialog('Molecule')
        .add(addDiv)
        .showModal(true);
    } catch (err) {
      const errStr = errorToConsole(err);
      console.error(errStr);
    }
  });
  $(formCanvas).on('mouseover', () => $(formCanvas).css('cursor', 'zoom-in'));
  $(formCanvas).on('mouseout', () => $(formCanvas).css('cursor', 'default'));
  const mol = getMolfileForImg(as, ss, useChirality, invertSS, invertAS, as2, invertAS2);
  await grok.functions.call('Chem:canvasMol', {
    x: 0, y: 0, w: formCanvas.width, h: formCanvas.height, canvas: formCanvas,
    molString: mol, scaffoldMolString: '',
    options: {normalizeDepiction: false, straightenDepiction: false}
  });
  moleculeImgDiv.append(formCanvas);
}

export function saveSdf(as: string, ss: string,
  oneEntity: boolean, useChirality: boolean,
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

export function getSdfTab(): HTMLDivElement {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  const moleculeImgDiv = ui.block([]);

  // inputs
  const inputColHeader = ui.h1('Sequences');
  const ssInput = ui.textInput('Sense Strand', '', () => { onInput.next(); });
  ssInput.root.style.color = 'red';
  const asInput = ui.textInput('Anti Sense', '', () => { onInput.next(); });
  const asInput2 = ui.textInput('Anti Sense 2', '', () => { onInput.next(); });
  const saveEntity = ui.boolInput('Save as one entity', true);
  const useChiral = ui.boolInput('Use chiral', true);

  const straight = '5 prime -> 3 prime';
  const inverse = '3 prime -> 5 prime';
  let ssInverse = false;
  let asInverse = false;
  const as2Inverse = false;

  DG.debounce<string>(onInput, 300).subscribe(() => {
    drawMolecule(moleculeImgDiv,
      ssInput.value, asInput.value, asInput2.value, ssInverse, asInverse, as2Inverse, useChiral.value!);
    // updateMolecule(ssInput.value, asInput.value, asInput2.value);
    // onSequenceChanged(sequence);
  });

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

  const body = ui.divV([ui.divH([ui.block([form1]), form2]), moleculeImgDiv]);
  $(form1).find('textarea').css('flex-grow', '1');
  $(form1).find('label').css('max-width', '140px');

  return body;
}
