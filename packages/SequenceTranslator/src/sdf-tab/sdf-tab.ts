import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/sdf-tab.css';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import {download} from '../utils/helpers';
import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {isValidSequence} from './sequence-codes-tools';
import {drawMolecule} from '../utils/structures-works/draw-molecule';
import {RichTextInput, demoColorizer} from '../utils/rich-text-input';

/** Data associated with strands */
type StrandData = {
  strand: string,
  invert: boolean
}

/** Get a molfile for a single strand */
function getMolfileForStrand(strand: string, invert: boolean): string {
  if (strand === '')
    return '';
  const validationOutput = isValidSequence(strand, null);
  const format = validationOutput.synthesizer![0];
  let molfile = '';
  try {
    molfile = sequenceToMolV3000(strand, invert, false, format!);
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error(errStr);
  }
  return molfile;
}

/** Get molfile for single strand or linked strands */
function getLinkedMolfile(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean
): string {
  const nonEmptyStrands = [ss, as, as2].filter((item) => item.strand !== '');
  if (nonEmptyStrands.length === 1) {
    return getMolfileForStrand(nonEmptyStrands[0].strand, nonEmptyStrands[0].invert);
  } else {
    const ssMol = getMolfileForStrand(ss.strand, ss.invert);
    const asMol = getMolfileForStrand(as.strand, as.invert);
    const as2Mol = getMolfileForStrand(as2.strand, as2.invert);

    // select only the non-empty anti-strands
    const antiStrands = [asMol, as2Mol].filter((item) => item !== '');
    const resultingMolfile = linkStrandsV3000({senseStrands: [ssMol], antiStrands: antiStrands}, useChiral);

    return resultingMolfile;
  }
}

/** Save sdf in case ss and as (and optionally as2) strands entered */
export function saveSdf(
  ss: StrandData, as: StrandData, as2: StrandData, useChiral: boolean,
  oneEntity: boolean
): void {
  const nonEmptyStrands = [ss.strand, as.strand, as2.strand].filter((item) => item !== '');
  if (
    nonEmptyStrands.length === 0 ||
    nonEmptyStrands.length === 1 && ss.strand === ''
  ) {
    grok.shell.warning('Enter SS and AS/AS2 to save SDF');
  } else {
    let result: string;
    if (oneEntity) {
      result = getLinkedMolfile(ss, as, as2, useChiral) + '\n$$$$\n';
    } else {
      const ssMol = getMolfileForStrand(ss.strand, ss.invert);
      const asMol = getMolfileForStrand(as.strand, as.invert);
      const as2Mol = getMolfileForStrand(as2.strand, as2.invert);
      result = ssMol + '\n' +
        `> <Sequence>\nSense Strand\n$$$$\n`;
      if (asMol) {
        result += asMol + '\n' +
        `> <Sequence>\nAnti Sense\n$$$$\n`;
      }
      if (as2Mol) {
        result += as2Mol + '\n' +
          `> <Sequence>\nAnti Sense 2\n$$$$\n`;
      }
    }

    // construct date-time in the form yyyy-mm-dd_hh-mm-ss
    const date = new Date();
    function pad(x: number): string {
      return (x >= 10) ? x.toString() : '0' + x.toString();
    }
    const dateString: string = date.getFullYear() + '-' + pad(date.getMonth() + 1) +
      '-' + pad(date.getDate()) + '_' + pad(date.getHours()) + '-' +
      pad(date.getMinutes()) + '-' + pad(date.getSeconds());

    download(`SequenceTranslator-${dateString}.sdf`, encodeURIComponent(result));
  }
}

/** UI of the SDF tab on the application's view */
export function getSdfTab(): HTMLDivElement {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  // default input values
  const straight = '5′ → 3′';
  const inverse = '3′ → 5′';
  let invertSS = false;
  let invertAS = false;
  let invertAS2 = false;

  // text inputs
  const ssInputBase = ui.textInput('', '', () => { onInput.next(); });
  const asInputBase = ui.textInput('', '', () => { onInput.next(); });
  const as2InputBase = ui.textInput('', '', () => { onInput.next(); });

  const invalidSubsequenceColorizer = function(input: string): HTMLSpanElement[] {
    // validate sequence
    const cutoff = isValidSequence(input, null).indexOfFirstNotValidChar;
    const isValid = cutoff < 0 || input === '';
    const greyTextSpan = ui.span([]);
    $(greyTextSpan).css('-webkit-text-fill-color', 'var(--grey-6)');
    const redTextSpan = ui.span([]);
    $(redTextSpan).css('-webkit-text-fill-color', 'red');

    if (!isValid) {
      greyTextSpan.innerHTML = input.slice(0, cutoff);
      redTextSpan.innerHTML = input.slice(cutoff);
    } else { greyTextSpan.innerHTML = input; }
    return [greyTextSpan, redTextSpan];
  };

  const ssRichInput = new RichTextInput(ssInputBase, invalidSubsequenceColorizer);
  const asRichInput = new RichTextInput(asInputBase, invalidSubsequenceColorizer);
  const as2RichInput = new RichTextInput(as2InputBase, demoColorizer);


  // bool inputs
  const saveEntity = ui.boolInput('Save as one entity', true);
  ui.tooltip.bind(saveEntity.root, 'Save SDF with all strands in one molfile');
  const useChiralInput = ui.boolInput('Use chiral', true);
  // todo: compose tooltip message:
  // ui.tooltip.bind(useChiralInput.root, '');

  // choice inputs
  const ssDirection = ui.choiceInput('SS direction', straight, [straight, inverse]);
  ssDirection.onChanged(() => {
    invertSS = ssDirection.value === inverse;
    onInput.next();
  });

  const asDirection = ui.choiceInput('AS direction', straight, [straight, inverse]);
  asDirection.onChanged(() => {
    invertAS = asDirection.value === inverse;
    onInput.next();
  });

  const as2Direction = ui.choiceInput('AS2 direction', straight, [straight, inverse]);
  as2Direction.onChanged(() => {
    invertAS2 = as2Direction.value === inverse;
    onInput.next();
  });

  // labels
  const ssLabel = ui.label('Sense Strand');
  const asLabel = ui.label('Anti Sense');
  const as2Label = ui.label('Anti Sense 2');

  // table layout
  const tableLayout = ui.table(
    ['ss', 'as1', 'as2'], (row, index) => {
      switch (row) {
      case 'ss':
        return [ssLabel, ssRichInput.root, ssDirection.root];
      case 'as1':
        return [asLabel, asRichInput.root, asDirection.root];
      case 'as2':
        return [as2Label, as2RichInput.root, as2Direction.root];
      }
    }, ['', '', '']
  );

  // text input label style
  for (const item of [ssLabel, asLabel, as2Label]) {
    item.parentElement!.classList.add('sdf-input-form', 'sdf-text-input-label');
    $(item.parentElement!).css('padding-top', '3px'); // otherwise overridden
  }

  // choice input label style
  for (const item of [ssDirection.root, asDirection.root, as2Direction.root])
    item.parentElement!.classList.add('sdf-input-form', 'sdf-choice-input-label');

  for (const item of [ssInputBase, asInputBase, as2InputBase]) {
    // text area's parent td
    item.root.parentElement!.classList.add('sdf-text-input-td');
  }

  // molecule image container
  const moleculeImgDiv = ui.block([]);
  $(moleculeImgDiv).addClass('sdf-mol-img');

  DG.debounce<string>(onInput, 300).subscribe(async () => {
    let molfile = '';
    try {
      molfile = getLinkedMolfile(
        {strand: ssInputBase.value, invert: invertSS},
        {strand: asInputBase.value, invert: invertAS},
        {strand: as2InputBase.value, invert: invertAS2}, useChiralInput.value!
      );
    } catch (err) {
      const errStr = errorToConsole(err);
      console.error(errStr);
    }
    // todo: calculate relative numbers
    const canvasWidth = 650;
    const canvasHeight = 150;
    await drawMolecule(moleculeImgDiv, canvasWidth, canvasHeight, molfile);
    // should the canvas be returned from the above function?
    $(moleculeImgDiv).find('canvas').css('float', 'inherit');
  });

  const saveButton = ui.buttonsInput([
    ui.bigButton('Save SDF', () =>
      saveSdf(
        {strand: ssInputBase.value, invert: invertSS},
        {strand: asInputBase.value, invert: invertAS},
        {strand: as2InputBase.value, invert: invertAS2},
        useChiralInput.value!,
        saveEntity.value!
      )
    )
  ]);

  const boolInputsAndButtonArray = [saveEntity.root, useChiralInput.root, saveButton];
  const boolInputsAndButton = ui.divV(boolInputsAndButtonArray);
  for (const item of boolInputsAndButtonArray)
    $(item).addClass('sdf-bool-button-block');

  const bottomDiv = ui.divH([boolInputsAndButton, moleculeImgDiv]);
  $(bottomDiv).addClass('sdf-bottom');

  const sdfTabBody = ui.divV([tableLayout, bottomDiv]);
  $(sdfTabBody).addClass('sdf-body');

  return sdfTabBody;
}
