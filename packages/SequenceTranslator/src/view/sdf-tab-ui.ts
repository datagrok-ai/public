/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// external modules dependencies
import * as rxjs from 'rxjs';
import $ from 'cash-dom';

// datagrok libraries dependencies
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

// inner dependencies
import {drawMolecule} from '../utils/structures-works/draw-molecule';
import {highlightInvalidSubsequence, demoPainter} from './input-painters';
import {getLinkedMolfile, saveSdf} from '../sdf-tab/sdf-tab';
import {ColoredTextInput} from '../utils/colored-text-input';

export class SdfTabUI {
  constructor() {
    this._htmlDivElement = getSdfTab();
  }

  private _htmlDivElement;

  get htmlDivElement() { return this._htmlDivElement; }
}

/** UI of the SDF tab on the application's view */
function getSdfTab(): HTMLDivElement {
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

  const ssColoredInput = new ColoredTextInput(ssInputBase, highlightInvalidSubsequence);
  const asColoredInput = new ColoredTextInput(asInputBase, demoPainter);
  const as2ColoredInput = new ColoredTextInput(as2InputBase, demoPainter);

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
        return [ssLabel, ssColoredInput.root, ssDirection.root];
      case 'as1':
        return [asLabel, asColoredInput.root, asDirection.root];
      case 'as2':
        return [as2Label, as2ColoredInput.root, as2Direction.root];
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
        {strand: ssInputBase.value.replace(/\s*/g, ''), invert: invertSS},
        {strand: asInputBase.value.replace(/\s*/g, ''), invert: invertAS},
        {strand: as2InputBase.value.replace(/\s*/g, ''), invert: invertAS2}, useChiralInput.value!
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
