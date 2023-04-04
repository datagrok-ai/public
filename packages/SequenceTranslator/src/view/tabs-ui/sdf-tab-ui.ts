/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {drawMolecule} from '../../utils/structures-works/draw-molecule';
import {highlightInvalidSubsequence} from '../input-painters';
import {getLinkedMolfile, saveSdf} from '../../sdf-tab/sdf-tab';
import {ColoredTextInput} from '../../utils/colored-text-input';
import {TabUI} from './tab-ui';

export class SdfTabUI implements TabUI {
  private constructor() { }

  public static init() {
    if (!SdfTabUI.instance)
      SdfTabUI.instance = new SdfTabUI();
    return SdfTabUI.instance;
  }

  private static instance: SdfTabUI;

  get htmlDivElement(): HTMLDivElement {
    const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

    // default input values
    const enum DIRECTION {
      STRAIGHT = '5′ → 3′',
      INVERSE = '3′ → 5′',
    };

    const strands = ['ss', 'as', 'as2'];
    const directionInversion = Object.fromEntries(
      strands.map((key) => [key, false])
    );
    console.log('directionInversion:', directionInversion);

    const inputBase = Object.fromEntries(
      strands.map(
        (key) => [key, ui.textInput('', '', () => { onInput.next(); })]
      )
    );

    const coloredInput = Object.fromEntries(
      strands.map(
        (key) => [key, new ColoredTextInput(inputBase[key], highlightInvalidSubsequence)]
      )
    );

    // bool inputs
    const saveEntity = ui.boolInput('Save as one entity', true);
    ui.tooltip.bind(saveEntity.root, 'Save SDF with all strands in one molfile');
    const useChiralInput = ui.boolInput('Use chiral', true);
    // todo: compose tooltip message:

    // choice inputs
    const ssDirection = ui.choiceInput('SS direction', DIRECTION.STRAIGHT, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]);
    ssDirection.onChanged(() => {
      directionInversion.ss = ssDirection.value === DIRECTION.INVERSE;
      onInput.next();
    });

    const asDirection = ui.choiceInput('AS direction', DIRECTION.STRAIGHT, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]);
    asDirection.onChanged(() => {
      directionInversion.as = asDirection.value === DIRECTION.INVERSE;
      onInput.next();
    });

    const as2Direction = ui.choiceInput('AS2 direction', DIRECTION.STRAIGHT, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]);
    as2Direction.onChanged(() => {
      directionInversion.as2 = as2Direction.value === DIRECTION.INVERSE;
      onInput.next();
    });

    // labels
    const ssLabel = ui.label('Sense Strand');
    const asLabel = ui.label('Anti Sense');
    const as2Label = ui.label('Anti Sense 2');

    const tableLayout = ui.table(
      ['ss', 'as1', 'as2'], (row) => {
        switch (row) {
        case 'ss':
          return [ssLabel, coloredInput.ss.root, ssDirection.root];
        case 'as1':
          return [asLabel, coloredInput.as.root, asDirection.root];
        case 'as2':
          return [as2Label, coloredInput.as2.root, as2Direction.root];
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

    for (const item of [inputBase.ss, inputBase.as, inputBase.as2]) {
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
          {strand: inputBase.ss.value.replace(/\s*/g, ''), invert: directionInversion.ss},
          {strand: inputBase.as.value.replace(/\s*/g, ''), invert: directionInversion.as},
          {strand: inputBase.as2.value.replace(/\s*/g, ''), invert: directionInversion.as2}, useChiralInput.value!
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
          {strand: inputBase.ss.value, invert: directionInversion.ss},
          {strand: inputBase.as.value, invert: directionInversion.as},
          {strand: inputBase.as2.value, invert: directionInversion.as2},
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
}
