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

    const enum DIRECTION {
      STRAIGHT = '5′ → 3′',
      INVERSE = '3′ → 5′',
    };

    const strands = ['ss', 'as', 'as2'];
    const directionInversion = Object.fromEntries(
      strands.map((key) => [key, false])
    );

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

    const saveEntity = ui.boolInput('Save as one entity', true);
    ui.tooltip.bind(saveEntity.root, 'Save SDF with all strands in one molfile');
    const useChiralInput = ui.boolInput('Use chiral', true);
    // todo: compose tooltip message:

    const directionChoiceInput = Object.fromEntries(
      strands.map(
        (key) => [key, ui.choiceInput(
          `${key.toUpperCase()} direction`, DIRECTION.STRAIGHT, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]
        )]
      )
    );

    strands.forEach((strand) => {
      directionChoiceInput[strand].onChanged(() => {
        directionInversion[strand] = directionChoiceInput[strand].value === DIRECTION.INVERSE;
        onInput.next();
      });
    });

    const labelNames = ['Sense Strand', 'Anti Sense', 'Anti Sense 2'];
    const labelNameMap = new Map(strands.map(
      (key, index) => [key, labelNames[index]]
    ));
    const label = Object.fromEntries(
      strands.map(
        (key) => [key, ui.label(labelNameMap.get(key)!)]
      )
    );

    const tableRows = strands.map((strand) => {
      return {
        label: label[strand],
        textInput: coloredInput[strand].root,
        choiceInput: directionChoiceInput[strand].root,
      };
    });
    const tableLayout = ui.table(
      tableRows, (item) => [item.label, item.textInput, item.choiceInput]);
    $(tableLayout).css('margin-top', '10px');

    for (const strand of strands) {
      let element = label[strand].parentElement!;
      element.classList.add('sdf-input-form', 'sdf-text-input-label');
      $(element).css('padding-top', '3px'); // otherwise overridden

      element = directionChoiceInput[strand].root.parentElement!;
      element.classList.add('sdf-input-form', 'sdf-choice-input-label');

      element = inputBase[strand].root.parentElement!;
      element.classList.add('sdf-text-input-td');
    }

    // molecule image container
    const moleculeImgDiv = ui.block([]);
    $(moleculeImgDiv).addClass('sdf-mol-img');

    DG.debounce<string>(onInput, 300).subscribe(async () => {
      let molfile = '';
      try {
        const strandData = Object.fromEntries(
          strands.map((strand) => [strand, {
            strand: inputBase[strand].value.replace(/\s*/g, ''),
            invert: directionInversion[strand]
          }])
        );
        molfile = getLinkedMolfile(strandData.ss, strandData.as, strandData.as2, useChiralInput.value!);
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
