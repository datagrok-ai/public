/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

// import {drawMolecule} from '../../utils/structures-works/draw-molecule';
import {highlightInvalidSubsequence} from '../input-painters';
import {getLinkedMolfile, saveSdf} from '../../sdf-tab/sdf-tab';
import {ColoredTextInput} from '../../utils/colored-text-input';
import {MoleculeImage} from '../molecule-img';
import {StrandData} from '../../sdf-tab/sdf-tab';

const enum DIRECTION {
  STRAIGHT = '5′ → 3′',
  INVERSE = '3′ → 5′',
};
const STRANDS = ['ss', 'as', 'as2'];

export class SdfTabUI {
  constructor() {
    this.onInput = new rxjs.Subject<string>();
    this.inputBase = Object.fromEntries(
      STRANDS.map(
        (key) => [key, ui.textInput('', '', () => { this.onInput.next(); })]
      )
    );
    this.useChiralInput = ui.boolInput('Use chiral', true);
    this.saveAllStrandsInput = ui.boolInput('Save as one entity', true);
    ui.tooltip.bind(this.saveAllStrandsInput.root, 'Save SDF with all strands in one molfile');
    this.directionInversion = Object.fromEntries(
      STRANDS.map((key) => [key, false])
    );
  }

  private onInput: rxjs.Subject<string>;
  private useChiralInput: DG.InputBase<boolean | null>;
  private saveAllStrandsInput: DG.InputBase<boolean | null>;
  private inputBase: {[key: string]: DG.InputBase<string>};
  private directionInversion: {[key: string]: boolean};

  get htmlDivElement(): HTMLDivElement {
    const tableLayout = this.getTableInput();

    // molecule image container
    const moleculeImgDiv = ui.block([]);
    $(moleculeImgDiv).addClass('st-sdf-mol-img');

    DG.debounce<string>(this.onInput, 300).subscribe(async () => {
      let molfile = '';
      try {
        const strandData = this.getStrandData();
        molfile = this.getMolfile(strandData.ss, strandData.as, strandData.as2);
      } catch (err) {
        const errStr = errorToConsole(err);
        console.error(errStr);
      }
      // todo: calculate relative numbers
      const canvasWidth = 650;
      const canvasHeight = 150;
      const molImgObj = new MoleculeImage(molfile);
      await molImgObj.drawMolecule(moleculeImgDiv, canvasWidth, canvasHeight);
      // await drawMolecule(moleculeImgDiv, canvasWidth, canvasHeight, molfile);
      // should the canvas be returned from the above function?
      $(moleculeImgDiv).find('canvas').css('float', 'inherit');
    });

    const saveButton = ui.buttonsInput([
      ui.bigButton('Save SDF', () => {
        const strandData = this.getStrandData();
        saveSdf(strandData.ss, strandData.as, strandData.as2,
          this.useChiralInput.value!, this.saveAllStrandsInput.value!);
      })
    ]);

    const boolInputsAndButtonArray = [this.saveAllStrandsInput.root, this.useChiralInput.root, saveButton];
    const boolInputsAndButton = ui.divV(boolInputsAndButtonArray);
    for (const item of boolInputsAndButtonArray)
      $(item).addClass('st-sdf-bool-button-block');

    const bottomDiv = ui.divH([boolInputsAndButton, moleculeImgDiv]);
    $(bottomDiv).addClass('st-sdf-bottom');

    const sdfTabBody = ui.divV([tableLayout, bottomDiv]);
    $(sdfTabBody).addClass('st-sdf-body');

    return sdfTabBody;
  }

  private getTableInput(): HTMLTableElement {
    const coloredInput = Object.fromEntries(
      STRANDS.map(
        (key) => [key, new ColoredTextInput(this.inputBase[key], highlightInvalidSubsequence)]
      )
    );

    // todo: compose tooltip message:

    const directionChoiceInput = Object.fromEntries(
      STRANDS.map(
        (key) => [key, ui.choiceInput(
          `${key.toUpperCase()} direction`, DIRECTION.STRAIGHT, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]
        )]
      )
    );

    STRANDS.forEach((strand) => {
      directionChoiceInput[strand].onChanged(() => {
        this.directionInversion[strand] = directionChoiceInput[strand].value === DIRECTION.INVERSE;
        this.onInput.next();
      });
    });

    const labelNames = ['Sense Strand', 'Anti Sense', 'Anti Sense 2'];
    const labelNameMap = new Map(STRANDS.map(
      (key, index) => [key, labelNames[index]]
    ));
    const label = Object.fromEntries(
      STRANDS.map(
        (key) => [key, ui.label(labelNameMap.get(key)!)]
      )
    );

    const tableRows = STRANDS.map((strand) => {
      return {
        label: label[strand],
        textInput: coloredInput[strand].root,
        choiceInput: directionChoiceInput[strand].root,
      };
    });
    const tableLayout = ui.table(
      tableRows, (item) => [item.label, item.textInput, item.choiceInput]);
    $(tableLayout).css('margin-top', '10px');

    for (const strand of STRANDS) {
      let element = label[strand].parentElement!;
      element.classList.add('st-sdf-input-form');
      // the following line is necessary because otherwise overridden by
      // d4-item-table class
      $(element).css('padding-top', '3px');

      element = directionChoiceInput[strand].root.parentElement!;
      element.classList.add('st-sdf-input-form', 'st-sdf-direction-choice');

      element = this.inputBase[strand].root.parentElement!;
      element.classList.add('st-sdf-text-input-td');
    }
    return tableLayout;
  }

  private getStrandData() {
    return Object.fromEntries(
      STRANDS.map((strand) => [strand, {
        strand: this.inputBase[strand].value.replace(/\s*/g, ''),
        invert: this.directionInversion[strand]
      }])
    );
  }

  private getMolfile(ss: StrandData, as: StrandData, as2: StrandData): string {
    return getLinkedMolfile(ss, as, as2, this.useChiralInput.value!);
  }
}
