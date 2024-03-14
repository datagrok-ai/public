/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import './style.css';
import $ from 'cash-dom';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {highlightInvalidSubsequence} from '../../common/ui-components/colored-input/input-painters';
import {getLinkedMolfile, saveSdf} from '../model/oligo-structure';
import {ColoredTextInput} from '../../common/ui-components/colored-input/colored-text-input';
import {MoleculeImage} from '../../common/ui-components/molecule-img';
import {StrandData} from '../model/oligo-structure';

const enum DIRECTION {
  STRAIGHT = '5′ → 3′',
  INVERSE = '3′ → 5′',
};
const STRANDS = ['ss', 'as', 'as2'] as const;

export class StructureAppLayout {
  constructor() {
    this.onInput = new rxjs.Subject<string>();
    this.onInvalidInput = new rxjs.Subject<string>();
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
    this.moleculeImgDiv = ui.block([]);
    $(this.moleculeImgDiv).addClass('st-structure-mol-img');

    DG.debounce<string>(this.onInput, 300).subscribe(async () => {
      await this.updateMoleculeImg();
    });

    DG.debounce<string>(this.onInvalidInput, 1000).subscribe(async () => {
      grok.shell.warning('Insert Sense strand');
    });
  }

  private onInput: rxjs.Subject<string>;
  private onInvalidInput: rxjs.Subject<string>;
  private useChiralInput: DG.InputBase<boolean | null>;
  private saveAllStrandsInput: DG.InputBase<boolean | null>;
  private inputBase: {[key: string]: DG.InputBase<string>};
  private directionInversion: {[key: string]: boolean};
  private moleculeImgDiv: HTMLDivElement;

  async getHtmlDivElement(): Promise<HTMLDivElement> {
    const tableLayout = this.getTableInput();
    const boolInputsAndButton = this.getBoolInputsAndButton();
    await this.updateMoleculeImg();
    const bottomDiv = ui.divH([boolInputsAndButton, this.moleculeImgDiv]);
    $(bottomDiv).addClass('st-structure-bottom');

    const layout = ui.divV([tableLayout, bottomDiv]);
    $(layout).addClass('st-structure-body');

    return layout;
  }

  private getBoolInputsAndButton(): HTMLDivElement {
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
      $(item).addClass('st-structure-bool-button-block');
    return boolInputsAndButton;
  }

  private getTableInput(): HTMLTableElement {
    const coloredInput = Object.fromEntries(
      STRANDS.map(
        (key) => [key, new ColoredTextInput(this.inputBase[key], highlightInvalidSubsequence)]
      )
    );

    const directionChoiceInput = Object.fromEntries(
      STRANDS.map(
        (key, idx) => {
          const selected = (idx === 0) ? DIRECTION.STRAIGHT : DIRECTION.INVERSE;
          return [key, ui.choiceInput(
            `${key.toUpperCase()} direction`, selected, [DIRECTION.STRAIGHT, DIRECTION.INVERSE]
          )];
        }
      )
    );

    STRANDS.forEach((strand, idx) => {
      directionChoiceInput[strand].onChanged(() => {
        let value = directionChoiceInput[strand].value === DIRECTION.INVERSE;
        // warning: the next line is necessary until the legacy notion of direction used in the molfile generation gets fixed
        if (idx > 0) value = !value;
        this.directionInversion[strand] = value;
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

    const clearBlock = Object.fromEntries(
      STRANDS.map(
        (key) => {
          const clearIcon = ui.icons.delete(() => { coloredInput[key].inputBase.value = ''; });
          const clearButton = ui.button(clearIcon, () => {});
          ui.tooltip.bind(clearButton, `Clear ${key.toUpperCase()}`);
          return [key, clearIcon];
        }
      ));

    const tableRows = STRANDS.map((strand) => {
      return {
        label: label[strand],
        textInput: coloredInput[strand].root,
        clear: clearBlock[strand],
        choiceInput: directionChoiceInput[strand].root,
      };
    });
    const tableLayout = ui.table(
      tableRows, (item) => [item.label, item.textInput, item.clear, item.choiceInput]);
    $(tableLayout).css('margin-top', '10px');

    for (const strand of STRANDS) {
      let element = label[strand].parentElement!;
      element.classList.add('st-structure-input-form');
      // the following line is necessary because otherwise overridden by
      // d4-item-table class
      $(element).css('padding-top', '3px');

      element = directionChoiceInput[strand].root.parentElement!;
      element.classList.add('st-structure-input-form', 'st-structure-direction-choice');

      element = this.inputBase[strand].root.parentElement!;
      element.classList.add('st-structure-text-input-td');
    }
    return tableLayout;
  }

  private getStrandData() {
    return Object.fromEntries(
      STRANDS.map((strand, idx) => {
        const invert = this.directionInversion[strand];
        return [strand, {
          strand: this.inputBase[strand].value.replace(/\s*/g, ''),
          invert: invert
        }];
      })
    );
  }

  private getMolfile(ss: StrandData, as: StrandData, as2: StrandData): string {
    // if (ss.strand === '' && (as.strand !== '' || as2.strand !== '')) {
    //   this.onInvalidInput.next();
    //   return '';
    // }

    return getLinkedMolfile(ss, as, as2, this.useChiralInput.value!);
  }

  private async updateMoleculeImg(): Promise<void> {
    let molfile = '';
    try {
      const strandData = this.getStrandData();
      if (Object.values(strandData).some((data) => data.strand !== ''))
        molfile = this.getMolfile(strandData.ss, strandData.as, strandData.as2);
    } catch (err) {
      const errStr = errorToConsole(err);
      console.error(errStr);
    }
    // todo: compute relative numbers
    const canvasWidth = 650;
    const canvasHeight = 150;
    const molImgObj = new MoleculeImage(molfile);
    await molImgObj.drawMolecule(this.moleculeImgDiv, canvasWidth, canvasHeight);
    // should the canvas be returned from the above function?
    $(this.moleculeImgDiv).find('canvas').css('float', 'inherit');
  }
}
