/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import './style.css';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import {ColoredTextInput} from '../../common/view/components/colored-input/colored-text-input';
import {MoleculeImage} from '../../common/view/components/molecule-img';
import {APP_NAME} from '../../common/view/const';
import {IsolatedAppUIBase} from '../../common/view/isolated-app-ui';
import {getLinkedMolfile, saveSdf, StrandData} from '../model/oligo-structure';
import {ITranslationHelper} from '../../../types';

import {_package} from '../../../package';

const enum DIRECTION {
  STRAIGHT = '5′ → 3′',
  INVERSE = '3′ → 5′',
};
const STRANDS = ['ss', 'as', 'as2'] as const;

class StructureAppLayout {
  private readonly th: ITranslationHelper;

  constructor() {
    this.th = _package;
    this.onInput = new rxjs.Subject<string>();
    this.onInvalidInput = new rxjs.Subject<string>();
    this.inputBase = Object.fromEntries(
      STRANDS.map(
        (key) => {
          const input = ui.input.textArea(key.toUpperCase(), {value: '', onValueChanged: () => {
            this.onInput.next();
            // WARNING: this fine tuning is necessary to fix layout within ui.form
            // js-api version ^1.21
            $(input.root.getElementsByTagName('div')).css('padding-left', '38px');
          }});
          return [key, input];
        }
      )
    );
    this.useChiralInput = ui.input.bool('Use chiral', {value: true});
    this.saveAllStrandsInput = ui.input.bool('Save as one entity', {value: true});
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
  private inputBase: { [key: string]: DG.InputBase<string> };
  private directionInversion: { [key: string]: boolean };
  private moleculeImgDiv: HTMLDivElement;

  async getHtmlDivElement(th: ITranslationHelper): Promise<HTMLDivElement> {
    const tableLayout = this.getTableInput(th);
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
          this.useChiralInput.value!, this.saveAllStrandsInput.value!, this.th);
      })
    ]);

    const boolInputsAndButtonArray = [this.saveAllStrandsInput.root, this.useChiralInput.root, saveButton];
    const boolInputsAndButton = ui.divV(boolInputsAndButtonArray);
    for (const item of boolInputsAndButtonArray)
      $(item).addClass('st-structure-bool-button-block');
    return boolInputsAndButton;
  }

  private getTableInput(th: ITranslationHelper): HTMLElement {
    const coloredInput = Object.fromEntries(
      STRANDS.map(
        (key) => [key, new ColoredTextInput(this.inputBase[key], th.highlightInvalidSubsequence)]
      )
    );

    const directionChoiceInput = Object.fromEntries(
      STRANDS.map(
        (key, idx) => {
          const selected = (idx === 0) ? DIRECTION.STRAIGHT : DIRECTION.INVERSE;
          return [key, ui.input.choice(
            `${key.toUpperCase()} direction`, {value: selected, items: [DIRECTION.STRAIGHT, DIRECTION.INVERSE]}
          )];
        }
      )
    );

    STRANDS.forEach((strand, idx) => {
      directionChoiceInput[strand].onChanged.subscribe((value) => {
        let inputValue = value === DIRECTION.INVERSE;
        // warning: the next line is necessary
        // until the legacy notion of direction used in the molfile generation gets fixed
        if (idx > 0) inputValue = !inputValue;
        this.directionInversion[strand] = inputValue;
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

    const sequenseInputs = ui.form([], 'st-structure-inputs');
    const clearButtons = ui.divV([], 'st-structure-clear-buttons');
    const directionInputs = ui.form([], 'st-direction-inputs');

    for (const strand of STRANDS) {
      sequenseInputs.append(this.inputBase[strand].root);
      clearButtons.append(clearBlock[strand]);
      directionInputs.append(directionChoiceInput[strand].root);
    }
    return ui.divH([sequenseInputs, clearButtons, directionInputs]);
  }

  private getStrandData() {
    return Object.fromEntries(
      STRANDS.map((strand) => {
        const invert = this.directionInversion[strand];
        return [strand, {
          strand: this.inputBase[strand].value.replace(/\s*/g, ''),
          invert: invert
        }];
      })
    );
  }

  private getMolfile(ss: StrandData, as: StrandData, as2: StrandData): string {
    return getLinkedMolfile(ss, as, as2, this.useChiralInput.value!, this.th);
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

export class OligoStructureUI extends IsolatedAppUIBase {
  constructor(
    private readonly th: ITranslationHelper
  ) {
    super(APP_NAME.STRUCTURE);
    this.layout = new StructureAppLayout();
  }

  private readonly layout: StructureAppLayout;

  protected getContent(): Promise<HTMLDivElement> {
    return this.layout.getHtmlDivElement(this.th);
  }
}
