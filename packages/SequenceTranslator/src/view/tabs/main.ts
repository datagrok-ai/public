/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import '../css/main-tab.css';
import $ from 'cash-dom';

import {highlightInvalidSubsequence} from '../utils/colored-input/input-painters';
import {ColoredTextInput} from '../utils/colored-input/colored-text-input';
import {INPUT_FORMATS} from '../../model/const';
import {SYNTHESIZERS as FORMAT} from '../../model/const';
import {SequenceToSmilesConverter} from '../../model/sequence-to-structure-utils/sequence-to-smiles';
import {SequenceToMolfileConverter} from '../../model/sequence-to-structure-utils/sequence-to-molfile';
import {convertSequence} from '../../model/format-translation/conversion-utils';
import {MoleculeImage} from '../utils/molecule-img';
import {download} from '../../model/helpers';
import {SEQUENCE_COPIED_MSG, SEQ_TOOLTIP_MSG, DEFAULT_INPUT} from '../const/main-tab';
import {FormatDetector} from '../../model/parsing-validation/format-detector';
import {SequenceValidator} from '../../model/parsing-validation/sequence-validator';
import {FormatConverter} from '../../model/format-translation/format-converter';

export class MainTabUI {
  constructor() {
    this.moleculeImgDiv = ui.block([]);
    this.outputTableDiv = ui.div([]);
    this.formatChoiceInput = ui.choiceInput('', INPUT_FORMATS.GCRS, Object.values(INPUT_FORMATS));
    this.formatChoiceInput.onInput(async () => {
      await this.updateLayout();
    });
    this.sequenceInputBase = ui.textInput('', DEFAULT_INPUT,
      (sequence: string) => {
        // Send event to DG.debounce()
        this.onInput.next(sequence);
      });

    this.init();

    DG.debounce<string>(this.onInput, 300).subscribe(async () => {
      this.init();
      await this.updateLayout();
    });
  }

  // todo: reduce # of state vars by further refactoring legacy code
  private onInput = new rxjs.Subject<string>();
  private moleculeImgDiv: HTMLDivElement;
  private outputTableDiv: HTMLDivElement;
  private formatChoiceInput: DG.InputBase;
  private sequenceInputBase: DG.InputBase;
  private molfile: string;
  private sequence: string;
  private format: FORMAT | null;

  async getHtmlElement(): Promise<HTMLDivElement> {
    const sequenceColoredInput = new ColoredTextInput(this.sequenceInputBase, highlightInvalidSubsequence);

    const downloadMolfileButton = ui.button(
      'Get Molfile',
      () => { this.saveMolfile(); },
      'Save sequence as Molfile V3000');

    const copySmilesButton = ui.button(
      'Copy SMILES',
      () => { this.copySmiles(); },
      'Copy SMILES for the sequence');

    const formatChoiceInput = ui.div([this.formatChoiceInput]);

    const inputTableRow = {
      format: formatChoiceInput,
      textInput: sequenceColoredInput.root,
    };
    const upperBlock = ui.table(
      [inputTableRow], (item) => [item.format, item.textInput]
    );
    upperBlock.classList.add('st-main-input-table');

    const outputTable = ui.block([
      this.outputTableDiv,
      downloadMolfileButton,
      copySmilesButton,
    ]);

    const mainTabBody = ui.box(
      ui.div([
        upperBlock,
        outputTable,
        this.moleculeImgDiv,
      ], {style: {paddingTop: '20px', paddingLeft: '20px'}})
    );

    await this.updateLayout();
    return mainTabBody;
  }

  private saveMolfile(): void {
    const result = (new SequenceToMolfileConverter(this.sequence, false,
      this.formatChoiceInput.value!)).convert();
    download(this.sequence + '.mol', encodeURIComponent(result));
  }

  private copySmiles(): void {
    const smiles = (new SequenceToSmilesConverter(this.sequence, false, this.formatChoiceInput.value!)).convert();
    navigator.clipboard.writeText(smiles).then(
      () => grok.shell.info(SEQUENCE_COPIED_MSG)
    );
  }

  private updateTable(): void {
    this.outputTableDiv.innerHTML = '';
    // todo: does not detect correctly (U-A)(U-A)
    const indexOfInvalidChar = (!this.format) ? 0 : (new SequenceValidator(this.sequence)).getInvalidCodeIndex(this.format!);
    const outputSequenceObj = convertSequence(this.sequence, indexOfInvalidChar, this.format);
    const tableRows = [];

    for (const key of Object.keys(outputSequenceObj).slice(1)) {
      const sequence = ('indexOfFirstInvalidChar' in outputSequenceObj) ?
        ui.divH([]) :
        ui.link(
          outputSequenceObj[key],
          () => navigator.clipboard.writeText(outputSequenceObj[key])
            .then(() => grok.shell.info(SEQUENCE_COPIED_MSG)),
          SEQ_TOOLTIP_MSG, ''
        );
      tableRows.push({
        format: key,
        sequence: sequence,
      });
    }
    const outputValues = ui.table(tableRows, (item) => [item.format, item.sequence], ['OUTPUT FORMAT', 'OUTPUT SEQUENCE']);
    outputValues.classList.add('st-main-output-table');

    this.outputTableDiv.append(ui.div([outputValues]));
    // this.outputTableDiv.classList.add()
  }

  private async updateMolImg(): Promise<void> {
    const canvasWidth = 500;
    const canvasHeight = 170;
    const molImgObj = new MoleculeImage(this.molfile);
    await molImgObj.drawMolecule(this.moleculeImgDiv, canvasWidth, canvasHeight);
    // should the canvas be returned from the above function?
    $(this.moleculeImgDiv).find('canvas').css('float', 'inherit');
  }

  // todo: sort mehtods
  private init(): void {
    this.sequence = this.getFormattedSequence();
    // this.format = this.getInputFormat();
    this.format = (new FormatDetector(this.sequence)).getFormat();

    // warning: getMolfile relies on this.format, so the order is important
    this.molfile = this.getMolfile();
    // this.molfile = '';
  }

  private getFormattedSequence(): string {
    return this.sequenceInputBase.value.replace(/\s/g, '');
  }

  private getMolfile(): string {
    if (!this.format)
      return '';
    if (this.format === FORMAT.HELM) {
      const gcrs = (new FormatConverter(this.sequence, this.format).convertTo(FORMAT.GCRS));
      return (new SequenceToMolfileConverter(gcrs, false, FORMAT.GCRS).convert());
    }
    const molfile = (new SequenceToMolfileConverter(this.sequence, false, this.format)).convert();
    return molfile;
  }

  private async updateLayout(): Promise<void> {
    this.formatChoiceInput.value = this.format;
    this.updateTable();
    await this.updateMolImg();
  }
}
