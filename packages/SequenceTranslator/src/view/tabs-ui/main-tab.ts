/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* external dependencies */
import * as rxjs from 'rxjs';
import $ from 'cash-dom';

/* datagrok dependencies */

/* internal dependencies */
import {highlightInvalidSubsequence} from '../utils/colored-input/input-painters';
import {ColoredTextInput} from '../utils/colored-input/colored-text-input';
import {INPUT_FORMATS} from '../../model/const';
import {sequenceToSmiles, sequenceToMolV3000} from '../../model/sequence-to-molfile-utils/from-monomers';
import {convertSequence, undefinedInputSequence, isValidSequence} from '../../model/code-converter/conversion-validation-tools';
import {drawMolecule} from '../utils/draw-molecule';
import {download} from '../../utils/helpers';
import {SEQUENCE_COPIED_MSG, SEQ_TOOLTIP_MSG, DEFAULT_INPUT} from '../../view/const/main-tab-const';

export class MainTabUI {
  constructor() {
    this.moleculeImgDiv = ui.block([]);
    this.outputTableDiv = ui.div([]);
    this.inputFormatChoiceInput = ui.choiceInput('', INPUT_FORMATS.GCRS, Object.values(INPUT_FORMATS));
    this.inputFormatChoiceInput.onInput(async () => {
      await this.updateLayout();
    });
    this.inputSequenceBase = ui.textInput('', DEFAULT_INPUT,
      (sequence: string) => {
        // Send event to DG.debounce()
        this.onInput.next(sequence);
      });

    this.init();

    DG.debounce<string>(this.onInput, 300).subscribe(async () => {
      this.init();
      await this.updateLayout();
      // onSequenceChanged(sequence);
    });
  }

  // todo: reduce # of state vars by further refactoring legacy code
  private onInput = new rxjs.Subject<string>();
  private moleculeImgDiv: HTMLDivElement;
  private outputTableDiv: HTMLDivElement;
  private inputFormatChoiceInput: DG.InputBase;
  private inputSequenceBase: DG.InputBase;
  private molfile: string;
  private sequence: string;
  private format: string;

  async getHtmlElement(): Promise<HTMLDivElement> {
    const sequenceColoredInput = new ColoredTextInput(this.inputSequenceBase, highlightInvalidSubsequence);

    const downloadMolfileButton = ui.button(
      'Get Molfile',
      () => { this.saveMolfile(); },
      'Save sequence as Molfile V3000');

    const copySmilesButton = ui.button(
      'Copy SMILES',
      () => { this.copySmiles(); },
      'Copy SMILES for the sequence');

    const formatChoiceInput = ui.div([this.inputFormatChoiceInput]);

    const tableRow = {
      format: formatChoiceInput,
      textInput: sequenceColoredInput.root,
    };
    const upperBlock = ui.table(
      [tableRow], (item) => [item.format, item.textInput]
    );

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
    const result = sequenceToMolV3000(this.sequence, false,
      this.inputFormatChoiceInput.value!);
    download(this.sequence + '.mol', encodeURIComponent(result));
  }

  private copySmiles(): void {
    navigator.clipboard.writeText(
      sequenceToSmiles(this.sequence, false, this.inputFormatChoiceInput.value!)
    ).then(
      () => grok.shell.info(SEQUENCE_COPIED_MSG)
    );
  }

  private updateTable(): void {
    this.outputTableDiv.innerHTML = '';
    const output = isValidSequence(this.sequence, null);
    const outputSequenceObj = convertSequence(this.sequence, output);
    const tableRows = [];

    for (const key of Object.keys(outputSequenceObj).slice(1)) {
      const sequence = ('indexOfFirstInvalidChar' in outputSequenceObj) ?
        ui.divH([]) :
        ui.link(
          //@ts-ignore // why ts-ignore? refactor to remove
          outputSequenceObj[key],
          //@ts-ignore // why ts-ignore? refactor to remove
          () => navigator.clipboard.writeText(outputSequenceObj[key])
            .then(() => grok.shell.info(SEQUENCE_COPIED_MSG)),
          SEQ_TOOLTIP_MSG, ''
        );
      tableRows.push({
        format: key,
        sequence: sequence,
      });
    }

    this.outputTableDiv.append(
      ui.div([
        ui.table(tableRows, (item) => [item.format, item.sequence], ['OUTPUT FORMAT', 'OUTPUT SEQUENCE'])
      ])
    );
  }

  private async updateMolImg(): Promise<void> {
    this.moleculeImgDiv.innerHTML = ''; // ??
    // todo: eliminate after refactoring legacy
    const output = isValidSequence(this.sequence, null);
    const outputSequenceObj = convertSequence(this.sequence, output);
    if (outputSequenceObj.type !== undefinedInputSequence && outputSequenceObj.Error !== undefinedInputSequence) {
      const formCanvasWidth = 500;
      const formCanvasHeight = 170;
      await drawMolecule(this.moleculeImgDiv, formCanvasWidth, formCanvasHeight, this.molfile);
    }
  }

  // todo: sort mehtods
  private init(): void {
    this.sequence = this.getFormattedSequence();
    this.format = this.getInputFormat();

    // warning: getMolfile relies on this.format, so the order is important
    this.molfile = this.getMolfile();
  }

  private getFormattedSequence(): string {
    return this.inputSequenceBase.value.replace(/\s/g, '');
  }

  private getMolfile(): string {
    return sequenceToMolV3000(this.sequence, false, this.format);
  }

  // todo: put synthesizers into an object/enum
  // todo: take into account possible exceptions
  private getInputFormat(): string {
    const vadimsOutput = isValidSequence(this.sequence, null);
    return vadimsOutput.synthesizer![0];
  }

  private async updateLayout(): Promise<void> {
    this.updateTable();
    await this.updateMolImg();
  }
}
