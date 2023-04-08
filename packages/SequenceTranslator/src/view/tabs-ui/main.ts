/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* external dependencies */
import * as rxjs from 'rxjs';
import $ from 'cash-dom';

/* datagrok dependencies */

/* internal dependencies */
import {highlightInvalidSubsequence} from '../input-painters';
import {ColoredTextInput} from '../../utils/colored-text-input';
import {INPUT_FORMATS} from '../../model/const';
import {sequenceToSmiles, sequenceToMolV3000} from '../../utils/structures-works/from-monomers';
import {convertSequence, undefinedInputSequence, isValidSequence} from '../../sdf-tab/sequence-codes-tools';
import {drawMolecule} from '../../utils/structures-works/draw-molecule';
import {download} from '../../utils/helpers';
import {SEQUENCE_COPIED_MSG, SEQ_TOOLTIP_MSG, DEFAULT_INPUT} from '../../view/const/main-tab-const';

export class MainTabUI {
  constructor() {
    this.moleculeImgDiv = ui.block([]);
    this.outputTableDiv = ui.div([]);
    this.inputFormatChoiceInput = ui.choiceInput('', INPUT_FORMATS.GCRS, Object.values(INPUT_FORMATS));
    this.inputFormatChoiceInput.onInput(async () => {
      await this.updateTableAndMolecule(this.inputSequenceBase.value.replace(/\s/g, ''));
    });
    this.inputSequenceBase = ui.textInput('', DEFAULT_INPUT,
      (sequence: string) => {
        // Send event to DG.debounce()
        this.onInput.next(sequence);
      });
  }

  private moleculeImgDiv: HTMLDivElement;
  private outputTableDiv: HTMLDivElement;
  private inputFormatChoiceInput: DG.InputBase;
  private inputSequenceBase: DG.InputBase;
  private onInput = new rxjs.Subject<string>();

  async getHtmlElement(): Promise<HTMLDivElement> {
    return await this.getMainTab((seq: string) => {});
  }

  private async updateTableAndMolecule(sequence: string): Promise<void> {
    this.moleculeImgDiv.innerHTML = '';
    this.outputTableDiv.innerHTML = '';
    const pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');

    try {
      sequence = sequence.replace(/\s/g, '');
      const output = isValidSequence(sequence, null);
      this.inputFormatChoiceInput.value = output.synthesizer![0];
      const outputSequenceObj = convertSequence(sequence, output);
      const tableRows = [];

      for (const key of Object.keys(outputSequenceObj).slice(1)) {
        const sequence = ('indexOfFirstInvalidChar' in outputSequenceObj) ?
          ui.divH([]) :
          ui.link(
            //@ts-ignore // why ts-ignore?
            outputSequenceObj[key],
            //@ts-ignore // why ts-ignore?
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

      if (outputSequenceObj.type !== undefinedInputSequence && outputSequenceObj.Error !== undefinedInputSequence) {
        const formCanvasWidth = 500;
        const formCanvasHeight = 170;
        const molfile = sequenceToMolV3000(
          this.inputSequenceBase.value.replace(/\s/g, ''), false, true,
          output.synthesizer![0]
        );
        await drawMolecule(this.moleculeImgDiv, formCanvasWidth, formCanvasHeight, molfile);
      } else {
        this.moleculeImgDiv.innerHTML = '';
      }
    } finally {
      pi.close();
    }
  }


  private async getMainTab(onSequenceChanged: (seq: string) => void): Promise<HTMLDivElement> {
    const sequenceColoredInput = new ColoredTextInput(this.inputSequenceBase, highlightInvalidSubsequence);

    DG.debounce<string>(this.onInput, 300).subscribe(async (sequence) => {
      await this.updateTableAndMolecule(sequence);
      onSequenceChanged(sequence);
    });

    const downloadMolfileButton = ui.button(
      'Get Molfile',
      async () => {
        const clearSequence = this.inputSequenceBase.value.replace(/\s/g, '');
        const result = sequenceToMolV3000(this.inputSequenceBase.value.replace(/\s/g, ''), false, false,
          this.inputFormatChoiceInput.value!);
        download(clearSequence + '.mol', encodeURIComponent(result));
      },
      'Save .mol file');
    // $(downloadMolfileButton).addClass(OUTLINE_BUTTON_CLASS);

    const copySmilesButton = ui.button(
      'Copy SMILES',
      () => {
        navigator.clipboard.writeText(
          sequenceToSmiles(this.inputSequenceBase.value.replace(/\s/g, ''), false, this.inputFormatChoiceInput.value!)
        ).then(() => grok.shell.info(SEQUENCE_COPIED_MSG));
      },
      'Copy SMILES string corresponding to the sequence');
    // $(copySmilesButton).addClass(OUTLINE_BUTTON_CLASS);

    const formatChoiceInput = ui.div([this.inputFormatChoiceInput], {style: {padding: '5px 0'}});

    // const sequenceInputLabel = ui.label('Sequence:');

    // const upperBlock = ui.divH([
    //   ui.label('Sequence'), sequenceColoredInput.root,
    //   formatChoiceInput,
    // ]);

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

    await this.updateTableAndMolecule(DEFAULT_INPUT);
    return mainTabBody;
  }
}
