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
  constructor(onInputChanged: (input: string) => void) {
    this._inputSequence = DEFAULT_INPUT;
    this._onInputChanged = onInputChanged;
  }

  private _inputSequence: string;
  private _onInputChanged: (input: string) => void;

  get inputSequence() { return this._inputSequence; }

  async getHtmlElement(): Promise<HTMLDivElement> {
    return await getMainTab(this._onInputChanged);
  }

  set inputSequence(newSequence: string) {
    // todo: validation of the inserted sequence
    this._inputSequence = newSequence;
  }
}

/** Produce HTML div for the 'main' tab */
export async function getMainTab(onSequenceChanged: (seq: string) => void): Promise<HTMLDivElement> {
  async function updateTableAndMolecule(sequence: string): Promise<void> {
    moleculeImgDiv.innerHTML = '';
    outputTableDiv.innerHTML = '';
    const pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');

    try {
      sequence = sequence.replace(/\s/g, '');
      const output = isValidSequence(sequence, null);
      inputFormatChoiceInput.value = output.synthesizer![0];
      const outputSequenceObj = convertSequence(sequence, output);
      const tableRows = [];

      for (const key of Object.keys(outputSequenceObj).slice(1)) {
        tableRows.push({
          format: key,
          sequence: ('indexOfFirstInvalidChar' in outputSequenceObj) ?
            ui.divH([]) : //@ts-ignore
            ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key])
              .then(() => grok.shell.info(SEQUENCE_COPIED_MSG)), SEQ_TOOLTIP_MSG, ''),
        });
      }

      outputTableDiv.append(
        ui.div([
          ui.table(tableRows, (item) => [item.format, item.sequence], ['OUTPUT FORMAT', 'OUTPUT SEQUENCE'])
        ])
        // ui.div([
        //   DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) =>
        //     [item.key, item.value], ['Format', 'Sequence']).root,
        // ])
      );

      if (outputSequenceObj.type !== undefinedInputSequence && outputSequenceObj.Error !== undefinedInputSequence) {
        const formCanvasWidth = 500;
        const formCanvasHeight = 170;
        const molfile = sequenceToMolV3000(
          inputSequenceBase.value.replace(/\s/g, ''), false, true,
          output.synthesizer![0]
        );
        await drawMolecule(moleculeImgDiv, formCanvasWidth, formCanvasHeight, molfile);
      } else {
        moleculeImgDiv.innerHTML = '';
      }
    } finally {
      pi.close();
    }
  }

  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  const inputFormatChoiceInput = ui.choiceInput('', INPUT_FORMATS.GCRS, Object.values(INPUT_FORMATS));
  inputFormatChoiceInput.onInput(async () => {
    await updateTableAndMolecule(inputSequenceBase.value.replace(/\s/g, ''));
  });
  const inputSequenceBase = ui.textInput('', DEFAULT_INPUT,
    (sequence: string) => {
      // Send event to DG.debounce()
      onInput.next(sequence);
    });

  const sequenceColoredInput = new ColoredTextInput(inputSequenceBase, highlightInvalidSubsequence);

  DG.debounce<string>(onInput, 300).subscribe(async (sequence) => {
    await updateTableAndMolecule(sequence);
    onSequenceChanged(sequence);
  });

  const downloadMolfileButton = ui.button(
    'Get Molfile',
    async () => {
      const clearSequence = inputSequenceBase.value.replace(/\s/g, '');
      const result = sequenceToMolV3000(inputSequenceBase.value.replace(/\s/g, ''), false, false,
        inputFormatChoiceInput.value!);
      download(clearSequence + '.mol', encodeURIComponent(result));
    },
    'Save .mol file');

  const copySmilesButton = ui.button(
    'Copy SMILES',
    () => {
      navigator.clipboard.writeText(
        sequenceToSmiles(inputSequenceBase.value.replace(/\s/g, ''), false, inputFormatChoiceInput.value!)
      ).then(() => grok.shell.info(SEQUENCE_COPIED_MSG));
    },
    'Copy SMILES string corresponding to the sequence');

  const formatChoiceInput = ui.div([inputFormatChoiceInput], {style: {padding: '5px 0'}});

  // const sequenceInputLabel = ui.label('Sequence:');

  // const upperBlock = ui.divH([
  //   ui.label('Sequence'), sequenceColoredInput.root,
  //   formatChoiceInput,
  // ]);

  // table layout
  const upperBlock = ui.table(
    ['ss'], (row) => {
      switch (row) {
      case 'ss':
        // return [sequenceInputLabel, sequenceColoredInput.root, formatChoiceInput];
        return [formatChoiceInput, sequenceColoredInput.root];
      }
    }, ['INPUT FORMAT', 'INPUT SEQUENCE']
  );

  const outputTableDiv = ui.div([]);
  const outputTable = ui.block([
    // ui.h1('Output'),
    outputTableDiv,
    downloadMolfileButton,
    copySmilesButton,
  ]);

  const moleculeImgDiv = ui.block([]);

  const mainTabBody = ui.box(
    // ui.splitH([
    //   ui.splitV([
    //     ui.panel([
    //       // appMainDescription,
    //       upperBlock,
    //       formatChoiceInput,
    //       outputTable,
    //       moleculeImgDiv,
    //     ], 'st-main-sequence'),
    //   ]),
    // ], {style: {height: '100%', width: '100%'}})
    ui.div([
      // appMainDescription,
      upperBlock,
      outputTable,
      moleculeImgDiv,
    ], {style: {paddingTop: '20px', paddingLeft: '20px'}})
  );

  await updateTableAndMolecule(DEFAULT_INPUT);
  return mainTabBody;
}
