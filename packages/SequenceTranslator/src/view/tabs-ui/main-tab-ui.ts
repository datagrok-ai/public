/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* external dependencies */
import * as rxjs from 'rxjs';

/* datagrok dependencies */

/* internal dependencies */
import {INPUT_FORMATS} from '../../model/const';
// import {map} from '../../hardcode-to-be-eliminated/map';// todo: elminate completely
// import {MODIFICATIONS} from '../../hardcode-to-be-eliminated/const';
// todo: unify with lib bio monomers works
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
        const indexOfFirstInvalidChar = ('indexOfFirstInvalidChar' in outputSequenceObj) ?
          JSON.parse(outputSequenceObj.indexOfFirstInvalidChar!).indexOfFirstInvalidChar :
          -1;

        tableRows.push({
          'key': key,
          'value': ('indexOfFirstInvalidChar' in outputSequenceObj) ?
            ui.divH([
              ui.divText(sequence.slice(0, indexOfFirstInvalidChar), {style: {color: 'grey'}}),
              ui.tooltip.bind(
                ui.divText(sequence.slice(indexOfFirstInvalidChar), {style: {color: 'red'}}),
                'Expected format: ' + JSON.parse(outputSequenceObj.indexOfFirstInvalidChar!).synthesizer +
                '. See tables with valid codes'
              ),
            ]) : //@ts-ignore
            ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key])
              .then(() => grok.shell.info(SEQUENCE_COPIED_MSG)), SEQ_TOOLTIP_MSG, ''),
        });
      }

      outputTableDiv.append(
        ui.div([
          DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) =>
            [item.key, item.value], ['Code', 'Sequence']).root,
        ])
      );

      if (outputSequenceObj.type !== undefinedInputSequence && outputSequenceObj.Error !== undefinedInputSequence) {
        const formCanvasWidth = 500;
        const formCanvasHeight = 170;
        const molfile = sequenceToMolV3000(
          inputSequenceField.value.replace(/\s/g, ''), false, true,
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

  const inputFormatChoiceInput = ui.choiceInput('Input format: ', 'Janssen GCRS Codes', Object.values(INPUT_FORMATS));
  inputFormatChoiceInput.onInput(async () => {
    await updateTableAndMolecule(inputSequenceField.value.replace(/\s/g, ''));
  });
  const inputSequenceField = ui.textInput('', DEFAULT_INPUT, (sequence: string) => {
    // Send event to DG.debounce()
    onInput.next(sequence);
  });

  DG.debounce<string>(onInput, 300).subscribe(async (sequence) => {
    await updateTableAndMolecule(sequence);
    onSequenceChanged(sequence);
  });

  const downloadMolfileButton = ui.bigButton(
    'Download Molfile',
    async () => {
      const clearSequence = inputSequenceField.value.replace(/\s/g, '');
      const result = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, false,
        inputFormatChoiceInput.value!);
      download(clearSequence + '.mol', encodeURIComponent(result));
    },
    'Save .mol file');

  const copySmilesButton = ui.bigButton(
    'Copy SMILES',
    () => {
      navigator.clipboard.writeText(
        sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''), false, inputFormatChoiceInput.value!)
      ).then(() => grok.shell.info(SEQUENCE_COPIED_MSG));
    },
    'Copy SMILES');

  const upperBlock = ui.div([
    downloadMolfileButton,
    copySmilesButton,
    ui.h1('Input sequence'),
    ui.div([], 'st-input-base'),
    inputSequenceField.root,
  ], 'st-main-input-sequence');

  const formatChoiceInput = ui.div([inputFormatChoiceInput], {style: {padding: '5px 0'}});

  const outputTableDiv = ui.div([]);
  const outputTable = ui.block([
    ui.h1('Output'),
    outputTableDiv,
  ]);

  const moleculeImgDiv = ui.block([]);

  const mainTabBody = ui.box(
    ui.splitH([
      ui.splitV([
        ui.panel([
          // appMainDescription,
          upperBlock,
          formatChoiceInput,
          outputTable,
          moleculeImgDiv,
        ], 'st-main-sequence'),
      ]),
    ], {style: {height: '100%', width: '100%'}})
  );

  await updateTableAndMolecule(DEFAULT_INPUT);
  return mainTabBody;
}
