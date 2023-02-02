import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';

// todo: elminate completely
import {map} from '../hardcode-to-be-eliminated/map';
import {MODIFICATIONS} from '../hardcode-to-be-eliminated/const';

// todo: unify with lib bio monomers works
import {sequenceToSmiles, sequenceToMolV3000} from '../utils/structures-works/from-monomers';

import {convertSequence, undefinedInputSequence, isValidSequence} from '../sdf-tab/sequence-codes-tools';
import {viewMonomerLib} from '../utils/monomer-lib-viewer';
import {drawMolecule} from '../utils/structures-works/draw-molecule';
import {download} from '../utils/helpers';

const DEFAULT_INPUT = 'fAmCmGmAmCpsmU';
const SEQUENCE_COPIED_MSG = 'Copied';
const SEQ_TOOLTIP_MSG = 'Copy sequence';

/** Produce HTML div for the 'main' tab */
export async function getMainTab(onSequenceChanged: (seq: string) => void): Promise<HTMLDivElement> {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

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
        const indexOfFirstNotValidChar = ('indexOfFirstNotValidChar' in outputSequenceObj) ?
          JSON.parse(outputSequenceObj.indexOfFirstNotValidChar!).indexOfFirstNotValidChar :
          -1;

        tableRows.push({
          'key': key,
          'value': ('indexOfFirstNotValidChar' in outputSequenceObj) ?
            ui.divH([
              ui.divText(sequence.slice(0, indexOfFirstNotValidChar), {style: {color: 'grey'}}),
              ui.tooltip.bind(
                ui.divText(sequence.slice(indexOfFirstNotValidChar), {style: {color: 'red'}}),
                'Expected format: ' + JSON.parse(outputSequenceObj.indexOfFirstNotValidChar!).synthesizer +
                '. See tables with valid codes on the right'
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

  const inputFormatChoiceInput = ui.choiceInput('Input format: ', 'Janssen GCRS Codes', Object.keys(map));
  inputFormatChoiceInput.onInput(() => {
    updateTableAndMolecule(inputSequenceField.value.replace(/\s/g, ''));
  });
  const moleculeImgDiv = ui.block([]);
  const outputTableDiv = ui.div([]);
  const inputSequenceField = ui.textInput('', DEFAULT_INPUT, (sequence: string) => {
    // Send event to DG.debounce()
    onInput.next(sequence);
  });

  DG.debounce<string>(onInput, 300).subscribe((sequence) => {
    updateTableAndMolecule(sequence);
    onSequenceChanged(sequence);
  });

  const downloadMolFileIcon = ui.iconFA('download', async () => {
    const clearSequence = inputSequenceField.value.replace(/\s/g, '');
    const result = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, false,
      inputFormatChoiceInput.value!);
    download(clearSequence + '.mol', encodeURIComponent(result));
  }, 'Save .mol file');

  const copySmilesIcon = ui.iconFA('copy', () => {
    navigator.clipboard.writeText(
      sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''), false, inputFormatChoiceInput.value!)
    ).then(() => grok.shell.info(SEQUENCE_COPIED_MSG));
  }, 'Copy SMILES');

  const mainTabBody = ui.box(
    ui.splitH([
      ui.splitV([
        ui.panel([
          // appMainDescription,
          ui.div([
            ui.h1('Input sequence'),
            ui.div([], 'input-base'),
            inputSequenceField.root,
          ], 'main-input-sequence'),
          ui.div([inputFormatChoiceInput], {style: {padding: '5px 0'}}),
          ui.block([
            ui.h1('Output'),
            outputTableDiv,
          ]),
          moleculeImgDiv,
        ], 'main-sequence'),
      ]),
      // codesTablesDiv,
    ], {style: {height: '100%', width: '100%'}})
  );

  return mainTabBody;
}
