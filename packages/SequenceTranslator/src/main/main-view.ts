import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {convertSequence, undefinedInputSequence, isValidSequence} from '../structures-works/sequence-codes-tools';
import {map} from '../structures-works/map';
import {MODIFICATIONS} from '../structures-works/const';
import {sequenceToSmiles, sequenceToMolV3000} from '../structures-works/from-monomers';
import $ from 'cash-dom';
import {download} from '../helpers';
import {extractAtomDataV3000} from '../structures-works/mol-transformations';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

const defaultInput = 'fAmCmGmAmCpsmU'; // todo: capitalize constants
const sequenceWasCopied = 'Copied'; // todo: wrap hardcoded literals into constants
const tooltipSequence = 'Copy sequence';

export async function mainView(onSequenceChanged: (seq: string) => void): Promise<HTMLDivElement> {
  const onInput: rxjs.Subject<string> = new rxjs.Subject<string>();

  async function updateTableAndMolecule(sequence: string, inputFormat: string): Promise<void> {
    moleculeSvgDiv.innerHTML = '';
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
                '. See tables with valid codes on the right',
              ),
            ]) : //@ts-ignore
            ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key])
              .then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, ''),
        });
      }

      outputTableDiv.append(
        ui.div([
          DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) =>
            [item.key, item.value], ['Code', 'Sequence']).root,
        ]),
      );

      if (outputSequenceObj.type != undefinedInputSequence && outputSequenceObj.Error != undefinedInputSequence) {
        const formCanvasWidth = 500;
        const formCanvasHeight = 170;
        const formCanvas = ui.canvas(
          formCanvasWidth * window.devicePixelRatio, formCanvasHeight * window.devicePixelRatio);
        formCanvas.style.width = `${formCanvasWidth}px`;
        formCanvas.style.height = `${formCanvasHeight}px`;

        formCanvas.addEventListener('click', async () => {
          try {
            const mol = sequenceToMolV3000(
              inputSequenceField.value.replace(/\s/g, ''), false, true,
              output.synthesizer![0],
            );
            console.log(mol);

            const addDiv = ui.div([], {style: {overflowX: 'scroll'}});

            // addDiv size required, but now available before dialog show()
            const coordinates = extractAtomDataV3000(mol);
            const cw: number = $(window).width() * 0.80; // addDiv.clientWidth
            const ch: number = $(window).height() * 0.70; // addDiv.clientHeight
            const molWidth: number = Math.max(...coordinates.x) - Math.min(...coordinates.x);
            const molHeight: number = Math.max(...coordinates.y) - Math.min(...coordinates.y);

            const wR: number = cw / molWidth;
            const hR: number = ch / molHeight;
            const r: number = hR; // Math.max(wR, hR);
            const dlgCanvasWidth = r * molWidth;
            const dlgCanvasHeight = r * molHeight;

            const dlgCanvas = ui.canvas(dlgCanvasWidth * window.devicePixelRatio, dlgCanvasHeight * window.devicePixelRatio);
            dlgCanvas.style.width = `${dlgCanvasWidth}px`;
            dlgCanvas.style.height = `${dlgCanvasHeight}px`;

            // // @ts-ignore
            // OCL.StructureView.drawMolecule(dlgCanvas, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
            // await grok.chem.canvasMol(0, 0, dlgCanvas.width, dlgCanvas.height, dlgCanvas, mol, null,
            //   {setNewCoords: false, normalizeDepiction: false, straightenDepiction: false});
            await grok.functions.call('Chem:canvasMol', {
              x: 0, y: 0, w: dlgCanvas.width, h: dlgCanvas.height, canvas: dlgCanvas,
              molString: mol, scaffoldMolString: '',
              options: {normalizeDepiction: false, straightenDepiction: false}
            });

            addDiv.appendChild(dlgCanvas);
            ui.dialog('Molecule: ' + inputSequenceField.value)
              .add(addDiv)
              .showModal(true);
          } catch (err) {
            const errStr = errorToConsole(err);
            console.error(errStr);
          }
        });
        $(formCanvas).on('mouseover', () => $(formCanvas).css('cursor', 'zoom-in'));
        $(formCanvas).on('mouseout', () => $(formCanvas).css('cursor', 'default'));
        const mol = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, true,
          output.synthesizer![0]);
        // // @ts-ignore
        // OCL.StructureView.drawMolecule(formCanvas, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
        await grok.functions.call('Chem:canvasMol', {
          x: 0, y: 0, w: formCanvas.width, h: formCanvas.height, canvas: formCanvas,
          molString: mol, scaffoldMolString: '',
          options: {normalizeDepiction: false, straightenDepiction: false}
        });
        moleculeSvgDiv.append(formCanvas);
      } else
        moleculeSvgDiv.innerHTML = '';
    } finally {
      pi.close();
    }
  }

  const inputFormatChoiceInput = ui.choiceInput('Input format: ', 'Janssen GCRS Codes', Object.keys(map));
  inputFormatChoiceInput.onInput((format: string) => {
    updateTableAndMolecule(inputSequenceField.value.replace(/\s/g, ''), format);
  });
  const moleculeSvgDiv = ui.block([]);
  const outputTableDiv = ui.div([]);
  const inputSequenceField = ui.textInput('', defaultInput, (sequence: string) => {
    // Send event to DG.debounce()
    onInput.next(sequence);
  });

  DG.debounce<string>(onInput, 300).subscribe((sequence) => {
    updateTableAndMolecule(sequence, inputFormatChoiceInput.value!);
    onSequenceChanged(sequence);
  });

  const asoGapmersGrid = DG.Viewer.grid(
    DG.DataFrame.fromObjects([
      {'Name': '2\'MOE-5Me-rU', 'BioSpring': '5', 'Janssen GCRS': 'moeT'},
      {'Name': '2\'MOE-rA', 'BioSpring': '6', 'Janssen GCRS': 'moeA'},
      {'Name': '2\'MOE-5Me-rC', 'BioSpring': '7', 'Janssen GCRS': 'moe5mC'},
      {'Name': '2\'MOE-rG', 'BioSpring': '8', 'Janssen GCRS': 'moeG'},
      {'Name': '5-Methyl-dC', 'BioSpring': '9', 'Janssen GCRS': '5mC'},
      {'Name': 'ps linkage', 'BioSpring': '*', 'Janssen GCRS': 'ps'},
      {'Name': 'dA', 'BioSpring': 'A', 'Janssen GCRS': 'A, dA'},
      {'Name': 'dC', 'BioSpring': 'C', 'Janssen GCRS': 'C, dC'},
      {'Name': 'dG', 'BioSpring': 'G', 'Janssen GCRS': 'G, dG'},
      {'Name': 'dT', 'BioSpring': 'T', 'Janssen GCRS': 'T, dT'},
      {'Name': 'rA', 'BioSpring': '', 'Janssen GCRS': 'rA'},
      {'Name': 'rC', 'BioSpring': '', 'Janssen GCRS': 'rC'},
      {'Name': 'rG', 'BioSpring': '', 'Janssen GCRS': 'rG'},
      {'Name': 'rU', 'BioSpring': '', 'Janssen GCRS': 'rU'},
    ])!, {showRowHeader: false, showCellTooltip: false, allowEdit: false},
  );

  const omeAndFluoroGrid = DG.Viewer.grid(
    DG.DataFrame.fromObjects([
      {'Name': '2\'-fluoro-U', 'BioSpring': '1', 'Axolabs': 'Uf', 'Janssen GCRS': 'fU'},
      {'Name': '2\'-fluoro-A', 'BioSpring': '2', 'Axolabs': 'Af', 'Janssen GCRS': 'fA'},
      {'Name': '2\'-fluoro-C', 'BioSpring': '3', 'Axolabs': 'Cf', 'Janssen GCRS': 'fC'},
      {'Name': '2\'-fluoro-G', 'BioSpring': '4', 'Axolabs': 'Gf', 'Janssen GCRS': 'fG'},
      {'Name': '2\'OMe-rU', 'BioSpring': '5', 'Axolabs': 'u', 'Janssen GCRS': 'mU'},
      {'Name': '2\'OMe-rA', 'BioSpring': '6', 'Axolabs': 'a', 'Janssen GCRS': 'mA'},
      {'Name': '2\'OMe-rC', 'BioSpring': '7', 'Axolabs': 'c', 'Janssen GCRS': 'mC'},
      {'Name': '2\'OMe-rG', 'BioSpring': '8', 'Axolabs': 'g', 'Janssen GCRS': 'mG'},
      {'Name': 'ps linkage', 'BioSpring': '*', 'Axolabs': 's', 'Janssen GCRS': 'ps'},
    ])!, {showRowHeader: false, showCellTooltip: false, allowEdit: false},
  );

  const overhangModificationsGrid = DG.Viewer.grid(
    DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Name', Object.keys(MODIFICATIONS)),
    ])!, {showRowHeader: false, showCellTooltip: false, allowEdit: false},
  );
  updateTableAndMolecule(defaultInput, inputFormatChoiceInput.value!);

  const codesTablesDiv = ui.splitV([
    ui.box(ui.h2('ASO Gapmers'), {style: {maxHeight: '40px'}}),
    asoGapmersGrid.root,
    ui.box(ui.h2('2\'-OMe and 2\'-F siRNA'), {style: {maxHeight: '40px'}}),
    omeAndFluoroGrid.root,
    ui.box(ui.h2('Overhang modifications'), {style: {maxHeight: '40px'}}),
    overhangModificationsGrid.root,
  ], {style: {maxWidth: '350px'}});

  const appMainDescription = ui.info([
    ui.divText('How to convert one sequence:', {style: {'font-weight': 'bolder'}}),
    ui.divText('Paste sequence into the text field below'),
    ui.divText('\n How to convert many sequences:', {style: {'font-weight': 'bolder'}}),
    ui.divText('1. Drag & drop an Excel or CSV file with sequences into Datagrok'),
    ui.divText('2. Right-click on the column header, then see the \'Convert\' menu'),
    ui.divText('This will add the result column to the right of the table'),
  ], 'Convert oligonucleotide sequences between Nucleotides, BioSpring, Axolabs, Mermade 12 and GCRS representations.');

  $(codesTablesDiv).hide();
  const switchInput = ui.switchInput('Codes', false, (v: boolean) => (v) ?
    $(codesTablesDiv).show() :
    $(codesTablesDiv).hide(),
  );

  const downloadMolFileIcon = ui.iconFA('download', async () => {
    const clearSequence = inputSequenceField.value.replace(/\s/g, '');
    const result = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, false,
      inputFormatChoiceInput.value!);
    download(clearSequence + '.mol', encodeURIComponent(result));
  }, 'Save .mol file');

  const copySmilesIcon = ui.iconFA('copy', () => {
    navigator.clipboard.writeText(
      sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''), false, inputFormatChoiceInput.value!),
    ).then(() => grok.shell.info(sequenceWasCopied));
  }, 'Copy SMILES');

  const topPanel = [
    downloadMolFileIcon,
    copySmilesIcon,
    switchInput.root,
  ];

  const v = grok.shell.v;
  const tabControl = grok.shell.sidebar;
  tabControl.onTabChanged.subscribe((_) => {
    v.setRibbonPanels([(tabControl.currentPane.name == 'MAIN') ? topPanel : []]);
  });
  v.setRibbonPanels([topPanel]);

  return ui.box(
    ui.splitH([
      ui.splitV([
        ui.panel([
          appMainDescription,
          ui.div([
            ui.h1('Input sequence'),
            ui.div([], 'input-base'),
            inputSequenceField.root,
          ], 'inputSequence'),
          ui.div([inputFormatChoiceInput], {style: {padding: '5px 0'}}),
          ui.block([
            ui.h1('Output'),
            outputTableDiv,
          ]),
          moleculeSvgDiv,
        ], 'sequence'),
      ]),
      codesTablesDiv,
    ], {style: {height: '100%', width: '100%'}}),
  );
}
