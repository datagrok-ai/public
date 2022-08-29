import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {convertSequence, undefinedInputSequence, isValidSequence} from '../structures-works/sequence-codes-tools';
import {map, MODIFICATIONS} from '../structures-works/map';
import {sequenceToSmiles, sequenceToMolV3000} from '../structures-works/from-monomers';
import $ from 'cash-dom';

const defaultInput = 'fAmCmGmAmCpsmU';
const sequenceWasCopied = 'Copied';
const tooltipSequence = 'Copy sequence';

export function mainView() {
  function updateTableAndMolecule(sequence: string, inputFormat: string): void {
    moleculeSvgDiv.innerHTML = '';
    outputTableDiv.innerHTML = '';
    const pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');
    let errorsExist = false;
    try {
      sequence = sequence.replace(/\s/g, '');
      const output = isValidSequence(sequence, null);
      output.synthesizer = [inputFormat];
      // inputFormatChoiceInput.value = output.synthesizer![0];
      const outputSequenceObj = convertSequence(sequence, output);
      const tableRows = [];

      for (const key of Object.keys(outputSequenceObj).slice(1)) {
        const indexOfFirstNotValidChar = ('indexOfFirstNotValidChar' in outputSequenceObj) ?
          JSON.parse(outputSequenceObj.indexOfFirstNotValidChar!).indexOfFirstNotValidChar :
          -1;
        if ('indexOfFirstNotValidChar' in outputSequenceObj) {
          const indexOfFirstNotValidChar = ('indexOfFirstNotValidChar' in outputSequenceObj) ?
            JSON.parse(outputSequenceObj.indexOfFirstNotValidChar!).indexOfFirstNotValidChar :
            -1;
          if (indexOfFirstNotValidChar != -1)
            errorsExist = true;
        }

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

      if (errorsExist) {
        const synthesizer = JSON.parse(outputSequenceObj.indexOfFirstNotValidChar!).synthesizer.slice(0, -6);
        asoGapmersGrid.onCellPrepare(function(gc) {
          gc.style.backColor = (gc.gridColumn.name == synthesizer) ? 0xFFF00000 : 0xFFFFFFFF;
        });
        omeAndFluoroGrid.onCellPrepare(function(gc) {
          gc.style.backColor = (gc.gridColumn.name == synthesizer) ? 0xFFF00000 : 0xFFFFFFFF;
        });
        switchInput.enabled = true;
      } else {
        asoGapmersGrid.onCellPrepare(function(gc) {gc.style.backColor = 0xFFFFFFFF;});
        omeAndFluoroGrid.onCellPrepare(function(gc) {gc.style.backColor = 0xFFFFFFFF;});
      }

      outputTableDiv.append(
        ui.div([
          DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) =>
            [item.key, item.value], ['Code', 'Sequence']).root,
        ]),
      );

      if (outputSequenceObj.type != undefinedInputSequence && outputSequenceObj.Error != undefinedInputSequence) {
        const canvas = ui.canvas(300, 170);
        canvas.addEventListener('click', () => {
          const canv = ui.canvas($(window).width(), $(window).height());
          const mol = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, true,
            output.synthesizer![0]);
          // @ts-ignore
          OCL.StructureView.drawMolecule(canv, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
          ui.dialog('Molecule: ' + inputSequenceField.value)
            .add(canv)
            .showModal(true);
        });
        $(canvas).on('mouseover', () => $(canvas).css('cursor', 'zoom-in'));
        $(canvas).on('mouseout', () => $(canvas).css('cursor', 'default'));
        const mol = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, true,
        output.synthesizer![0]);
        // @ts-ignore
        OCL.StructureView.drawMolecule(canvas, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
        moleculeSvgDiv.append(canvas);
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
    updateTableAndMolecule(sequence, inputFormatChoiceInput.value!);
  });

  const asoDf = DG.DataFrame.fromObjects([
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
  ])!;
  const asoGapmersGrid = DG.Viewer.grid(asoDf, {showRowHeader: false, showCellTooltip: false});

  asoDf.onCurrentCellChanged.subscribe((_) => {
    navigator.clipboard.writeText(asoDf.currentCell.value).then(() => grok.shell.info('Copied'));
  });

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
      ])!, {showRowHeader: false, showCellTooltip: false},
  );

  const overhangModificationsGrid = DG.Viewer.grid(
      DG.DataFrame.fromColumns([
        DG.Column.fromStrings('Name', Object.keys(MODIFICATIONS)),
      ])!, {showRowHeader: false, showCellTooltip: false},
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

  const topPanel = [
    ui.iconFA('download', () => {
      const result = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''), false, false,
        inputFormatChoiceInput.value!);
      const element = document.createElement('a');
      element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
      element.setAttribute('download', inputSequenceField.value.replace(/\s/g, '') + '.mol');
      element.click();
    }, 'Save .mol file'),
    ui.iconFA('copy', () => {
      navigator.clipboard.writeText(
        sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''), false, inputFormatChoiceInput.value!))
        .then(() => grok.shell.info(sequenceWasCopied));
    }, 'Copy SMILES'),
    switchInput.root,
  ];

  const v = grok.shell.v;
  const tabControl = grok.shell.sidebar;
  tabControl.onTabChanged.subscribe((_) =>
    v.setRibbonPanels([(tabControl.currentPane.name == 'MAIN') ? topPanel : []]));
  v.setRibbonPanels([topPanel]);

  return ui.box(
    ui.splitH([
      ui.splitV([
        ui.panel([
          appMainDescription,
          ui.div([
            ui.h1('Input sequence'),
            ui.div([
              inputSequenceField.root,
            ], 'input-base'),
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
