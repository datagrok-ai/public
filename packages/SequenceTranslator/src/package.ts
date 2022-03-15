/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import $ from 'cash-dom';
import {defineAxolabsPattern} from './defineAxolabsPattern';
import {saveSenseAntiSense} from './structures-works/save-sense-antisense';
import {sequenceToSmiles, sequenceToMolV3000} from './structures-works/from-monomers';
import {convertSequence, undefinedInputSequence} from './structures-works/sequence-codes-tools';
import {SALTS_CSV} from './salts';
import {USERS_CSV} from './users';
import {ICDS} from './ICDs';
import {SOURCES} from './sources';

export const _package = new DG.Package();

const defaultInput = 'A';//'AGGTCCTCTTGACTTAGGCC';
const sequenceWasCopied = 'Copied';
const tooltipSequence = 'Copy sequence';

//name: Sequence Translator
//tags: app
export function sequenceTranslator(): void {
  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  function updateTableAndMolecule(sequence: string): void {
    moleculeSvgDiv.innerHTML = '';
    outputTableDiv.innerHTML = '';
    const pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');
    let errorsExist = false;
    try {
      const outputSequenceObj = convertSequence(sequence);
      const tableRows = [];

      for (const key of Object.keys(outputSequenceObj).slice(1)) {
        const indexOfFirstNotValidCharacter = ('indexOfFirstNotValidCharacter' in outputSequenceObj) ?
          JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).indexOfFirstNotValidCharacter :
          -1;
        if ('indexOfFirstNotValidCharacter' in outputSequenceObj) {
          const indexOfFirstNotValidCharacter = ('indexOfFirstNotValidCharacter' in outputSequenceObj) ?
            JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).indexOfFirstNotValidCharacter :
            -1;
          if (indexOfFirstNotValidCharacter != -1)
            errorsExist = true;
        }

        tableRows.push({
          'key': key,
          'value': ('indexOfFirstNotValidCharacter' in outputSequenceObj) ?
            ui.divH([
              ui.divText(sequence.slice(0, indexOfFirstNotValidCharacter), {style: {color: 'grey'}}),
              ui.tooltip.bind(
                ui.divText(sequence.slice(indexOfFirstNotValidCharacter), {style: {color: 'red'}}),
                'Expected format: ' + JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).expectedSynthesizer +
                '. See tables with valid codes on the right',
              ),
            ]) : //@ts-ignore
            ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key])
              .then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, ''),
        });
      }

      if (errorsExist) {
        const expectedSynthesizer = JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!)
          .expectedSynthesizer.slice(0, -6);
        asoGapmersGrid.onCellPrepare(function(gc) {
          gc.style.backColor = (gc.gridColumn.name == expectedSynthesizer) ? 0xFFF00000 : 0xFFFFFFFF;
        });
        omeAndFluoroGrid.onCellPrepare(function(gc) {
          gc.style.backColor = (gc.gridColumn.name == expectedSynthesizer) ? 0xFFF00000 : 0xFFFFFFFF;
        });
        switchInput.enabled = true;
      } else {
        asoGapmersGrid.onCellPrepare(function(gc) {
          gc.style.backColor = 0xFFFFFFFF;
        });
        omeAndFluoroGrid.onCellPrepare(function(gc) {
          gc.style.backColor = 0xFFFFFFFF;
        });
      }

      outputTableDiv.append(
        ui.div([
          DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) =>
            [item.key, item.value], ['Code', 'Sequence']).root,
        ], 'table'),
      );
      semTypeOfInputSequence.textContent = 'Detected input type: ' + outputSequenceObj.type;

      if (outputSequenceObj.type != undefinedInputSequence && outputSequenceObj.Error != undefinedInputSequence) {
        const canvas = ui.canvas(300, 170);
        canvas.addEventListener('click', () => {
          const canv = ui.canvas($(window).width(), $(window).height());
          const mol = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''));
          // @ts-ignore
          OCL.StructureView.drawMolecule(canv, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
          ui.dialog('Molecule: ' + inputSequenceField.value)
            .add(canv)
            .showModal(true);
        });
        $(canvas).on('mouseover', () => $(canvas).css('cursor', 'zoom-in'));
        $(canvas).on('mouseout', () => $(canvas).css('cursor', 'default'));
        const mol = sequenceToMolV3000(inputSequenceField.value.replace(/\s/g, ''));
        // @ts-ignore
        OCL.StructureView.drawMolecule(canvas, OCL.Molecule.fromMolfile(mol), {suppressChiralText: true});
        moleculeSvgDiv.append(canvas);
      } else
        moleculeSvgDiv.innerHTML = '';
    } finally {
      pi.close();
    }
  }

  const semTypeOfInputSequence = ui.divText('');
  const moleculeSvgDiv = ui.block([]);
  const outputTableDiv = ui.div([], 'table');
  const inputSequenceField = ui.textInput('', defaultInput, (sequence: string) => updateTableAndMolecule(sequence));

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
    DG.DataFrame.fromObjects([
      {'Name': '(invabasic)'},
      {'Name': '(GalNAc-2-JNJ)'},
    ])!, {showRowHeader: false, showCellTooltip: false},
  );
  updateTableAndMolecule(defaultInput);

  const appMainDescription = ui.info([
    ui.divText('How to convert one sequence:', {style: {'font-weight': 'bolder'}}),
    ui.divText('Paste sequence into the text field below'),
    ui.divText('\n How to convert many sequences:', {style: {'font-weight': 'bolder'}}),
    ui.divText('1. Drag & drop an Excel or CSV file with sequences into Datagrok'),
    ui.divText('2. Right-click on the column header, then see the \'Convert\' menu'),
    ui.divText('This will add the result column to the right of the table'),
  ], 'Convert oligonucleotide sequences between Nucleotides, BioSpring, Axolabs, Mermade 12 and GCRS representations.');

  const codesTablesDiv = ui.splitV([
    ui.box(ui.h2('ASO Gapmers'), {style: {maxHeight: '40px'}}),
    asoGapmersGrid.root,
    ui.box(ui.h2('2\'-OMe and 2\'-F siRNA'), {style: {maxHeight: '40px'}}),
    omeAndFluoroGrid.root,
    ui.box(ui.h2('Overhang modifications'), {style: {maxHeight: '40px'}}),
    overhangModificationsGrid.root,
  ], {style: {maxWidth: '350px'}});

  const tabControl = ui.tabControl({
    'MAIN': ui.box(
      ui.splitH([
        ui.splitV([
          ui.panel([
            appMainDescription,
            ui.div([
              ui.h1('Input sequence'),
              ui.div([
                inputSequenceField.root,
              ], 'input-base'),
            ], 'sequenceInput'),
            semTypeOfInputSequence,
            ui.block([
              ui.h1('Output'),
              outputTableDiv,
            ]),
            moleculeSvgDiv,
          ], 'sequence'),
        ]),
        codesTablesDiv,
      ], {style: {height: '100%', width: '100%'}}),
    ),
    'AXOLABS': defineAxolabsPattern(),
    'SDF': saveSenseAntiSense(),
  });

  const v = grok.shell.newView('Sequence Translator', [tabControl]);
  v.box = true;

  const switchInput = ui.switchInput('Codes', true, (v: boolean) => (v) ?
    $(codesTablesDiv).show() :
    $(codesTablesDiv).hide(),
  );

  const topPanel = [
    ui.iconFA('download', () => {
      const smiles = sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''));
      const result = `${OCL.Molecule.fromSmiles(smiles).toMolfile()}\n`;
      const element = document.createElement('a');
      element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
      element.setAttribute('download', inputSequenceField.value.replace(/\s/g, '') + '.mol');
      element.click();
    }, 'Save .mol file'),
    ui.iconFA('copy', () => {
      navigator.clipboard.writeText(sequenceToSmiles(inputSequenceField.value.replace(/\s/g, '')))
        .then(() => grok.shell.info(sequenceWasCopied));
    }, 'Copy SMILES'),
    switchInput.root,
  ];

  tabControl.onTabChanged.subscribe((_) =>
    v.setRibbonPanels([(tabControl.currentPane.name == 'MAIN') ? topPanel : []]));
  v.setRibbonPanels([topPanel]);

  $('.sequence')
    .children().css('padding', '5px 0');
  $('.sequenceInput .input-base')
    .css('margin', '0');
  $('.sequenceInput textarea')
    .css('resize', 'none')
    .css('min-height', '50px')
    .css('width', '100%')
    .attr('spellcheck', 'false');
  $('.sequenceInput select')
    .css('width', '100%');
}

async function saveTableAsSdFile(table: DG.DataFrame) {
  if (!table.columns.contains('Compound Name')) {
    grok.shell.warning(
      'File saved without columns \'Compound Name\', \'Compound Components\', \'Cpd MW\', \'Salt mass\', \'Batch MW\'');
  }
  const structureColumn = table.columns.byName('Sequence');
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    try {
      const smiles = sequenceToSmiles(structureColumn.get(i));
      const mol = OCL.Molecule.fromSmiles(smiles);
      result += `\n${mol.toMolfile()}\n`;
      for (const col of table.columns)
        result += `>  <${col.name}>\n${col.get(i)}\n\n`;
      result += '$$$$';
    } catch (error) {
      console.error(error);
    }
  }
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

//tags: autostart
export function autostartOligoSdFileSubscription() {
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type == 'TableView' && v.dataFrame.columns.contains('Type'))
      oligoSdFile(v.dataFrame);
  });
}

export function oligoSdFile(table: DG.DataFrame) {
  const saltsDf = DG.DataFrame.fromCsv(SALTS_CSV);
  const usersDf = DG.DataFrame.fromCsv(USERS_CSV);
  const sourcesDf = DG.DataFrame.fromCsv(SOURCES);
  const icdsDf = DG.DataFrame.fromCsv(ICDS);
  function addColumns(t: DG.DataFrame, saltsDf: DG.DataFrame) {
    if (t.columns.contains('Compound Name'))
      return grok.shell.error('Columns already exist!');

    table.col('Source')?.init('Johnson and Johnson Pharma');
    table.col('ICD')?.init('No Contract');

    const sequence = t.col('Sequence')!;
    const salt = t.col('Salt')!;
    const equivalents = t.col('Equivalents')!;

    t.columns.addNewString('Compound Name').init((i: number) => sequence.get(i));
    t.columns.addNewString('Compound Comments').init((i: number) => (i > 0 && i % 2 == 0) ?
      sequence.getString(i) + '; duplex of SS: ' + sequence.getString(i - 2) + ' and AS: ' + sequence.getString(i - 1) :
      sequence.getString(i),
    );
    const chargeCol = saltsDf.col('CHARGE')!.toList();
    const saltNames = saltsDf.col('DISPLAY')!.toList();
    const molWeight = saltsDf.col('MOLWEIGHT')!.toList();
    t.columns.addNewFloat('Cpd MW').init((i: number) => ((i + 1) % 3 == 0) ? DG.FLOAT_NULL : molWeight[i]);
    t.columns.addNewFloat('Salt mass').init((i: number) => {
      const v = chargeCol[saltNames.indexOf(salt.get(i))];
      const n = (v == null) ? 0 : chargeCol[saltNames.indexOf(salt.get(i))];
      return n * equivalents.get(i);
    });
    t.columns.addNewCalculated('Batch MW', '${Cpd MW} + ${Salt mass}', DG.COLUMN_TYPE.FLOAT, false);

    addColumnsPressed = true;
    return newDf = t;
  }

  const columnsOrder = ['Chemistry', 'Number', 'Type', 'Chemistry Name', 'Internal compound ID',
    'IDP', 'Sequence', 'Compound Name', 'Compound Comments', 'Salt', 'Equivalents', 'Purity', 'Cpd MW', 'Salt mass',
    'Batch MW', 'Source', 'ICD', 'Owner'];
  let newDf: DG.DataFrame;
  let addColumnsPressed = false;

  const d = ui.div([
    ui.icons.edit(() => {
      d.innerHTML = '';
      d.append(
        ui.link('Add Columns', () => {
          addColumns(table, saltsDf);
          grok.shell.tableView(table.name).grid.columns.setOrder(columnsOrder);
        }, 'Add columns: Compound Name, Compound Components, Cpd MW, Salt mass, Batch MW', ''),
        ui.button('Save SD file', () => saveTableAsSdFile(addColumnsPressed ? newDf : table)),
      );
      const view = grok.shell.getTableView(table.name);
      const saltCol = view.grid.col('Salt')!;
      const typeCol = view.grid.col('Type')!;
      const ownerCol = view.grid.col('Owner')!;
      const sourcesCol = view.grid.col('Source')!;
      const icdsCol = view.grid.col('ICD')!;
      saltCol.cellType = 'html';
      typeCol.cellType = 'html';
      ownerCol.cellType = 'html';
      sourcesCol.cellType = 'html';
      icdsCol.cellType = 'html';
      view.grid.onCellPrepare(function(gc: DG.GridCell) {
        console.log('Start');
        if (gc.isTableCell) {
          console.log(gc.gridColumn.name);
          if (gc.gridColumn.name == 'Type')
            gc.style.element = ui.choiceInput('', gc.cell.value, ['AS', 'SS', 'Duplex']).root;
          else if (gc.gridColumn.name == 'Owner') {
            gc.style.element = ui.choiceInput('', gc.cell.value, usersDf.columns.byIndex(0).toList(), () => {
              view.dataFrame.col('Owner')!.set(gc.gridRow, '');
            }).root;
          } else if (gc.gridColumn.name == 'Salt') {
            gc.style.element = ui.choiceInput('', gc.cell.value, saltsDf.columns.byIndex(1).toList(), () => {
              view.dataFrame.col('Salt')!.set(gc.gridRow, '');
            }).root;
          } else if (gc.gridColumn.name == 'Source') {
            gc.style.element = ui.choiceInput('', gc.cell.value, sourcesDf.columns.byIndex(0).toList(), () => {
              view.dataFrame.col('Source')!.set(gc.gridRow, '');
            }).root;
          } else if (gc.gridColumn.name == 'ICD') {
            gc.style.element = ui.choiceInput('', gc.cell.value, icdsDf.columns.byIndex(0).toList(), () => {
              view.dataFrame.col('ICD')!.set(gc.gridRow, '');
            }).root;
          }
        }
      });

      table.onDataChanged.subscribe((_) => {
        if (table.currentCol.name == 'IDP' && typeof table.currentCell.value != 'number')
          grok.shell.error('Value should be numeric');
      });
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
