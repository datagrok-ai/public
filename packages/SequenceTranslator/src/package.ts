/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import $ from 'cash-dom';
import {defineAxolabsPattern} from './defineAxolabsPattern';
import {saveSenseAntiSense} from './structures-works/save-sense-antisense';
import {sequenceToSmiles, sequenceToMolV3000} from './structures-works/from-monomers';
import {convertSequence, undefinedInputSequence, isValidSequence, getFormat} from
  './structures-works/sequence-codes-tools';
import {map, COL_NAMES, MODIFICATIONS} from './structures-works/map';
import {siRnaAxolabsToGcrs} from './structures-works/converters';
import {SALTS_CSV} from './salts';
import {USERS_CSV} from './users';
import {ICDS} from './ICDs';
import {SOURCES} from './sources';
import {IDPS} from './IDPs';

export const _package = new DG.Package();

const defaultInput = 'fAmCmGmAmCpsmU';
const sequenceWasCopied = 'Copied';
const tooltipSequence = 'Copy sequence';

//name: Sequence Translator
//tags: app
export function sequenceTranslator(): void {
  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  function updateTableAndMolecule(sequence: string, inputFormat: string, isSet: boolean): void {
    moleculeSvgDiv.innerHTML = '';
    outputTableDiv.innerHTML = '';
    const pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');
    let errorsExist = false;
    try {
      sequence = sequence.replace(/\s/g, '');
      const output = isValidSequence(sequence, null);
      if (isSet)
        output.synthesizer = [inputFormat];
      inputFormatChoiceInput.value = output.synthesizer![0];
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

  const inputFormatChoiceInput = ui.choiceInput(
    'Input format: ', 'Janssen GCRS Codes', Object.keys(map), (format: string) => {
      updateTableAndMolecule(inputSequenceField.value.replace(/\s/g, ''), format, true);
    });
  const moleculeSvgDiv = ui.block([]);
  const outputTableDiv = ui.div([]);
  const inputSequenceField = ui.textInput('', defaultInput, (sequence: string) => updateTableAndMolecule(sequence,
    inputFormatChoiceInput.value!, false));

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
  updateTableAndMolecule(defaultInput, inputFormatChoiceInput.value!, true);

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
    ),
    'AXOLABS': defineAxolabsPattern(),
    'SDF': saveSenseAntiSense(),
  });

  $(codesTablesDiv).hide();

  const v = grok.shell.newView('Sequence Translator', [tabControl]);
  v.box = true;

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

  tabControl.onTabChanged.subscribe((_) =>
    v.setRibbonPanels([(tabControl.currentPane.name == 'MAIN') ? topPanel : []]));
  v.setRibbonPanels([topPanel]);
}

async function saveTableAsSdFile(table: DG.DataFrame) {
  if (!table.columns.contains('Compound Name')) {
    grok.shell.warning(
      'File saved without columns \'' +
      [COL_NAMES.COMPOUND_NAME, COL_NAMES.COMPOUND_COMMENTS, COL_NAMES.CPD_MW,
        COL_NAMES.SALT_MASS, COL_NAMES.BATCH_MW].join('\', \''),
    );
  }
  const structureColumn = table.col(COL_NAMES.SEQUENCE)!;
  const typeColumn = table.col(COL_NAMES.TYPE)!;
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    const format = getFormat(structureColumn.get(i));
    result += (typeColumn.get(i) == 'SS') ?
      sequenceToMolV3000(structureColumn.get(i), false, true, format!) + '\n' + `>  <Sequence>\nSense Strand\n\n` :
      sequenceToMolV3000(structureColumn.get(i), true, true, format!) + '\n' + `>  <Sequence>\nAnti Sense\n\n`;
    for (const col of table.columns) {
      if (col.name != COL_NAMES.SEQUENCE)
        result += `>  <${col.name}>\n${col.get(i)}\n\n`;
    }
    result += '$$$$\n\n';
  }
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

const weightsObj: {[code: string]: number} = {};
for (const synthesizer of Object.keys(map)) {
  for (const technology of Object.keys(map[synthesizer])) {
    for (const code of Object.keys(map[synthesizer][technology]))
      weightsObj[code] ?? map[synthesizer][technology][code].weight;
  }
}
for (const [key, value] of Object.entries(MODIFICATIONS))
  weightsObj[key] = value.molecularWeight;

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}

function stringifyItems(items: string[]): string {
  return '["' + items.join('", "') + '"]';
}

function molecularWeight(sequence: string, weightsObj: {[index: string]: number}): number {
  const codes = sortByStringLengthInDescendingOrder(Object.keys(weightsObj)).concat(Object.keys(MODIFICATIONS));
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += weightsObj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode!.length;
  }
  return weight - 61.97;
}

//tags: autostart
export function autostartOligoSdFileSubscription() {
  let alreadyAdded = false;
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type == 'TableView') {
      if (v.dataFrame.columns.contains(COL_NAMES.TYPE))
        oligoSdFile(v.dataFrame);
      grok.events.onContextMenu.subscribe((args) => {
        for (const col of v.dataFrame.columns) {
          if (!alreadyAdded && DG.Detector.sampleCategories(col, (s) => /^[fsACGUacgu]{6,}$/.test(s))) {
            alreadyAdded = true;
            args.args.menu.item('Convert to GCRS', () => {
              const seqCol = args.args.context.table.currentCol;
              args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
                return siRnaAxolabsToGcrs(seqCol.get(i));
              });
            });
            break;
          };
        }
      });
    }
  });
}

export function oligoSdFile(table: DG.DataFrame) {
  const saltsDf = DG.DataFrame.fromCsv(SALTS_CSV);
  const usersDf = DG.DataFrame.fromCsv(USERS_CSV);
  const sourcesDf = DG.DataFrame.fromCsv(SOURCES);
  const icdsDf = DG.DataFrame.fromCsv(ICDS);
  const idpsDf = DG.DataFrame.fromCsv(IDPS);

  function addColumns(t: DG.DataFrame, saltsDf: DG.DataFrame) {
    if (t.columns.contains(COL_NAMES.COMPOUND_NAME))
      return grok.shell.error('Columns already exist!');

    const sequence = t.col(COL_NAMES.SEQUENCE)!;
    const salt = t.col(COL_NAMES.SALT)!;
    const equivalents = t.col(COL_NAMES.EQUIVALENTS)!;

    t.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => sequence.get(i));
    t.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => (i > 0 && i % 2 == 0) ?
      sequence.getString(i) + '; duplex of SS: ' + sequence.getString(i - 2) + ' and AS: ' + sequence.getString(i - 1) :
      sequence.getString(i),
    );
    const molWeightCol = saltsDf.col('MOLWEIGHT')!;
    const saltNamesList = saltsDf.col('DISPLAY')!.toList();
    t.columns.addNewFloat(COL_NAMES.CPD_MW)
      .init((i: number) => molecularWeight(sequence.get(i), weightsObj));
    t.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
      const saltRowIndex = saltNamesList.indexOf(salt.get(i));
      const mw = molWeightCol.get(saltRowIndex);
      return mw * equivalents.get(i);
    });
    t.columns.addNewCalculated(COL_NAMES.BATCH_MW,
      '${' + COL_NAMES.CPD_MW + '} + ${' + COL_NAMES.SALT_MASS + '}', DG.COLUMN_TYPE.FLOAT, false,
    );

    addColumnsPressed = true;
    return newDf = t;
  }

  let newDf: DG.DataFrame;
  let addColumnsPressed = false;

  const d = ui.div([
    ui.icons.edit(() => {
      d.innerHTML = '';
      if (table.col(COL_NAMES.IDP)!.type != DG.COLUMN_TYPE.STRING)
        table.changeColumnType(COL_NAMES.IDP, DG.COLUMN_TYPE.STRING);
      d.append(
        ui.link('Add Columns', () => {
          addColumns(table, saltsDf);
          grok.shell.tableView(table.name).grid.columns.setOrder(Object.values(COL_NAMES));
        }, 'Add columns: \'' + [COL_NAMES.COMPOUND_NAME, COL_NAMES.COMPOUND_COMMENTS, COL_NAMES.CPD_MW,
          COL_NAMES.SALT_MASS, COL_NAMES.BATCH_MW].join('\', \''), ''),
        ui.button('Save SD file', () => saveTableAsSdFile(addColumnsPressed ? newDf : table)),
      );

      const view = grok.shell.getTableView(table.name);

      view.table!.col(COL_NAMES.TYPE)!.setTag(DG.TAGS.CHOICES, '["AS", "SS", "Duplex"]');
      view.table!.col(COL_NAMES.OWNER)!.setTag(DG.TAGS.CHOICES, stringifyItems(usersDf.columns.byIndex(0).toList()));
      view.table!.col(COL_NAMES.SALT)!.setTag(DG.TAGS.CHOICES, stringifyItems(saltsDf.columns.byIndex(0).toList()));
      view.table!.col(COL_NAMES.SOURCE)!.setTag(DG.TAGS.CHOICES, stringifyItems(sourcesDf.columns.byIndex(0).toList()));
      view.table!.col(COL_NAMES.ICD)!.setTag(DG.TAGS.CHOICES, stringifyItems(icdsDf.columns.byIndex(0).toList()));
      view.table!.col(COL_NAMES.IDP)!.setTag(DG.TAGS.CHOICES, stringifyItems(idpsDf.columns.byIndex(0).toList()));

      grok.events.onContextMenu.subscribe((args) => {
        if ([COL_NAMES.TYPE, COL_NAMES.OWNER, COL_NAMES.SALT, COL_NAMES.SOURCE, COL_NAMES.ICD, COL_NAMES.IDP]
          .includes(args.args.context.table.currentCol.name)) {
          args.args.menu.item('Fill Column With Value', () => {
            const v = args.args.context.table.currentCell.value;
            args.args.context.table.currentCell.column.init(v);
          });
        }
      });
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
