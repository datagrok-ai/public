import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {weightsObj, individualBases, nearestNeighbour, SYNTHESIZERS} from './map';
import {validate} from './validation';
import {deleteWord, saveAsCsv, sortByStringLengthInDescOrder, mergeOptions, normalizeSequence} from './helpers';
import {COL_NAMES, CURRENT_USER, STORAGE_NAME, ADMIN_USERS, getAdditionalModifications, addModificationButton,
  deleteAdditionalModification} from './additional-modifications';

export const _package = new DG.Package();

const NAME_OF_COLUMN_WITH_SEQUENCES = 'Sequence';


//name: opticalDensity
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['NMole', 'Milligrams', 'Micrograms']}
//output: double opticalDensity
export async function opticalDensity(sequence: string, amount: number, outputUnits: string, extCoefsObj:
  {[index: string]: number}): Promise<number> {
  const ec = await extinctionCoefficient(sequence, extCoefsObj);
  if (outputUnits == 'Milligrams' || outputUnits == 'Micrograms' || outputUnits == 'mg' || outputUnits == 'µg')
    return (outputUnits == 'Milligrams' ? 1 : 0.001) * amount * ec / molecularWeight(sequence);
  if (outputUnits == 'OD')
    return amount;
  const coefficient = (outputUnits == 'NMole') ? 1000000 : (outputUnits == 'Milligrams') ? 1 : 1000;
  return amount * ec / coefficient;
}

//name: nMole
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micrograms']}
//output: double nMole
export async function nMole(sequence: string, amount: number, outputUnits: string, extinctionCoefficientsObj:
  {[index: string]: number}, weightsObj: {[index: string]: number}): Promise<number> {
  const ec = await extinctionCoefficient(sequence, extinctionCoefficientsObj);
  return (outputUnits == 'Optical Density') ?
    1000000 * amount / ec :
    1000 * amount / molecularWeight(sequence, weightsObj);
}

//name: molecularMass
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micromoles', 'Millimoles']}
//output: double molecularMass
export async function molecularMass(sequence: string, amount: number, outputUnits: string): Promise<number> {
  const additionalModificationsDf = await getAdditionalModifications();
  const additionalAbbreviations = additionalModificationsDf.col(COL_NAMES.ABBREVIATION)!.toList();
  const extinctionCoefficients = additionalModificationsDf.col(COL_NAMES.EXTINCTION_COEFFICIENT)!.toList();
  const additionalWeights = additionalModificationsDf.col(COL_NAMES.MOLECULAR_WEIGHT)!.toList();
  const extinctionCoefficientsObj: {[index: string]: number} = {};
  const additionalWeightsObj: {[index: string]: number} = {};
  additionalAbbreviations.forEach((key, i) => additionalWeightsObj[key] = additionalWeights[i]);
  additionalAbbreviations.forEach((key, i) => {
    if (extinctionCoefficients[i] != 'Base')
      extinctionCoefficientsObj[key] = extinctionCoefficients[i];
  });
  const ec = await extinctionCoefficient(sequence, extinctionCoefficientsObj);
  const od = await opticalDensity(sequence, amount, outputUnits, extinctionCoefficientsObj);
  const nm = await nMole(sequence, amount, outputUnits, extinctionCoefficientsObj, additionalWeightsObj);
  if (outputUnits == 'Optical Density' || outputUnits == 'OD') {
    return (ec == 0) ?
      amount * molecularWeight(sequence, additionalWeightsObj) :
      1000 * amount * molecularWeight(sequence, additionalWeightsObj) / ec;
  }
  const coefficient = (outputUnits == 'Milligrams' || outputUnits == 'Micromoles') ? 1 : 1000;
  return amount / ec * molecularWeight(sequence) * coefficient * od / nm;
}

export function molecularWeight(sequence: string, additionalWeightsObj?: {[index: string]: number}): number {
  const codes = (additionalWeightsObj == null) ?
    sortByStringLengthInDescOrder(Object.keys(weightsObj)) :
    sortByStringLengthInDescOrder(Object.keys(weightsObj).concat(Object.keys(additionalWeightsObj)));
  const obj = (additionalWeightsObj != null) ? mergeOptions(weightsObj, additionalWeightsObj) : weightsObj;
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += obj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode!.length;
  }
  return weight - 61.97;
}

export async function extinctionCoefficient(sequence: string, extCoefsObj?: {[i: string]: number}): Promise<number> {
  const additionalModificationsDf = await getAdditionalModifications();
  // const additionalCodes = additionalModificationsDf.col(COL_NAMES.ABBREVIATION)!.categories;
  // const output = isValidSequence(sequence, additionalCodes);
  let ns = normalizeSequence(sequence, SYNTHESIZERS.GCRS, null, additionalModificationsDf);
  let nearestNeighbourSum = 0;
  let individualBasisSum = 0;
  let modificationsSum = 0;
  if (extCoefsObj != null) {
    for (const modif of Object.keys(extCoefsObj)) {//@ts-ignore
      if (extCoefsObj[modif] != 'Base') {//@ts-ignore
        modificationsSum += (sequence.match(new RegExp(modif, 'g')) || []).length * parseFloat(extCoefsObj[modif]);
        ns = deleteWord(ns, modif);
      }
    }
  }
  for (let i = 0; i < ns.length - 2; i += 2) {
    nearestNeighbourSum += (ns[i] == ns[i + 2]) ?
      nearestNeighbour[ns.slice(i, i + 2)][ns.slice(i + 2, i + 4)] :
      (
        nearestNeighbour['r' + ((ns[i + 1] == 'T') ? 'U' : ns[i + 1])]['r' + ((ns[i + 3] == 'T') ? 'U' : ns[i + 3])] +
        nearestNeighbour['d' + ((ns[i + 1] == 'U') ? 'T' : ns[i + 1])]['d' + ((ns[i + 3] == 'U') ? 'T' : ns[i + 3])]
      ) / 2;
  }
  for (let i = 2; i < ns.length - 2; i += 2)
    individualBasisSum += individualBases[ns.slice(i, i + 2)];
  return nearestNeighbourSum - individualBasisSum + modificationsSum;
}

//name: Oligo Batch Calculator
//tags: app
export async function OligoBatchCalculatorApp(): Promise<void> {
  const additionalModsDf = await getAdditionalModifications();
  const additionalCodes = additionalModsDf.col(COL_NAMES.ABBREVIATION)!.categories;
  const additionalAbbreviations = additionalModsDf.col(COL_NAMES.ABBREVIATION)!.toList();
  const additionalWeights = additionalModsDf.col(COL_NAMES.MOLECULAR_WEIGHT)!.toList();
  const extinctionCoefficients = additionalModsDf.col(COL_NAMES.EXTINCTION_COEFFICIENT)!.toList();
  const additionalWeightsObj: {[index: string]: number} = {};
  const extinctionCoeffsObj: {[index: string]: number} = {};
  additionalAbbreviations.forEach((key, i) => additionalWeightsObj[key] = additionalWeights[i]);
  additionalAbbreviations.forEach((key, i) => {
    if (extinctionCoefficients[i] != 'Base')
      extinctionCoeffsObj[key] = (extinctionCoefficients[i] == null) ? 1 : extinctionCoefficients[i];
  });
  const mainGrid = DG.Viewer.grid(DG.DataFrame.create(), {'showRowHeader': false});

  async function render(text: string): Promise<void> {
    const sequences = text.split('\n')
      .map((s) => s.replace(/\s/g, ''))
      .filter((item) => item);

    const indicesOfFirstNotValidCharacter = Array(sequences.length);
    const normalizedSequences = Array(sequences.length);
    const molecularWeights = Array(sequences.length);
    const extinctionCoefficients = Array(sequences.length);
    const nMoles = Array(sequences.length);
    const opticalDensities = Array(sequences.length);
    const molecularMasses = Array(sequences.length);
    const reasonsOfError = Array(sequences.length);

    for (const [i, sequence] of sequences.entries()) {
      indicesOfFirstNotValidCharacter[i] = validate(sequence, additionalCodes);
      if (indicesOfFirstNotValidCharacter[i] < 0) {
        normalizedSequences[i] = normalizeSequence(sequence, SYNTHESIZERS.GCRS, null, additionalModsDf);
        if (normalizedSequences[i].length > 2) {
          try {
            molecularWeights[i] = molecularWeight(sequence, additionalWeightsObj);
            extinctionCoefficients[i] = await extinctionCoefficient(normalizedSequences[i], extinctionCoeffsObj);
            nMoles[i] = await nMole(sequence, yieldAmount.value, units.value, extinctionCoeffsObj,
              additionalWeightsObj);
            opticalDensities[i] = await opticalDensity(sequence, yieldAmount.value, units.value, extinctionCoeffsObj);
            molecularMasses[i] = await molecularMass(sequence, yieldAmount.value, units.value);
          } catch (e) {
            reasonsOfError[i] = 'Unknown error, please report it to Datagrok team';
            indicesOfFirstNotValidCharacter[i] = 0;
            grok.shell.error(String(e));
          }
        } else {
          reasonsOfError[i] = 'Sequence should contain at least two nucleotides';
          indicesOfFirstNotValidCharacter[i] = 0;
        }
      // }
      // else if (output.synthesizer == null)
      // reasonsOfError[i] = 'Not valid input';
      } else {
        reasonsOfError[i] = 'Sequence is expected to be in synthesizer \'' + SYNTHESIZERS.GCRS +
          '\', please see table below to see list of valid codes';
      }
    };

    const moleName1 = (units.value == 'µmole' || units.value == 'mg') ? 'µmole' : 'nmole';
    const moleName2 = (units.value == 'µmole') ? 'µmole' : 'nmole';
    const massName = (units.value == 'µmole') ? 'mg' : (units.value == 'mg') ? units.value : 'µg';
    const c = (units.value == 'mg' || units.value == 'µmole') ? 1000 : 1;

    mainGrid.dataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'Item', Array(...Array(sequences.length + 1).keys()).slice(1)),
      DG.Column.fromStrings(NAME_OF_COLUMN_WITH_SEQUENCES, sequences),
      DG.Column.fromList('int', 'Length', normalizedSequences.map((s) => s.length / 2)),
      DG.Column.fromList('double', 'OD 260', opticalDensities),
      DG.Column.fromList('double', moleName1, nMoles),
      DG.Column.fromList('double', 'Mass [' + massName + ']', molecularMasses),
      DG.Column.fromList('double', moleName2 + '/OD', nMoles.map(function(n, i) {return c * n / opticalDensities[i];})),
      DG.Column.fromList('double', 'µg/OD', molecularMasses.map(function(n, i) {return c * n / opticalDensities[i];})),
      DG.Column.fromList('double', 'MW', molecularWeights),
      DG.Column.fromList('double', 'Ext. Coefficient', extinctionCoefficients),
    ]);

    const col = mainGrid.col(NAME_OF_COLUMN_WITH_SEQUENCES)!;
    col.cellType = 'html';
    mainGrid.onCellPrepare(function(gc) {
      if (gc.isTableCell && gc.gridColumn.name == NAME_OF_COLUMN_WITH_SEQUENCES) {
        const items = (indicesOfFirstNotValidCharacter[gc.gridRow] < 0) ?
          [ui.divText(gc.cell.value, {style: {color: 'grey'}})] :
          [
            ui.divText(gc.cell.value.slice(0, indicesOfFirstNotValidCharacter[gc.gridRow]), {style: {color: 'grey'}}),
            ui.tooltip.bind(
              ui.divText(gc.cell.value.slice(indicesOfFirstNotValidCharacter[gc.gridRow]), {style: {color: 'red'}}),
              reasonsOfError[gc.gridRow],
            ),
          ];
        gc.style.element = ui.divH(items, {style: {margin: '6px 0 0 6px'}});
      }
    });
  }

  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const defaultInput = 'fAmCmGmAmCpsmU\nmApsmApsfGmAmUmCfGfAfC\nmAmUfGmGmUmCmAfAmGmA';

  const inputSequences = ui.textInput('', defaultInput, (txt: string) => render(txt));
  const yieldAmount = ui.floatInput('', 1, () => render(inputSequences.value));
  const units = ui.choiceInput('', 'OD', ['OD', 'µg', 'mg', 'µmole', 'nmole'], () => render(inputSequences.value));

  await render(defaultInput);

  const title = ui.panel([ui.h2('Oligo Properties')], 'ui-panel ui-box');
  title.style.maxHeight = '40px';
  $(title).children('h2').css('margin', '0px');

  const asoGapmersGrid = DG.Viewer.grid(
    DG.DataFrame.fromObjects([
      {'Name': '2\'MOE-5Me-rU', 'BioSpring': '5', 'Janssen GCRS': 'moeT', 'Weight': 378.27},
      {'Name': '2\'MOE-rA', 'BioSpring': '6', 'Janssen GCRS': 'moeA', 'Weight': 387.29},
      {'Name': '2\'MOE-5Me-rC', 'BioSpring': '7', 'Janssen GCRS': 'moe5mC', 'Weight': 377.29},
      {'Name': '2\'MOE-rG', 'BioSpring': '8', 'Janssen GCRS': 'moeG', 'Weight': 403.28},
      {'Name': '5-Methyl-dC', 'BioSpring': '9', 'Janssen GCRS': '5mC', 'Weight': 303.21},
      {'Name': 'ps linkage', 'BioSpring': '*', 'Janssen GCRS': 'ps', 'Weight': 16.07},
      {'Name': 'dA', 'BioSpring': 'A', 'Janssen GCRS': 'A, dA', 'Weight': 313.21},
      {'Name': 'dC', 'BioSpring': 'C', 'Janssen GCRS': 'C, dC', 'Weight': 289.18},
      {'Name': 'dG', 'BioSpring': 'G', 'Janssen GCRS': 'G, dG', 'Weight': 329.21},
      {'Name': 'dT', 'BioSpring': 'T', 'Janssen GCRS': 'T, dT', 'Weight': 304.2},
      {'Name': 'rA', 'BioSpring': '', 'Janssen GCRS': 'rA', 'Weight': 329.21},
      {'Name': 'rC', 'BioSpring': '', 'Janssen GCRS': 'rC', 'Weight': 305.18},
      {'Name': 'rG', 'BioSpring': '', 'Janssen GCRS': 'rG', 'Weight': 345.21},
      {'Name': 'rU', 'BioSpring': '', 'Janssen GCRS': 'rU', 'Weight': 306.17},
    ])!, {showRowHeader: false, showCellTooltip: false},
  );

  const omeAndFluoroGrid = DG.Viewer.grid(
    DG.DataFrame.fromObjects([
      {'Name': '2\'-fluoro-U', 'BioSpring': '1', 'Axolabs': 'Uf', 'Janssen GCRS': 'fU', 'Weight': 308.16},
      {'Name': '2\'-fluoro-A', 'BioSpring': '2', 'Axolabs': 'Af', 'Janssen GCRS': 'fA', 'Weight': 331.2},
      {'Name': '2\'-fluoro-C', 'BioSpring': '3', 'Axolabs': 'Cf', 'Janssen GCRS': 'fC', 'Weight': 307.18},
      {'Name': '2\'-fluoro-G', 'BioSpring': '4', 'Axolabs': 'Gf', 'Janssen GCRS': 'fG', 'Weight': 347.19},
      {'Name': '2\'OMe-rU', 'BioSpring': '5', 'Axolabs': 'u', 'Janssen GCRS': 'mU', 'Weight': 320.2},
      {'Name': '2\'OMe-rA', 'BioSpring': '6', 'Axolabs': 'a', 'Janssen GCRS': 'mA', 'Weight': 343.24},
      {'Name': '2\'OMe-rC', 'BioSpring': '7', 'Axolabs': 'c', 'Janssen GCRS': 'mC', 'Weight': 319.21},
      {'Name': '2\'OMe-rG', 'BioSpring': '8', 'Axolabs': 'g', 'Janssen GCRS': 'mG', 'Weight': 359.24},
      {'Name': 'ps linkage', 'BioSpring': '*', 'Axolabs': 's', 'Janssen GCRS': 'ps', 'Weight': 16.07},
    ])!, {showRowHeader: false, showCellTooltip: false},
  );

  const additionaModifsGrid = DG.Viewer.grid(additionalModsDf, {showRowHeader: false, showCellTooltip: true});
  additionaModifsGrid.col(COL_NAMES.LONG_NAMES)!.width = 110;
  additionaModifsGrid.col(COL_NAMES.ABBREVIATION)!.width = 80;
  additionaModifsGrid.col(COL_NAMES.MOLECULAR_WEIGHT)!.width = 105;
  additionaModifsGrid.col(COL_NAMES.BASE_MODIFICATION)!.width = 110;
  additionaModifsGrid.col(COL_NAMES.EXTINCTION_COEFFICIENT)!.width = 100;

  // Hide 'CHANGE_LOGS' column, display its content in tooltip
  additionaModifsGrid.columns.setVisible(additionalModsDf.columns.names().slice(0, -1));
  additionalModsDf.col(COL_NAMES.CHANGE_LOGS)!.name = '~' + COL_NAMES.CHANGE_LOGS;
  additionaModifsGrid.onCellTooltip(function(cell, x, y) {
    if (cell.isTableCell) {
      const v = additionalModsDf.col('~' + COL_NAMES.CHANGE_LOGS)!.get(cell.gridRow).split('; ').slice(0, -1);
      ui.tooltip.show(ui.divText(v), x, y);
      return true;
    }
  });

  const codesTablesDiv = ui.splitV([
    ui.box(ui.h2('Additional modifications'), {style: {maxHeight: '40px'}}),
    additionaModifsGrid.root,
    ui.box(ui.h2('ASO Gapmers'), {style: {maxHeight: '40px'}}),
    asoGapmersGrid.root,
    ui.box(ui.h2('2\'-OMe and 2\'-F modifications'), {style: {maxHeight: '40px'}}),
    omeAndFluoroGrid.root,
  ], {style: {maxWidth: '600px'}});

  additionalModsDf.col(COL_NAMES.BASE_MODIFICATION)!
    .setTag(DG.TAGS.CHOICES, '["NO", "rU", "rA", "rC", "rG", "dA", "dC", "dG", "dT"]');

  const view = grok.shell.newView('Oligo Batch Calculator', [
    ui.splitH([
      ui.splitV([
        ui.box(
          ui.panel([
            ui.h2('Yield Amount & Units'),
            ui.divH([
              yieldAmount.root,
              units.root,
            ]),
            ui.h2('Input Sequences'),
            ui.div([
              inputSequences.root,
            ], 'inputSequences'),
          ]), {style: {maxHeight: '230px'}},
        ),
        ui.splitV([
          title,
          mainGrid.root,
        ]),
      ]),
      codesTablesDiv,
    ]),
  ]);
  view.box = true;
  view.path = '/apps/OligoBatchCalculator/';
  view.setRibbonPanels([[
    ui.iconFA('redo', () => inputSequences.value = ''),
    ui.iconFA('plus', () => addModificationButton(additionalModsDf)),
    ui.iconFA('arrow-to-bottom', () => saveAsCsv(mainGrid.dataFrame!)),
    ui.switchInput('Codes', true, (v: boolean) => (v) ? $(codesTablesDiv).show() : $(codesTablesDiv).hide()).root,
  ]]);

  const col = additionaModifsGrid.col(COL_NAMES.ACTION)!;
  col.cellType = 'html';
  additionaModifsGrid.onCellPrepare(function(gc) {
    if (gc.isTableCell && gc.gridColumn.name == COL_NAMES.ACTION) {
      const icon = ui.iconFA('trash-alt');
      gc.style.element = ui.button(icon, () => deleteAdditionalModification(additionalModsDf, gc.gridRow));
    }
  });

  let tempValue = '';
  additionalModsDf.onCurrentCellChanged.subscribe(() => {
    tempValue = additionalModsDf.currentCell.value;
  });

  DG.debounce(additionalModsDf.onValuesChanged, 10).subscribe(async (_) => {
    grok.dapi.users.current().then((user) => {
      if (!ADMIN_USERS.includes(user.firstName + ' ' + user.lastName))
        return grok.shell.warning('You don\'t have permission for this action');
    });
    if (additionalModsDf.currentCol.name == COL_NAMES.ABBREVIATION) {
      const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
      if (additionalModsDf.currentCell.value.length > 100)
        return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
      if (additionalModsDf.currentCell.value in entries) {
        additionalModsDf.set(additionalModsDf.currentCol.name, additionalModsDf.currentRowIdx, tempValue);
        return grok.shell.warning('Abbreviation ' + additionalModsDf.currentCell.value + ' already exists');
      }
    }
    if (additionalModsDf.currentCol.name == COL_NAMES.LONG_NAMES && additionalModsDf.currentCell.value.length > 300)
      return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');

    const rowIndex = additionalModsDf.currentCell.rowIndex;
    if (additionalModsDf.currentCol.name == COL_NAMES.BASE_MODIFICATION) {
      if (additionalModsDf.currentCell.value == 'NO') {
        const extCoefChoiceInput = ui.floatInput('', 0);
        ui.dialog('Enter Extinction Coefficient Value')
          .add(extCoefChoiceInput)
          .onOK(() => {
            const col = additionalModsDf.getCol(COL_NAMES.EXTINCTION_COEFFICIENT);
            col.set(rowIndex, String(extCoefChoiceInput.value), false);
            additionaModifsGrid.invalidate();
          })
          .show();
      } else {
        const col = additionalModsDf.getCol(COL_NAMES.EXTINCTION_COEFFICIENT);
        col.set(rowIndex, 'Base', false);
        additionaModifsGrid.invalidate();
      }
    }

    await grok.dapi.userDataStorage.postValue(
      STORAGE_NAME,
      additionalModsDf.col(COL_NAMES.ABBREVIATION)!.get(rowIndex),
      JSON.stringify({
        longName: additionalModsDf.col(COL_NAMES.LONG_NAMES)!.get(rowIndex),
        abbreviation: additionalModsDf.col(COL_NAMES.ABBREVIATION)!.get(rowIndex),
        molecularWeight: additionalModsDf.col(COL_NAMES.MOLECULAR_WEIGHT)!.get(rowIndex),
        extinctionCoefficient: additionalModsDf.col(COL_NAMES.EXTINCTION_COEFFICIENT)!.get(rowIndex),
        baseModification: additionalModsDf.col(COL_NAMES.BASE_MODIFICATION)!.get(rowIndex),
        changeLogs: additionalModsDf.col('~' + COL_NAMES.CHANGE_LOGS)!.get(rowIndex),
      }),
      CURRENT_USER,
    );
  });
}
