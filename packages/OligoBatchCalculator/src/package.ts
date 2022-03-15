import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {map, individualBases, nearestNeighbour} from './map';
import {isValidSequence, getAllCodesOfSynthesizer} from './validation';

export const _package = new DG.Package();

const CURRENT_USER = false;
const STORAGE_NAME = 'oligo-batch-calculator-storage';
const OVERHANG_COL_NAMES = {
  LONG_NAMES: 'Long name',
  ABBREVIATION: 'Abbreviation',
  MOLECULAR_WEIGHT: 'Molecular weight',
  BASE_MODIFICATION: 'Base modification',
  EXTINCTION_COEFFICIENT: 'Ext. coefficient',
  ACTION: 'Action',
  CHANGE_LOGS: 'Change logs',
};
const ADMIN_USERS = ['Baozhong Zhao', 'Sijin Guo', 'Saika Siddiqui', 'Vadym Kovadlo'];
const NAME_OF_COLUMN_WITH_SEQUENCES = 'Sequence';

let weightsObj: {[code: string]: number} = {};
const normalizedObj: {[code: string]: string} = {};
for (const synthesizer of Object.keys(map)) {
  for (const technology of Object.keys(map[synthesizer])) {
    for (const code of Object.keys(map[synthesizer][technology])) {
      weightsObj[code] = map[synthesizer][technology][code].weight;
      normalizedObj[code] = map[synthesizer][technology][code].normalized;
    }
  }
}

function saveAsCsv(table: DG.DataFrame): void {
  const link = document.createElement('a');
  link.setAttribute('href', 'data:text/csv;charset=utf-8,\uFEFF' + encodeURI(table.toCsv()));
  link.setAttribute('download', 'Oligo Properties.csv');
  link.click();
}

export function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}

export async function getOverhangModificationsDf(): Promise<DG.DataFrame> {
  const modifications: any[] = [];
  const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
  if (entries != null && Object.keys(entries).length == 0)
    grok.shell.info('Storage is empty. Try to post something to the storage');
  else {
    for (const key of Object.keys(entries))
      modifications.push(JSON.parse(entries[key]));
  }
  const molWeightList = modifications.map((e) => (e.molecularWeight == undefined) ? 0 : e.molecularWeight);
  const extinctionCoefList = modifications.map((e) => String(e.extinctionCoefficient));
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(OVERHANG_COL_NAMES.LONG_NAMES, modifications.map((e) => e.longName)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.ABBREVIATION, modifications.map((e) => e.abbreviation)), // @ts-ignore
    DG.Column.fromFloat32Array(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, molWeightList),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.BASE_MODIFICATION, modifications.map((e) => e.baseModification)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, extinctionCoefList),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.ACTION, Array(modifications.length)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.CHANGE_LOGS, modifications.map((e) => e.changeLogs)),
  ])!;
}

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
  {[index: string]: number}): Promise<number> {
  const ec = await extinctionCoefficient(sequence, extinctionCoefficientsObj);
  return (outputUnits == 'Optical Density') ? 1000000 * amount / ec : 1000 * amount / molecularWeight(sequence);
}

//name: molecularMass
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micromoles', 'Millimoles']}
//output: double molecularMass
export async function molecularMass(sequence: string, amount: number, outputUnits: string): Promise<number> {
  const overhangModificationsDf = await getOverhangModificationsDf();
  const overhangsAbbreviations = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.toList();
  const extinctionCoefficients = overhangModificationsDf.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.toList();
  const extinctionCoefficientsObj: {[index: string]: number} = {};
  overhangsAbbreviations.forEach((key, i) => extinctionCoefficientsObj[key] = extinctionCoefficients[i]);
  const ec = await extinctionCoefficient(sequence, extinctionCoefficientsObj);
  const od = await opticalDensity(sequence, amount, outputUnits, extinctionCoefficientsObj);
  const nm = await nMole(sequence, amount, outputUnits, extinctionCoefficientsObj);
  if (outputUnits == 'Optical Density' || outputUnits == 'OD')
    return (ec == 0) ? amount * molecularWeight(sequence) : 1000 * amount * molecularWeight(sequence) / ec;
  const coefficient = (outputUnits == 'Milligrams' || outputUnits == 'Micromoles') ? 1 : 1000;
  return amount / ec * molecularWeight(sequence) * coefficient * od / nm;
}

//name: molecularWeight
//input: string sequence
//output: double molecularWeight
export function molecularWeight(sequence: string, overhangsWeightsObj?: {[index: string]: number}): number {
  const codes = (overhangsWeightsObj == null) ?
    sortByStringLengthInDescendingOrder(Object.keys(weightsObj)) :
    sortByStringLengthInDescendingOrder(Object.keys(weightsObj).concat(Object.keys(overhangsWeightsObj)));
  if (overhangsWeightsObj != null)
    weightsObj = mergeOptions(weightsObj, overhangsWeightsObj);
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += weightsObj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode!.length;
  }
  return weight - 61.97;
}

function deleteWord(sequence: string, searchTerm: string): string {
  let n = sequence.search(searchTerm);
  while (sequence.search(searchTerm) > -1) {
    n = sequence.search(searchTerm);
    sequence = sequence.substring(0, n - 1) + sequence.substring(n + searchTerm.length - 1, sequence.length);
  }
  return sequence;
}

//name: extinctionCoefficient
//input: string sequence
//output: double extinctionCoefficient
export async function extinctionCoefficient(sequence: string, extCoefsObj?: {[i: string]: number}): Promise<number> {
  const overhangModificationsDf = await getOverhangModificationsDf();
  const overhangCodes = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.categories;
  const output = isValidSequence(sequence, overhangCodes);
  sequence = normalizeSequence(sequence, output.expectedSynthesizer, output.expectedTechnology);
  let nearestNeighbourSum = 0;
  let individualBasisSum = 0;
  let modificationsSum = 0;
  if (extCoefsObj != null) {
    for (const modif of Object.keys(extCoefsObj)) {
      modificationsSum += (sequence.match(new RegExp(modif, 'g')) || []).length * extCoefsObj[modif];
      sequence = deleteWord(sequence, modif);
    }
  }
  for (let i = 0; i < sequence.length - 2; i += 2) {
    nearestNeighbourSum += (sequence[i] == sequence[i + 2]) ?
      nearestNeighbour[sequence.slice(i, i + 2)][sequence.slice(i + 2, i + 4)] :
      (
        nearestNeighbour['r' + ((sequence[i + 1] == 'T') ? 'U' : sequence[i + 1])]['r' + ((sequence[i + 3] == 'T') ?
          'U' : sequence[i + 3])] +
        nearestNeighbour['d' + ((sequence[i + 1] == 'U') ? 'T' : sequence[i + 1])]['d' + ((sequence[i + 3] == 'U') ?
          'T' : sequence[i + 3])]
      ) / 2;
  }
  for (let i = 2; i < sequence.length - 2; i += 2)
    individualBasisSum += individualBases[sequence.slice(i, i + 2)];
  return nearestNeighbourSum - individualBasisSum + modificationsSum;
}

function normalizeSequence(sequence: string, synthesizer: string | null, technology: string | null): string {
  const codes = (technology == null) ?
    getAllCodesOfSynthesizer(synthesizer!) :
    Object.keys(map[synthesizer!][technology]);
  const sortedCodes = sortByStringLengthInDescendingOrder(codes);
  const regExp = new RegExp('(' + sortedCodes.join('|') + ')', 'g');
  return sequence.replace(regExp, function(code) {return normalizedObj[code];});
}

async function addModificationButton(modificationsDf: DG.DataFrame): Promise<void> {
  grok.dapi.users.current().then((user) => {
    if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
      const longName = ui.stringInput(OVERHANG_COL_NAMES.LONG_NAMES, '');
      ui.tooltip.bind(longName.root, 'Examples: \'Inverted Abasic\', \'Cyanine 3 CPG\', \'5-Methyl dC\'');
      const abbreviation = ui.stringInput(OVERHANG_COL_NAMES.ABBREVIATION, '');
      ui.tooltip.bind(abbreviation.root, 'Examples: \'invabasic\', \'Cy3\', \'5MedC\'');
      const molecularWeight = ui.floatInput(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, 0);
      const baseModification = ui.choiceInput(OVERHANG_COL_NAMES.BASE_MODIFICATION,
        'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
          if (v != 'NO')
            extCoefficient.value = 'Base';
          extCoefficient.enabled = (v == 'NO');
        });
      const extCoefficient = ui.stringInput(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, '');
      ui.dialog('Add Modification')
        .add(ui.block([
          longName.root,
          abbreviation.root,
          molecularWeight.root,
          baseModification.root,
          extCoefficient.root,
        ]))
        .onOK(async () => {
          if (longName.value.length > 300)
            return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');
          if (abbreviation.value.length > 100)
            return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
          const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
          if (abbreviation.value in entries)
            return grok.shell.warning('Abbreviation ' + abbreviation.value + ' already exists');
          const modifiedLogs = Date() + ' by ' + user.firstName + ' ' + user.lastName + '; ';
          grok.dapi.userDataStorage.postValue(
            STORAGE_NAME,
            abbreviation.value,
            JSON.stringify({
              longName: longName.value,
              abbreviation: abbreviation.value,
              molecularWeight: molecularWeight.value,
              extinctionCoefficient: extCoefficient.value,
              baseModification: baseModification.value,
              changeLogs: modifiedLogs,
            }),
            CURRENT_USER,
          ).then(() => grok.shell.info('Posted'));
          modificationsDf.rows.addNew([
            longName.value, abbreviation.value, molecularWeight.value,
            extCoefficient.value, baseModification.value, modifiedLogs,
          ]);
        })
        .show();
    } else
      grok.shell.info('You don\'t have permission for this action');
  });
}

function deleteOverhangModification(overhangModificationsDf: DG.DataFrame, rowIndex: number): void {
  ui.dialog(
    'Do you want to delete ' + overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.getString(rowIndex) + ' ?',
  )
    .onOK(() => {
      grok.dapi.users.current().then(async (user) => {
        if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
          overhangModificationsDf.rows.removeAt(rowIndex, 1, true);
          const keyToDelete = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.get(rowIndex);
          await grok.dapi.userDataStorage.remove(STORAGE_NAME, keyToDelete, CURRENT_USER);
        } else
          grok.shell.info('You don\'t have permission for this action');
      });
    })
    .show();
}

function editOverhangModification(overhangModificationsDf: DG.DataFrame, rowIndex: number): void {
  grok.dapi.users.current().then(async (user) => {
    if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
      const longName = ui.stringInput(OVERHANG_COL_NAMES.LONG_NAMES,
        overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.get(rowIndex));
      ui.tooltip.bind(longName.root, 'Examples: \'Inverted Abasic\', \'Cyanine 3 CPG\', \'5-Methyl dC\'');
      const oldAbbreviation = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.get(rowIndex);
      const abbreviation = ui.stringInput(OVERHANG_COL_NAMES.ABBREVIATION, oldAbbreviation);
      ui.tooltip.bind(abbreviation.root, 'Examples: \'invabasic\', \'Cy3\', \'5MedC\'');
      const molecularWeight = ui.floatInput(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT,
        overhangModificationsDf.col(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT)!.get(rowIndex));
      const baseModification = ui.choiceInput(OVERHANG_COL_NAMES.BASE_MODIFICATION,
        'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
          if (v != 'NO')
            extinctionCoefficient.value = 'Base';
          extinctionCoefficient.enabled = (v == 'NO');
        });
      const extinctionCoefficient = ui.stringInput(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT,
        overhangModificationsDf.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.get(rowIndex));
      const changeLogsCol = overhangModificationsDf.col(OVERHANG_COL_NAMES.CHANGE_LOGS)!;
      ui.dialog('Add Modification')
        .add(ui.block([
          longName.root,
          abbreviation.root,
          molecularWeight.root,
          baseModification.root,
          extinctionCoefficient.root,
        ]))
        .onOK(async () => {
          if (longName.value.length > 300)
            return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');
          if (abbreviation.value.length > 100)
            return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
          const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
          if (abbreviation.value in entries)
            return grok.shell.warning('Abbreviation ' + abbreviation.value + ' already exists');
          const newLog = changeLogsCol.get(rowIndex) + Date() + ' by ' + user.firstName + ' ' + user.lastName + '; ';
          await grok.dapi.userDataStorage.postValue(
            STORAGE_NAME,
            abbreviation.value,
            JSON.stringify({
              longName: longName.value,
              abbreviation: abbreviation.value,
              molecularWeight: molecularWeight.value,
              extinctionCoefficient: extinctionCoefficient.value,
              baseModification: baseModification.value,
              changeLogs: newLog,
            }),
            CURRENT_USER,
          );
          if (oldAbbreviation != abbreviation.value)
            await grok.dapi.userDataStorage.remove(STORAGE_NAME, oldAbbreviation, CURRENT_USER);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.LONG_NAMES, rowIndex, longName.value);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.ABBREVIATION, rowIndex, abbreviation.value);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, rowIndex, molecularWeight.value);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, rowIndex, extinctionCoefficient.value);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.BASE_MODIFICATION, rowIndex, baseModification.value);
          overhangModificationsDf.set(OVERHANG_COL_NAMES.CHANGE_LOGS, rowIndex, newLog);
        })
        .show();
    } else
      grok.shell.info('You don\'t have permission for this action');
  });
}

function mergeOptions(obj1: {[index: string]: number}, obj2: {[index: string]: number}): {[index: string]: number} {
  const obj3: {[index: string]: number} = {};
  for (const attrname in obj1) {
    if (Object.prototype.hasOwnProperty.call(obj1, attrname))
      obj3[attrname] = obj1[attrname];
  }
  for (const attrname in obj2) {
    if (Object.prototype.hasOwnProperty.call(obj1, attrname))
      obj3[attrname] = obj2[attrname];
  }
  return obj3;
}

//name: Oligo Batch Calculator
//tags: app
export async function OligoBatchCalculatorApp(): Promise<void> {
  const overhangModificationsDf = await getOverhangModificationsDf();
  const overhangCodes = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.categories;
  const overhangsAbbreviations = overhangModificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.toList();
  const overhangWeights = overhangModificationsDf.col(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT)!.toList();
  const extinctionCoefficients = overhangModificationsDf.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.toList();
  const overhangsWeightsObj: {[index: string]: number} = {};
  const extinctionCoeffsObj: {[index: string]: number} = {};
  overhangsAbbreviations.forEach((key, i) => overhangsWeightsObj[key] = overhangWeights[i]);
  overhangsAbbreviations.forEach((key, i) => extinctionCoeffsObj[key] = extinctionCoefficients[i]);

  async function render(text: string): Promise<void> {
    gridDiv.innerHTML = '';

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

    for (const sequence of sequences) {
      const i = sequences.indexOf(sequence);
      const output = isValidSequence(sequence, overhangCodes);
      indicesOfFirstNotValidCharacter[i] = output.indexOfFirstNotValidCharacter;
      if (indicesOfFirstNotValidCharacter[i] < 0) {
        normalizedSequences[i] = normalizeSequence(sequence, output.expectedSynthesizer, output.expectedTechnology);
        if (normalizedSequences[i].length > 2) {
          try {
            molecularWeights[i] = molecularWeight(sequence, overhangsWeightsObj);
            extinctionCoefficients[i] = await extinctionCoefficient(normalizedSequences[i], extinctionCoeffsObj);
            nMoles[i] = await nMole(sequence, yieldAmount.value, units.value, extinctionCoeffsObj);
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
      } else if (output.expectedSynthesizer == null)
        reasonsOfError[i] = 'Not valid input';
      else {
        reasonsOfError[i] = 'Sequence is expected to be in synthesizer \'' + output.expectedSynthesizer +
          '\', please see table below to see list of valid codes';
      }
    }

    const moleName1 = (units.value == 'µmole' || units.value == 'mg') ? 'µmole' : 'nmole';
    const moleName2 = (units.value == 'µmole') ? 'µmole' : 'nmole';
    const massName = (units.value == 'µmole') ? 'mg' : (units.value == 'mg') ? units.value : 'µg';
    const c = (units.value == 'mg' || units.value == 'µmole') ? 1000 : 1;

    table = DG.DataFrame.fromColumns([
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

    const grid = DG.Viewer.grid(table, {'showRowHeader': false});
    const col = grid.col(NAME_OF_COLUMN_WITH_SEQUENCES);
    col!.cellType = 'html';
    grid.onCellPrepare(function(gc) {
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
    gridDiv.append(grid.root);
  }

  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const defaultInput = 'fAmCmGmAmCpsmU\nmApsmApsfGmAmUmCfGfAfC\nmAmUfGmGmUmCmAfAmGmA';
  let table = DG.DataFrame.create();
  const gridDiv = ui.box();

  const inputSequences = ui.textInput('', defaultInput, (txt: string) => render(txt));
  const yieldAmount = ui.floatInput('', 1, () => render(inputSequences.value));
  const units = ui.choiceInput('', 'OD', ['OD', 'µg', 'mg', 'µmole', 'nmole'], () => render(inputSequences.value));

  await render(defaultInput);

  const title = ui.panel([ui.h2('Oligo Properties')], 'ui-panel ui-box');
  title.style.maxHeight = '40px';
  $(title).children('h2').css('margin', '0px');

  const asoGapmersDf = DG.DataFrame.fromObjects([
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
  ]);
  const asoGapmersGrid = DG.Viewer.grid(asoGapmersDf!, {showRowHeader: false, showCellTooltip: false});

  const omeAndFluoroDf = DG.DataFrame.fromObjects([
    {'Name': '2\'-fluoro-U', 'BioSpring': '1', 'Axolabs': 'Uf', 'Janssen GCRS': 'fU', 'Weight': 308.16},
    {'Name': '2\'-fluoro-A', 'BioSpring': '2', 'Axolabs': 'Af', 'Janssen GCRS': 'fA', 'Weight': 331.2},
    {'Name': '2\'-fluoro-C', 'BioSpring': '3', 'Axolabs': 'Cf', 'Janssen GCRS': 'fC', 'Weight': 307.18},
    {'Name': '2\'-fluoro-G', 'BioSpring': '4', 'Axolabs': 'Gf', 'Janssen GCRS': 'fG', 'Weight': 347.19},
    {'Name': '2\'OMe-rU', 'BioSpring': '5', 'Axolabs': 'u', 'Janssen GCRS': 'mU', 'Weight': 320.2},
    {'Name': '2\'OMe-rA', 'BioSpring': '6', 'Axolabs': 'a', 'Janssen GCRS': 'mA', 'Weight': 343.24},
    {'Name': '2\'OMe-rC', 'BioSpring': '7', 'Axolabs': 'c', 'Janssen GCRS': 'mC', 'Weight': 319.21},
    {'Name': '2\'OMe-rG', 'BioSpring': '8', 'Axolabs': 'g', 'Janssen GCRS': 'mG', 'Weight': 359.24},
    {'Name': 'ps linkage', 'BioSpring': '*', 'Axolabs': 's', 'Janssen GCRS': 'ps', 'Weight': 16.07},
  ]);
  const omeAndFluoroGrid = DG.Viewer.grid(omeAndFluoroDf!, {showRowHeader: false, showCellTooltip: false});

  const overhangModifsGrid = DG.Viewer.grid(overhangModificationsDf, {showRowHeader: false, showCellTooltip: false});
  overhangModifsGrid.col(OVERHANG_COL_NAMES.LONG_NAMES)!.width = 110;
  overhangModifsGrid.col(OVERHANG_COL_NAMES.ABBREVIATION)!.width = 80;
  overhangModifsGrid.col(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT)!.width = 105;
  overhangModifsGrid.col(OVERHANG_COL_NAMES.BASE_MODIFICATION)!.width = 110;
  overhangModifsGrid.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.width = 100;

  const codesTablesDiv = ui.splitV([
    ui.box(ui.h2('ASO Gapmers'), {style: {maxHeight: '40px'}}),
    asoGapmersGrid.root,
    ui.box(ui.h2('2\'-OMe and 2\'-F modifications'), {style: {maxHeight: '40px'}}),
    omeAndFluoroGrid.root,
    ui.box(ui.h2('Additional modifications'), {style: {maxHeight: '40px'}}),
    overhangModifsGrid.root,
  ], {style: {maxWidth: '600px'}});

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
            ], 'inputSequence'),
          ]), {style: {maxHeight: '230px'}},
        ),
        ui.splitV([
          title,
          gridDiv,
        ]),
      ]),
      codesTablesDiv,
    ]),
  ]);
  view.box = true;

  const switchInput = ui.switchInput('Codes', true, (v: boolean) => (v) ?
    $(codesTablesDiv).show() :
    $(codesTablesDiv).hide(),
  );

  const col = overhangModifsGrid.col(OVERHANG_COL_NAMES.ACTION)!;
  col.cellType = 'html';
  overhangModifsGrid.onCellPrepare(function(gc) {
    if (gc.isTableCell && gc.gridColumn.name == OVERHANG_COL_NAMES.ACTION) {
      gc.style.element = ui.divH([
        ui.button(ui.iconFA('trash-alt'), () => deleteOverhangModification(overhangModificationsDf, gc.gridRow)),
        ui.button(ui.iconFA('edit'), () => editOverhangModification(overhangModificationsDf, gc.gridRow)),
      ]);
    }
  });

  view.setRibbonPanels([[
    ui.iconFA('redo', () => inputSequences.value = ''),
    ui.iconFA('plus', () => addModificationButton(overhangModificationsDf)),
    ui.iconFA('arrow-to-bottom', () => saveAsCsv(table)),
    switchInput.root,
  ]]);

  $('.inputSequence textarea')
    .css('resize', 'none')
    .css('min-height', '70px')
    .css('width', '100%')
    .css('font-family', 'monospace')
    .attr('spellcheck', 'false');
}
