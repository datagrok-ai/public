import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {WEIGHTS, INDIVIDUAL_BASES, NEAREST_NEIGHBOUR, SYNTHESIZERS, CURRENT_USER, STORAGE_NAME, DEFAULT_INPUT,
  ADDITIONAL_MODS_COL_NAMES, MAIN_COL_NAMES, UNITS, EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION} from './constants';
import {isValidSequence, validate} from './validation';
import {deleteWord, saveAsCsv, sortByStringLengthInDescOrder, mergeOptions, normalizeSequence,
  isCurrentUserAppAdmin} from './helpers';
import {addModification, editModification} from './additional-modifications';
import {opticalDensityCalc, molecularMassCalc, nMoleCalc} from './calculations-simplified';

export const _package = new DG.Package();
const windows = grok.shell.windows;
windows.showProperties = false;
windows.showToolbox = false;
windows.showHelp = false;

//name: getUnits
//output: list<string> units
export function getUnits(): string[] {
  return Object.values(UNITS);
}

//name: opticalDensity
//input: string sequence
//input: double amount
//input: string outputUnits {choices: OligoBatchCalculator: getUnits}
//output: double opticalDensity
export async function opticalDensity(sequence: string, amount: number, outputUnits: string,
  extCoefsObj: {[index: string]: number}): Promise<number> {
  const ec = await extinctionCoefficient(sequence, extCoefsObj);
  if (outputUnits == UNITS.MILLI_GRAM && outputUnits == UNITS.MICRO_GRAM)
    return (outputUnits == UNITS.MICRO_GRAM ? 1 : 0.001) * amount * ec / molecularWeight(sequence);
  if (outputUnits == UNITS.OPTICAL_DENSITY)
    return amount;
  const coefficient = (outputUnits == UNITS.NANO_MOLE) ? 1000000 : (outputUnits == UNITS.MILLI_GRAM) ? 1 : 1000;
  return amount * ec / coefficient;
}

//name: nMole
//input: string sequence
//input: double amount
//input: string outputUnits {choices: OligoBatchCalculator: getUnits}
//output: double nMole
export async function nMole(sequence: string, amount: number, outputUnits: string, extinctionCoefficientsObj:
  {[index: string]: number}, weightsObj: {[index: string]: number}): Promise<number> {
  const ec = await extinctionCoefficient(sequence, extinctionCoefficientsObj);
  return (outputUnits == UNITS.OPTICAL_DENSITY) ?
    1000000 * amount / ec :
    1000 * amount / molecularWeight(sequence, weightsObj);
}

//name: molecularMass
//input: string sequence
//input: double amount
//input: string outputUnits {choices: OligoBatchCalculator: getUnits}
//output: double molecularMass
export async function molecularMass(sequence: string, amount: number, outputUnits: string): Promise<number> {
  const additionalWeightsObj: {[index: string]: number} = {};
  const extinctionCoeffsObj: {[index: string]: number} = {};
  const modifications: any[] = [];
  const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
  const invalidKeys = [
    'baseModification', 'extinctionCoefficient', 'molecularWeight', 'abbreviation', 'longName', 'changeLogs',
  ];
  for (const key of Object.keys(entries)) {
    if (!invalidKeys.includes(key))
      modifications.push(JSON.parse(entries[key]));
  }
  const additionalAbbreviations = modifications.map((e) => e.abbreviation);
  const additionalWeights = modifications.map((e) => (e.molecularWeight == undefined) ? 0 : e.molecularWeight);
  const extinctionCoefficients = modifications.map((e) => e.extinctionCoefficient);
  additionalAbbreviations.forEach((key, i) => {
    additionalWeightsObj[key] = additionalWeights[i];
    if (extinctionCoefficients[i] != EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION)
      extinctionCoeffsObj[key] = extinctionCoefficients[i] ?? 1;
  });
  const ec = await extinctionCoefficient(sequence, extinctionCoeffsObj);
  const od = await opticalDensity(sequence, amount, outputUnits, extinctionCoeffsObj);
  const nm = await nMole(sequence, amount, outputUnits, extinctionCoeffsObj, additionalWeightsObj);
  if (outputUnits == UNITS.OPTICAL_DENSITY) {
    return (ec == 0) ?
      amount * molecularWeight(sequence, additionalWeightsObj) :
      1000 * amount * molecularWeight(sequence, additionalWeightsObj) / ec;
  }
  const coefficient = (outputUnits == UNITS.MILLI_GRAM) ? 1 : 1000;
  return amount / ec * molecularWeight(sequence) * coefficient * od / nm;
}

//name: molecularWeight
//input: string sequence
//input: string additionalWeightsObj
//output: double molWeight
export function molecularWeight(sequence: string, additionalWeightsObj?: {[index: string]: number}): number {
  const codes = (additionalWeightsObj == null) ?
    sortByStringLengthInDescOrder(Object.keys(WEIGHTS)) :
    sortByStringLengthInDescOrder(Object.keys(WEIGHTS).concat(Object.keys(additionalWeightsObj)));
  const obj = (additionalWeightsObj != null) ? mergeOptions(WEIGHTS, additionalWeightsObj) : WEIGHTS;
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += obj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode.length;
  }
  return weight - 61.97;
}

export async function extinctionCoefficient(sequence: string, extCoefsObj?: {[i: string]: number}): Promise<number> {
  const modifications: any[] = [];
  const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
  const invalidKeys = [
    'baseModification', 'extinctionCoefficient', 'molecularWeight', 'abbreviation', 'longName', 'changeLogs',
  ];
  for (const key of Object.keys(entries)) {
    if (!invalidKeys.includes(key))
      modifications.push(JSON.parse(entries[key]));
  }
  const molWeightList = modifications.map((e) => (e.molecularWeight == undefined) ? '0' : String(e.molecularWeight));
  const extinctionCoefList = modifications.map((e) => String(e.extinctionCoefficient));
  const additionalModsDf = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES, modifications.map((e) => e.longName)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION, modifications.map((e) => e.abbreviation)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT, molWeightList),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION, modifications.map((e) => e.baseModification)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT, extinctionCoefList),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.ACTION, Array(modifications.length)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS, modifications.map((e) => e.changeLogs)),
  ]);
  let ns = normalizeSequence(sequence, SYNTHESIZERS.GCRS, null, additionalModsDf);
  let nearestNeighbourSum = 0;
  let individualBasisSum = 0;
  let modificationsSum = 0;
  if (extCoefsObj != null) {
    for (const modif of Object.keys(extCoefsObj)) {
      if (
        String(extCoefsObj[modif]) != EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION &&
        extCoefsObj[modif] != undefined &&
        !isNaN(Number(extCoefsObj[modif]))
      ) {
        modificationsSum += (sequence.match(new RegExp(modif, 'g')) || []).length * Number(extCoefsObj[modif]);
        ns = deleteWord(ns, modif);
      }
    }
  }
  for (let i = 0; i < ns.length - 2; i += 2) {
    nearestNeighbourSum += (ns[i] == ns[i + 2]) ?
      NEAREST_NEIGHBOUR[ns.slice(i, i + 2)][ns.slice(i + 2, i + 4)] :
      (
        NEAREST_NEIGHBOUR['r' + ((ns[i + 1] == 'T') ? 'U' : ns[i + 1])]['r' + ((ns[i + 3] == 'T') ? 'U' : ns[i + 3])] +
        NEAREST_NEIGHBOUR['d' + ((ns[i + 1] == 'U') ? 'T' : ns[i + 1])]['d' + ((ns[i + 3] == 'U') ? 'T' : ns[i + 3])]
      ) / 2;
  }
  for (let i = 2; i < ns.length - 2; i += 2)
    individualBasisSum += INDIVIDUAL_BASES[ns.slice(i, i + 2)];
  return nearestNeighbourSum - individualBasisSum + modificationsSum;
}

//name: Oligo Batch Calculator
//tags: app
export async function OligoBatchCalculatorApp(): Promise<void> {
  const additionalWeightsObj: {[index: string]: number} = {};
  const extinctionCoeffsObj: {[index: string]: number} = {};
  const modifications: any[] = [];
  const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
  const invalidKeys = [
    'baseModification', 'extinctionCoefficient', 'molecularWeight', 'abbreviation', 'longName', 'changeLogs',
  ];
  for (const key of Object.keys(entries)) {
    if (!invalidKeys.includes(key))
      modifications.push(JSON.parse(entries[key]));
  }
  const molWeightList = modifications.map((e) => (e.molecularWeight == undefined) ? '0' : String(e.molecularWeight));
  const extinctionCoefList = modifications.map((e) => String(e.extinctionCoefficient));
  const additionalModsDf = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES, modifications.map((e) => e.longName)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION, modifications.map((e) => e.abbreviation)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT, molWeightList),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION, modifications.map((e) => e.baseModification)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT, extinctionCoefList),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.ACTION, Array(modifications.length)),
    DG.Column.fromStrings(ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS, modifications.map((e) => e.changeLogs)),
  ]);
  const additionalAbbreviations = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).toList();
  const additionalWeights = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT).toList();
  const extinctionCoefficients = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT).toList();
  additionalAbbreviations.forEach((key, i) => {
    additionalWeightsObj[key] = additionalWeights[i];
    if (extinctionCoefficients[i] != EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION)
      extinctionCoeffsObj[key] = extinctionCoefficients[i] ?? 1;
  });

  const mainGrid = DG.Viewer.grid(DG.DataFrame.create(), {
    showRowHeader: false,
    allowEdit: false,
    showCellTooltip: true,
  });

  async function render(text: string): Promise<void> {
    const additionalAbbreviations = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).toList();
    const additionalWeights = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT).toList();
    const extinctionCoefficients2 = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT).toList();
    additionalAbbreviations.forEach((key, i) => {
      additionalWeightsObj[key] = additionalWeights[i];
      if (extinctionCoefficients2[i] != EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION)
        extinctionCoeffsObj[key] = extinctionCoefficients2[i] ?? 1;
    });

    const sequences = text.split('\n')
      .map((s) => s.replace(/\s/g, ''))
      .filter((item) => item);

    const indicesOfFirstNotValidCharacter = Array(sequences.length);
    const normalizedSequences = Array(sequences.length);
    const molecularWeights = new Float32Array(sequences.length);
    const extinctionCoefficients = new Float32Array(sequences.length);
    const nMoles = new Float32Array(sequences.length);
    const opticalDensities = new Float32Array(sequences.length);
    const molecularMasses = new Float32Array(sequences.length);
    const reasonsOfError = Array(sequences.length);

    for (const [i, sequence] of sequences.entries()) {
      indicesOfFirstNotValidCharacter[i] = validate(sequence, additionalAbbreviations);
      if (isValidSequence(indicesOfFirstNotValidCharacter[i])) {
        normalizedSequences[i] = normalizeSequence(sequence, SYNTHESIZERS.GCRS, null, additionalModsDf);
        if (normalizedSequences[i].length > 2) {
          try {
            molecularWeights[i] = molecularWeight(sequence, additionalWeightsObj);
            extinctionCoefficients[i] = await extinctionCoefficient(normalizedSequences[i], extinctionCoeffsObj);
            nMoles[i] = nMoleCalc(yieldAmount.value!, units.value!, molecularWeights[i], extinctionCoefficients[i]);
            opticalDensities[i] = opticalDensityCalc(sequence, yieldAmount.value!, units.value!,
              extinctionCoefficients[i]);
            molecularMasses[i] = molecularMassCalc(sequence, yieldAmount.value!, units.value!,
              extinctionCoefficients[i], opticalDensities[i], nMoles[i], additionalWeightsObj);
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
        normalizedSequences[i] = [];
        molecularWeights[i] = NaN;
        extinctionCoefficients[i] = NaN;
        nMoles[i] = NaN;
        opticalDensities[i] = NaN;
        molecularMasses[i] = NaN;
        reasonsOfError[i] = 'Sequence is expected to be in synthesizer \'' + SYNTHESIZERS.GCRS +
          '\', please see table below to see list of valid codes';
      }
    };

    const moleColumnName = (units.value == UNITS.MICRO_MOLE || units.value == UNITS.MILLI_GRAM) ?
      UNITS.MICRO_MOLE : UNITS.NANO_MOLE;
    const moleName2 = (units.value == UNITS.MICRO_MOLE) ? UNITS.MICRO_MOLE : UNITS.NANO_MOLE;
    const massName = (units.value == UNITS.MICRO_MOLE) ?
      UNITS.MILLI_GRAM :
      (units.value == UNITS.MILLI_GRAM) ?
        units.value :
        UNITS.MICRO_GRAM;
    const c = (units.value == UNITS.MILLI_GRAM || units.value == UNITS.MICRO_MOLE) ? 1000 : 1;

    mainGrid.dataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, MAIN_COL_NAMES.ITEM,
        Array(...Array(sequences.length + 1).keys()).slice(1)),
      DG.Column.fromStrings(MAIN_COL_NAMES.SEQUENCE, sequences),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, MAIN_COL_NAMES.LENGTH, normalizedSequences.map((s) => s.length / 2)),
      DG.Column.fromFloat32Array(MAIN_COL_NAMES.OPTICAL_DENSITY, opticalDensities),
      DG.Column.fromFloat32Array(moleColumnName, nMoles),
      DG.Column.fromFloat32Array(`Mass [${massName}]`, molecularMasses),
      DG.Column.fromFloat32Array(`${moleName2}/OD`, nMoles.map(function(n, i) {return c * n / opticalDensities[i];})),
      DG.Column.fromFloat32Array(MAIN_COL_NAMES.MASS_OD_RATIO,
        molecularMasses.map(function(n, i) {return c * n / opticalDensities[i];})),
      DG.Column.fromFloat32Array(MAIN_COL_NAMES.MOLECULAR_WEIGHT, molecularWeights),
      DG.Column.fromFloat32Array(MAIN_COL_NAMES.EXTINCTION_COEFFICIENT, extinctionCoefficients),
    ]);

    const col = mainGrid.col(MAIN_COL_NAMES.SEQUENCE)!;
    col.cellType = 'html';
    mainGrid.onCellPrepare(function(gc) {
      if (gc.isTableCell && gc.gridColumn.name == MAIN_COL_NAMES.SEQUENCE) {
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

  const inputSequences = ui.textInput('', DEFAULT_INPUT, (txt: string) => render(txt));
  const yieldAmount = ui.floatInput('', 1, () => render(inputSequences.value));
  const units = ui.choiceInput('', UNITS.OPTICAL_DENSITY, Object.values(UNITS), () => render(inputSequences.value));

  await render(DEFAULT_INPUT);

  const downloadIcon = ui.iconFA('download', () => saveAsCsv(mainGrid.dataFrame), 'Save as CSV file');
  $(downloadIcon).css('margin-left', '5px');

  const title = ui.panel([
    ui.divH([
      ui.h2('Oligo Properties'),
      downloadIcon,
    ], {style: {'display': 'flex', 'align-items': 'center'}}),
  ], 'ui-panel ui-box');
  title.style.maxHeight = '40px';
  $(title).children('h2').css('margin', '0px');

  const additionaModifsGrid = DG.Viewer.grid(additionalModsDf, {
    showRowHeader: false,
    showCellTooltip: true,
    allowEdit: (await isCurrentUserAppAdmin()),
  });

  const addModificationIcon = ui.iconFA('plus', () => addModification(additionalModsDf), 'Add new modidfication');
  $(addModificationIcon).css('margin-left', '5px');
  $(addModificationIcon).css('margin-top', '12px');

  const codesTablesDiv = ui.splitV([
    ui.box(
      ui.divH([
        ui.h2('Additional modifications'),
        addModificationIcon,
      ]), {style: {maxHeight: '40px'}},
    ),
    additionaModifsGrid.root,
  ], {style: {maxWidth: '600px'}});

  const clearIcon = ui.iconFA('redo', () => inputSequences.value = '', 'Clear input field');
  $(clearIcon).css('margin-left', '5px');
  $(clearIcon).css('margin-top', '12px');

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
            ui.divH([ui.h2('Input Sequences'), clearIcon]),
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
    ui.switchInput('Show additional modifications', true, (v: boolean) => {
      (v) ? $(codesTablesDiv).show() : $(codesTablesDiv).hide();
    }).root,
  ]]);

  editModification(additionalModsDf, additionaModifsGrid);
}
