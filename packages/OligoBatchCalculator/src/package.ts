/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {map} from "./map";

export let _package = new DG.Package();

const STORAGE_NAME = 'oligo-batch-calculator-storage';
const OVERHANG_COL_NAMES = {
  LONG_NAMES: 'Long name',
  ABBREVIATION: 'Abbreviation',
  MOLECULAR_WEIGHT: 'Molecular weight',
  BASE_MODIFICATION: 'Base modification',
  EXTINCTION_COEFFICIENT: 'Ext. coefficient',
  ACTION: 'Action'
}

let weightsObj: {[code: string]: number} = {};
let normalizedObj: {[code: string]: string} = {};
for (let synthesizer of Object.keys(map))
  for (let technology of Object.keys(map[synthesizer]))
    for (let code of Object.keys(map[synthesizer][technology])) {
      weightsObj[code] = map[synthesizer][technology][code].weight;
      normalizedObj[code] = map[synthesizer][technology][code].normalized;
    }

function saveAsCsv(table: DG.DataFrame) {
  let link = document.createElement("a");
  link.setAttribute("href", "data:text/csv;charset=utf-8,\uFEFF" + encodeURI(table.toCsv()));
  link.setAttribute("download", "Oligo Properties.csv");
  link.click();
}

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) { return b.length - a.length; });
}

async function getModificationsDf() {
  // //@ts-ignore
  // await grok.dapi.userDataStorage.remove(STORAGE_NAME, null, false);
  let modifications: any[] = [];
  let entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, false);
  if (entries !== null && Object.keys(entries).length === 0)
    grok.shell.info('Storage is empty. Try to post something to the storage');
  else {
    Object.keys(entries).forEach((key) => {
      modifications.push(JSON.parse(entries[key]));
    });
  }
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(OVERHANG_COL_NAMES.LONG_NAMES, modifications.map((e) => e.longName)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.ABBREVIATION, modifications.map((e) => e.abbreviation)),    // @ts-ignore
    DG.Column.fromFloat32Array(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, modifications.map((e) => (e.molecularWeight == undefined) ? 0 : e.molecularWeight)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.BASE_MODIFICATION, modifications.map((e) => e.baseModification)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, modifications.map((e) => e.extinctionCoefficient)),
    DG.Column.fromStrings(OVERHANG_COL_NAMES.ACTION, Array(modifications.length).fill(''))
  ])!;
}

//name: Oligo Batch Calculator
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['NMole', 'Milligrams', 'Micrograms', 'Optical Density']}
//output: double opticalDensity
//output: double nMole
//output: double molecularMass
//output: double molecularWeight
//output: double extinctionCoefficient
export function OligoBatchCalculator(sequence: string, amount: number, outputUnits: string) {
  return {
    opticalDensity: opticalDensity(sequence, amount, outputUnits),
    nMole: nMole(sequence, amount, outputUnits),
    molecularMass: molecularMass(sequence, amount, outputUnits),
    molecularWeight: molecularWeight(sequence),
    extinctionCoefficient: extinctionCoefficient(sequence)
  };
}

//name: opticalDensity
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['NMole', 'Milligrams', 'Micrograms']}
//output: double opticalDensity
export function opticalDensity(sequence: string, amount: number, outputUnits: string): number {
  if (outputUnits == 'Milligrams' || outputUnits == 'Micrograms' || outputUnits == 'mg' || outputUnits == 'µg')
    return (outputUnits == 'Milligrams' ? 1 : 0.001) * amount * extinctionCoefficient(sequence) / molecularWeight(sequence);
  if (outputUnits == 'OD')
    return amount;
  const coefficient = (outputUnits == 'NMole') ? 1000000 : (outputUnits == 'Milligrams') ? 1 : 1000;
  return amount * extinctionCoefficient(sequence) / coefficient;
}

//name: nMole
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micrograms']}
//output: double nMole
export function nMole(sequence: string, amount: number, outputUnits: string): number {
  return (outputUnits == 'Optical Density') ? 1000000 * amount / extinctionCoefficient(sequence) : 1000 * amount / molecularWeight(sequence);
}

//name: molecularMass
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micromoles', 'Millimoles']}
//output: double molecularMass
export function molecularMass(sequence: string, amount: number, outputUnits: string): number {
  if (outputUnits == 'Optical Density' || outputUnits == 'OD') {
    let ec = extinctionCoefficient(sequence);
    return (ec == 0) ? amount * molecularWeight(sequence) : 1000 * amount * molecularWeight(sequence) / ec;
  }
  const coefficient = (outputUnits == 'Milligrams' || outputUnits == 'Micromoles') ? 1 : 1000;
  return amount / extinctionCoefficient(sequence) * molecularWeight(sequence) * coefficient * opticalDensity(sequence, amount, outputUnits) / nMole(sequence, amount, outputUnits);
}

//name: molecularWeight
//input: string sequence
//output: double molecularWeight
export function molecularWeight(sequence: string): number {
  const codes = sortByStringLengthInDescendingOrder(Object.keys(weightsObj));
  let weight = 0, i = 0;
  while (i < sequence.length) {
    let matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += weightsObj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode!.length;
  }
  return weight - 61.97;
}

//name: extinctionCoefficient
//input: string sequence
//output: double extinctionCoefficient
export function extinctionCoefficient(sequence: string): number {
  let output = isValid(sequence);
  sequence = normalizeSequence(sequence, output.expectedSynthesizer!);
  const individualBases: {[index: string]: number} = {
    'dA': 15400, 'dC': 7400, 'dG': 11500, 'dT': 8700, 'rA': 15400, 'rC': 7200, 'rG': 11500, 'rU': 9900
  },
  nearestNeighbour: {[firstBase: string]: {[secondBase: string]: number}} = {
    'dA': {'dA': 27400, 'dC': 21200, 'dG': 25000, 'dT': 22800, 'rA': 27400, 'rC': 21000, 'rG': 25000, 'rU': 24000},
    'dC': {'dA': 21200, 'dC': 14600, 'dG': 18000, 'dT': 15200, 'rA': 21000, 'rC': 14200, 'rG': 17800, 'rU': 16200},
    'dG': {'dA': 25200, 'dC': 17600, 'dG': 21600, 'dT': 20000, 'rA': 25200, 'rC': 17400, 'rG': 21600, 'rU': 21200},
    'dT': {'dA': 23400, 'dC': 16200, 'dG': 19000, 'dT': 16800, 'rA': 24600, 'rC': 17200, 'rG': 20000, 'rU': 19600},
    'rA': {'rA': 27400, 'rC': 21000, 'rG': 25000, 'rU': 24000, 'dA': 27400, 'dC': 21200, 'dG': 25000, 'dT': 22800},
    'rC': {'rA': 21000, 'rC': 14200, 'rG': 17800, 'rU': 16200, 'dA': 21200, 'dC': 14600, 'dG': 18000, 'dT': 15200},
    'rG': {'rA': 25200, 'rC': 17400, 'rG': 21600, 'rU': 21200, 'dA': 25200, 'dC': 17600, 'dG': 21600, 'dT': 20000},
    'rU': {'rA': 24600, 'rC': 17200, 'rG': 20000, 'rU': 19600, 'dA': 23400, 'dC': 16200, 'dG': 19000, 'dT': 16800}
  };
  let ec1 = 0, ec2 = 0;
  for (let i = 0; i < sequence.length - 2; i += 2)
    ec1 += (sequence[i] == sequence[i + 2]) ?
      nearestNeighbour[sequence.slice(i, i + 2)][sequence.slice(i + 2, i + 4)] :
      (
        nearestNeighbour['r' + ((sequence[i + 1] == 'T') ? 'U' : sequence[i + 1])]['r' + ((sequence[i + 3] == 'T') ? 'U' : sequence[i + 3])]
        +
        nearestNeighbour['d' + ((sequence[i + 1] == 'U') ? 'T' : sequence[i + 1])]['d' + ((sequence[i + 3] == 'U') ? 'T' : sequence[i + 3])]
      ) / 2;
  for (let i = 2; i < sequence.length - 2; i += 2)
    ec2 += individualBases[sequence.slice(i, i + 2)];
  return ec1 - ec2;
}

function normalizeSequence(sequence: string, synthesizer: string): string {
  const codes = sortByStringLengthInDescendingOrder(getAllCodesOfSynthesizer(synthesizer));
  const re = new RegExp('(' + codes.join('|') + ')', 'g');
  return sequence.replace(re, function (code) {return normalizedObj[code]});
}

function getAllCodesOfSynthesizer(synthesizer: string) {
  let codes: string[] = [];
  for (let technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes;
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string): string[] {
  let synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    if (codes.some((s) => s == sequence.slice(0, s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

function isValid(sequence: string) {
  let possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence);
  if (possibleSynthesizers.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedRepresentation: null };

  let outputIndices = Array(possibleSynthesizers.length).fill(0);

  const firstUniqueCharacters = ['r', 'd'], nucleotides = ["A", "U", "T", "C", "G"];

  possibleSynthesizers.forEach((synthesizer, synthesizerIndex) => {
    let codes = getAllCodesOfSynthesizer(synthesizer);
    while (outputIndices[synthesizerIndex] < sequence.length) {

      let matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[synthesizerIndex], outputIndices[synthesizerIndex] + c.length));

      if (matchedCode == null)
        break;

      if (  // for mistake pattern 'rAA'
        outputIndices[synthesizerIndex] > 1 &&
        nucleotides.includes(sequence[outputIndices[synthesizerIndex]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[synthesizerIndex] - 2])
      ) break;

      if (  // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[synthesizerIndex] + 1]) &&
        nucleotides.includes(sequence[outputIndices[synthesizerIndex]])
      ) {
        outputIndices[synthesizerIndex]++;
        break;
      }

      outputIndices[synthesizerIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedSythesizer = Math.max.apply(Math, outputIndices);
  const indexOfFirstNotValidCharacter = (indexOfExpectedSythesizer == sequence.length) ? -1 : indexOfExpectedSythesizer;

  return {
    indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
    expectedSynthesizer: possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)]
  };
}

//name: Oligo Batch Calculator
//tags: app
export async function OligoBatchCalculatorApp() {

  function updateTable(text: string) {
    gridDiv.innerHTML = '';

    let sequences = text.split('\n')
      .map((s) => s.replace(/\s/g, ''))
      .filter(item => item);

    let indicesOfFirstNotValidCharacter = Array(sequences.length),
      normalizedSequences = Array(sequences.length),
      molecularWeights = Array(sequences.length),
      extinctionCoefficients = Array(sequences.length),
      nMoles = Array(sequences.length),
      opticalDensities = Array(sequences.length),
      molecularMasses = Array(sequences.length),
      reasonsOfError = Array(sequences.length),
      expectedSynthesizers = Array(sequences.length);

    sequences.forEach((sequence, i) => {
      let output = isValid(sequence);
      indicesOfFirstNotValidCharacter[i] = output.indexOfFirstNotValidCharacter;
      expectedSynthesizers[i] = output.expectedSynthesizer;
      if (indicesOfFirstNotValidCharacter[i] < 0) {
        normalizedSequences[i] = normalizeSequence(sequence, expectedSynthesizers[i]);
        if (normalizedSequences[i].length > 2) {
          try {
            molecularWeights[i] = molecularWeight(sequence);
            extinctionCoefficients[i] = extinctionCoefficient(normalizedSequences[i]);
            nMoles[i] = nMole(sequence, yieldAmount.value, units.value);
            opticalDensities[i] = opticalDensity(sequence, yieldAmount.value, units.value);
            molecularMasses[i] = molecularMass(sequence, yieldAmount.value, units.value);
          } catch (e) {
            reasonsOfError[i] = 'Unknown error, please report it to Datagrok team';
            indicesOfFirstNotValidCharacter[i] = 0;
            grok.shell.error(String(e));
          }
        } else {
          reasonsOfError[i] = 'Sequence should contain at least two nucleotides';
          indicesOfFirstNotValidCharacter[i] = 0;
        }
      } else if (expectedSynthesizers[i] == null)
        reasonsOfError[i] = "Not valid input";
      else
        reasonsOfError[i] = "Sequence is expected to be in synthesizer format '" +  expectedSynthesizers[i] +
          "', please see table below to see list of valid codes";
    });

    const moleName1 = (units.value == 'µmole' || units.value == 'mg') ? 'µmole' : 'nmole',
      moleName2 = (units.value == 'µmole') ? 'µmole' : 'nmole',
      massName = (units.value == 'µmole') ? 'mg' : (units.value == 'mg') ? units.value : 'µg',
      coefficient = (units.value == 'mg' || units.value == 'µmole') ? 1000 : 1,
      nameOfColumnWithSequences = 'Sequence';

    table = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'Item', Array(...Array(sequences.length + 1).keys()).slice(1)),
      DG.Column.fromList('string', nameOfColumnWithSequences, sequences),
      DG.Column.fromList('int', 'Length', normalizedSequences.map((s) => s.length / 2)),
      DG.Column.fromList('double', 'OD 260', opticalDensities),
      DG.Column.fromList('double', moleName1, nMoles),
      DG.Column.fromList('double', 'Mass (' + massName + ')', molecularMasses),
      DG.Column.fromList('double', moleName2 + '/OD', nMoles.map(function(n, i) {return coefficient * n / opticalDensities[i]})),
      DG.Column.fromList('double', 'µg/OD', molecularMasses.map(function(n, i) {return coefficient * n / opticalDensities[i]})),
      DG.Column.fromList('double', 'MW', molecularWeights),
      DG.Column.fromList('int', 'Ext. Coefficient', extinctionCoefficients)
    ]);

    let grid = DG.Viewer.grid(table, { 'showRowHeader': false });
    let col = grid.columns.byName(nameOfColumnWithSequences);
    col!.cellType = 'html';

    grid.onCellPrepare(function (gc) {
      if (gc.isTableCell && gc.gridColumn.name == nameOfColumnWithSequences) {
        let items = (indicesOfFirstNotValidCharacter[gc.gridRow] < 0) ?
          [ui.divText(gc.cell.value, {style: {color: "grey"}})] :
          [
            ui.divText(gc.cell.value.slice(0, indicesOfFirstNotValidCharacter[gc.gridRow]), {style: {color: "grey"}}),
            ui.tooltip.bind(
              ui.divText(gc.cell.value.slice(indicesOfFirstNotValidCharacter[gc.gridRow]), {style: {color: "red"}}),
              reasonsOfError[gc.gridRow]
            )
          ];
        gc.style.element = ui.divH(items, {style: {margin: '6px 0 0 6px'}});
      }
    });

    gridDiv.append(grid.root);
  }

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const defaultInput = 'fAmCmGmAmCpsmU\nmApsmApsfGmAmUmCfGfAfC\nmAmUfGmGmUmCmAfAmGmA';
  let table = DG.DataFrame.create();
  let gridDiv = ui.box();

  let inputSequences = ui.textInput("", defaultInput, async (txt: string) => updateTable(txt));
  let yieldAmount = ui.floatInput('', 1, () => updateTable(inputSequences.value));
  let units = ui.choiceInput('', 'OD', ['OD', 'µg', 'mg', 'µmole', 'nmole'], () => updateTable(inputSequences.value));

  updateTable(defaultInput);

  let tables = ui.divV([]);
  for (let synthesizer of Object.keys(map)) {
    for (let technology of Object.keys(map[synthesizer])) {
      let tableRows = [];
      for (let [key, value] of Object.entries(map[synthesizer][technology]))
        tableRows.push({'name': value.name, 'code': key, 'weight': value['weight']});
      tables.append(
        DG.HtmlTable.create(
          tableRows,
          (item: {name: string; code: string; weight: number}) => [item['name'], item['code'], item['weight']],
          [synthesizer + ' ' + technology, 'Code', 'Weight']
        ).root,
        ui.div([], {style: {height: '30px'}})
      );
    }
  }

  async function addModificationButton() {
    grok.dapi.users.current().then((user) => {
      if (adminUsers.includes(user.firstName + ' ' + user.lastName)) {
        let longName = ui.stringInput(OVERHANG_COL_NAMES.LONG_NAMES, '');
        ui.tooltip.bind(longName.root, "Examples: 'Inverted Abasic', 'Cyanine 3 CPG', '5-Methyl dC'");
        let abbreviation = ui.stringInput(OVERHANG_COL_NAMES.ABBREVIATION, '');
        ui.tooltip.bind(abbreviation.root, "Examples: 'invabasic', 'Cy3', '5MedC'");
        let molecularWeight = ui.floatInput(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, 0);
        let baseModification = ui.choiceInput(OVERHANG_COL_NAMES.BASE_MODIFICATION, 'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
          if (v != 'NO')
            extinctionCoefficient.value = 'Base';
          extinctionCoefficient.enabled = (v == 'NO');
        });
        let extinctionCoefficient = ui.stringInput(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, '');
        ui.dialog('Add Modification')
          .add(ui.block([
            longName.root,
            abbreviation.root,
            molecularWeight.root,
            baseModification.root,
            extinctionCoefficient.root,
          ]))
          .onOK(() => {
            if (longName.value.length > 300)
              return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters')
            if (abbreviation.value.length > 100)
              return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters')
            grok.dapi.userDataStorage.postValue(
              STORAGE_NAME,
              longName.value,
              JSON.stringify({
                longName: longName.value,
                abbreviation: abbreviation.value,
                molecularWeight: molecularWeight.value,
                extinctionCoefficient: extinctionCoefficient.value,
                baseModification: baseModification.value,
              }),
              false
            ).then(() => grok.shell.info('Posted'));
            modificationsDf.rows.addNew([
              longName.value, abbreviation.value, molecularWeight.value, extinctionCoefficient.value, baseModification.value
            ]);
          })
          .show();
      } else
        grok.shell.info('You don\'t have permission for this action');
    });
  }

  let title = ui.panel([ui.h2('Oligo Properties')], 'ui-panel ui-box');
  title.style.maxHeight = '40px';
  $(title).children('h2').css('margin', '0px');

  const asoGapmersDf = DG.DataFrame.fromObjects([
    { Name: "2'MOE-5Me-rU", BioSpring: "5", "Janssen GCRS": "moeT", Weight: 378.27 },
    { Name: "2'MOE-rA", BioSpring: "6", "Janssen GCRS": "moeA", Weight: 387.29 },
    { Name: "2'MOE-5Me-rC", BioSpring: "7", "Janssen GCRS": "moe5mC", Weight: 377.29 },
    { Name: "2'MOE-rG", BioSpring: "8", "Janssen GCRS": "moeG", Weight: 403.28 },
    { Name: "5-Methyl-dC", BioSpring: "9", "Janssen GCRS": "5mC", Weight: 303.21 },
    { Name: "ps linkage", BioSpring: "*", "Janssen GCRS": "ps", Weight: 16.07 },
    { Name: "dA", BioSpring: "A", "Janssen GCRS": "A, dA", Weight: 313.21 },
    { Name: "dC", BioSpring: "C", "Janssen GCRS": "C, dC", Weight: 289.18 },
    { Name: "dG", BioSpring: "G", "Janssen GCRS": "G, dG", Weight: 329.21 },
    { Name: "dT", BioSpring: "T", "Janssen GCRS": "T, dT", Weight: 304.2 },
    { Name: "rA", BioSpring: "", "Janssen GCRS": "rA", Weight: 329.21 },
    { Name: "rC", BioSpring: "", "Janssen GCRS": "rC", Weight: 305.18 },
    { Name: "rG", BioSpring: "", "Janssen GCRS": "rG", Weight: 345.21 },
    { Name: "rU", BioSpring: "", "Janssen GCRS": "rU", Weight: 306.17 }
  ]);
  const omeAndFluoroDf = DG.DataFrame.fromObjects([
    { Name: "2'-fluoro-U", BioSpring: "1", Axolabs: "Uf", "Janssen GCRS": "fU", Weight: 308.16 },
    { Name: "2'-fluoro-A", BioSpring: "2", Axolabs: "Af", "Janssen GCRS": "fA", Weight: 331.2 },
    { Name: "2'-fluoro-C", BioSpring: "3", Axolabs: "Cf", "Janssen GCRS": "fC", Weight: 307.18 },
    { Name: "2'-fluoro-G", BioSpring: "4", Axolabs: "Gf", "Janssen GCRS": "fG", Weight: 347.19 },
    { Name: "2'OMe-rU", BioSpring: "5", Axolabs: "u", "Janssen GCRS": "mU", Weight: 320.2 },
    { Name: "2'OMe-rA", BioSpring: "6", Axolabs: "a", "Janssen GCRS": "mA", Weight: 343.24 },
    { Name: "2'OMe-rC", BioSpring: "7", Axolabs: "c", "Janssen GCRS": "mC", Weight: 319.21 },
    { Name: "2'OMe-rG", BioSpring: "8", Axolabs: "g", "Janssen GCRS": "mG", Weight: 359.24 },
    { Name: "ps linkage", BioSpring: "*", Axolabs: "s", "Janssen GCRS": "ps", Weight: 16.07 }
  ]);
  const modificationsDf = await getModificationsDf();
  let modifGrid = DG.Viewer.grid(modificationsDf, { showRowHeader: false, showCellTooltip: false });
  modifGrid.col(OVERHANG_COL_NAMES.LONG_NAMES)!.width = 110;
  modifGrid.col(OVERHANG_COL_NAMES.ABBREVIATION)!.width = 80;
  modifGrid.col(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT)!.width = 105;
  modifGrid.col(OVERHANG_COL_NAMES.BASE_MODIFICATION)!.width = 110;
  modifGrid.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.width = 100;
  let view = grok.shell.newView('Oligo Batch Calculator', [
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
              inputSequences.root
            ],'inputSequence')
          ]), { style: { maxHeight: '230px' } }
        ),
        ui.splitV([
          title,
          gridDiv
        ])
      ]),
      ui.splitV([
        ui.box(ui.h2('ASO Gapmers'), { style: {maxHeight: '40px'} }),
        DG.Viewer.grid(asoGapmersDf!, { showRowHeader: false, showCellTooltip: false }).root,
        ui.box(ui.h2("2'-OMe and 2'-F modifications"), { style: {maxHeight: '40px'} }),
        DG.Viewer.grid(omeAndFluoroDf!, { showRowHeader: false, showCellTooltip: false }).root,
        ui.box(ui.h2('Overhang modifications'), { style: {maxHeight: '40px'} }),
        modifGrid.root
      ], { style: { maxWidth: '600px' } })
    ])
  ]);
  const adminUsers = ['Baozhong Zhao', 'Sijin Guo', 'Saika Siddiqui', 'Vadym Kovadlo'];
  view.box = true;
  let col = modifGrid.col('Action')!;
  col.cellType = 'html';
  modifGrid.onCellPrepare(function (gc) {
    if (gc.isTableCell && gc.gridColumn.name == OVERHANG_COL_NAMES.ACTION) {
      gc.style.element = ui.divH([
        ui.button(ui.iconFA('trash-alt'), () => {
          grok.dapi.users.current().then(async(user) => {
            if (adminUsers.includes(user.firstName + ' ' + user.lastName)) {
              modificationsDf.rows.removeAt(gc.gridRow, 1, true);
              await grok.dapi.userDataStorage.remove(STORAGE_NAME, modificationsDf.col(OVERHANG_COL_NAMES.LONG_NAMES)!.get(gc.gridRow), false);
            } else
              grok.shell.info('You don\'t have permission for this action');
          });
        }),
        ui.button(ui.iconFA('edit'), () => {
          grok.dapi.users.current().then((user) => {
            if (adminUsers.includes(user.firstName + ' ' + user.lastName)) {
              let oldLongName = modificationsDf.col(OVERHANG_COL_NAMES.LONG_NAMES)!.get(gc.gridRow);
              let longName = ui.stringInput(OVERHANG_COL_NAMES.LONG_NAMES, oldLongName);
              ui.tooltip.bind(longName.root, "Examples: 'Inverted Abasic', 'Cyanine 3 CPG', '5-Methyl dC'");
              let abbreviation = ui.stringInput(OVERHANG_COL_NAMES.ABBREVIATION, modificationsDf.col(OVERHANG_COL_NAMES.ABBREVIATION)!.get(gc.gridRow));
              ui.tooltip.bind(abbreviation.root, "Examples: 'invabasic', 'Cy3', '5MedC'");
              let molecularWeight = ui.floatInput(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, modificationsDf.col(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT)!.get(gc.gridRow));
              let baseModification = ui.choiceInput(OVERHANG_COL_NAMES.BASE_MODIFICATION, 'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
                if (v != 'NO')
                  extinctionCoefficient.value = 'Base';
                extinctionCoefficient.enabled = (v == 'NO');
              });
              let extinctionCoefficient = ui.stringInput(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, modificationsDf.col(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT)!.get(gc.gridRow));
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
                    return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters')
                  if (abbreviation.value.length > 100)
                    return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters')
                  await grok.dapi.userDataStorage.postValue(
                    STORAGE_NAME,
                    longName.value,
                    JSON.stringify({
                      longName: longName.value,
                      abbreviation: abbreviation.value,
                      molecularWeight: molecularWeight.value,
                      extinctionCoefficient: extinctionCoefficient.value,
                      baseModification: baseModification.value,
                    }),
                    false
                  );
                  if (oldLongName != longName.value)
                    await grok.dapi.userDataStorage.remove(STORAGE_NAME, oldLongName, false);
                  modificationsDf.set(OVERHANG_COL_NAMES.LONG_NAMES, gc.gridRow, longName.value);
                  modificationsDf.set(OVERHANG_COL_NAMES.ABBREVIATION, gc.gridRow, abbreviation.value);
                  modificationsDf.set(OVERHANG_COL_NAMES.MOLECULAR_WEIGHT, gc.gridRow, molecularWeight.value);
                  modificationsDf.set(OVERHANG_COL_NAMES.EXTINCTION_COEFFICIENT, gc.gridRow, extinctionCoefficient.value);
                  modificationsDf.set(OVERHANG_COL_NAMES.BASE_MODIFICATION, gc.gridRow, baseModification.value);
                })
                .show();
            } else
              grok.shell.info('You don\'t have permission for this action');
          });
        })
      ]);
    }
  });

  view.setRibbonPanels([[
    ui.iconFA('redo', () => inputSequences.value = ''),
    ui.iconFA('plus', () => addModificationButton()),
    ui.iconFA('arrow-to-bottom', () => saveAsCsv(table))
  ]]);

  $('.inputSequence textarea')
    .css('resize','none')
    .css('min-height','70px')
    .css('width','100%')
    .css('font-family','monospace')
    .attr("spellcheck", "false");
}
