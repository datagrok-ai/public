/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

import {map} from "./map";

export let _package = new DG.Package();

let molecularWeightsObj: {[index: string]: number} = {};
for (let representation of Object.values(map))
  for (let code of Object.keys(representation))
    molecularWeightsObj[code] = representation[code].molecularWeight;

function sortByStringLengthInDescendingOrderToCheckForMatchWithLongerCodesFirst(array: string[]): string[] {
  return array.sort(function(a, b) { return b.length - a.length; });
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
  if (outputUnits == 'Milligrams' || outputUnits == 'Micrograms' || outputUnits == 'mg' || outputUnits == 'µg') {
    const coefficient = outputUnits == 'Milligrams' ? 1 : 0.001;
    return coefficient * amount * extinctionCoefficient(sequence) / molecularWeight(sequence);
  } else if (outputUnits == 'OD') {
    return amount;
  }
  let coefficient = (outputUnits == 'NMole') ? 1000000 : (outputUnits == 'Milligrams') ? 1 : 1000;
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
export function molecularWeight(sequence: string, representation?: string): number {
  const codes = (representation) ?
    Object.keys(molecularWeightsObj[representation]) :
    sortByStringLengthInDescendingOrderToCheckForMatchWithLongerCodesFirst(Object.keys(molecularWeightsObj));
  let molecularWeight = 0;
  let i = 0;
  while (i < sequence.length) {
    let matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length));
    molecularWeight += molecularWeightsObj[sequence.slice(i, i + matchedCode!.length)];
    i += matchedCode!.length;
  }
  return molecularWeight - 61.97;
}

//name: extinctionCoefficient
//input: string sequence
//output: double extinctionCoefficient
export function extinctionCoefficient(sequence: string): number {
  sequence = normalizeSequence(sequence);
  const individualBases: {[index: string]: number} = {
      'dA': 15400, 'dC': 7400, 'dG': 11500, 'dT': 8700,
      'rA': 15400, 'rC': 7200, 'rG': 11500, 'rU': 9900
    },
    nearestNeighbour: {[index: string]: {[index: string]: number}} = {
      'dA': {'dA': 27400, 'dC': 21200, 'dG': 25000, 'dT': 22800},
      'dC': {'dA': 21200, 'dC': 14600, 'dG': 18000, 'dT': 15200},
      'dG': {'dA': 25200, 'dC': 17600, 'dG': 21600, 'dT': 20000},
      'dT': {'dA': 23400, 'dC': 16200, 'dG': 19000, 'dT': 16800},
      'rA': {'rA': 27400, 'rC': 21000, 'rG': 25000, 'rU': 24000},
      'rC': {'rA': 21000, 'rC': 14200, 'rG': 17800, 'rU': 16200},
      'rG': {'rA': 25200, 'rC': 17400, 'rG': 21600, 'rU': 21200},
      'rU': {'rA': 24600, 'rC': 17200, 'rG': 20000, 'rU': 19600}
    };
  let ec1 = 0, ec2 = 0;
  for (let i = 0; i < sequence.length - 2; i += 2)
    if (sequence[i] == sequence[i + 2])
      ec1 += nearestNeighbour[sequence.slice(i, i + 2)][sequence.slice(i + 2, i + 4)];
    else
      ec1 += (
        nearestNeighbour['r' + (sequence[i + 1] == 'T') ? 'U' : sequence[i + 1]]['r' + (sequence[i + 3] == 'T') ? 'U' : sequence[i + 3]]
        +
        nearestNeighbour['d' + (sequence[i + 1] == 'U') ? 'T' : sequence[i + 1]]['d' + (sequence[i + 3] == 'U') ? 'T' : sequence[i + 3]]
      ) / 2;
  for (let i = 2; i < sequence.length - 2; i += 2)
    ec2 += individualBases[sequence.slice(i, i + 2)];
  return ec1 - ec2;
}

function normalizeSequence(sequence: string): string {
  //TODO: get normalization strings from map.ts (problem: some keys are in several representations at the same time, example: "8", "A", etc.)
  //TODO: create searchValue for replace functions from map.ts
  //TODO: ask conventional name of this operation(instead of 'normalize') and export this function
  //TODO: pass representation argument
  const isNormalized = /^[rdAUGCT]+$/.test(sequence);
  const isRna = /^[AUGC]+$/.test(sequence);
  const isDna = /^[ATGC]+$/.test(sequence);
  if (isNormalized && !isDna && !isRna)
    return sequence;
  const isSiRnaAxolabs = /^[fAUGCuacgs]+$/.test(sequence);
  const isGCRS = /^[fmpsACGU]+$/.test(sequence);
  const isGcrsGapmers = /^.*moe.+$/.test(sequence) || /^.*5mC+$/.test(sequence);
  const obj: {[index: string]: string} = isRna ?
    {"A": "rA", "U": "rU", "G": "rG", "C": "rC"} :
    isDna ?
      {"A": "dA", "T": "dT", "G": "dG", "C": "dC"} :
      isSiRnaAxolabs ?
        {"Af": "rA", "Uf": "rU", "Gf": "rG", "Cf": "rC", "u": "rU", "a": "rA", "c": "rC", "g": "rG", "s": "", "fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG"} :
        isGCRS ?
          {"fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG", "mU": "rU", "mA": "rA", "mC": "rC", "mG": "rG", "ps": "", "s": ""} :
          isGcrsGapmers ?
            {"moeA": "rA", "(5m)moeC": "rC", "moe5mC": "rC", "moeG": "rG", "moeT": "rU", "(5m)C": "rC", "5mC": "rC", "U": "rU", "T": "rU", "A": "rA", "C": "rC", "G": "rG", "fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG", "mU": "rU", "mA": "rA", "mC": "rC", "mG": "rG", "ps": "", "s": ""} :
            {"ps": "", "mA": "rA", "mU": "rU", "mG": "rG", "mC": "rC", "fA": "rA", "fU": "rU", "fG": "rG", "fC": "rC"};

  return isRna ?
    sequence.replace(/[AUGC]/g, function (x) {return obj[x]}) :
    isDna ?
      sequence.replace(/[ATGC]/g, function (x) {return obj[x]}) :
      isSiRnaAxolabs ?
        sequence.replace(/(Uf|Af|Cf|Gf|fU|fA|fC|fG|u|a|c|g|s)/g, function (x) {return obj[x]}) :
        isGCRS ?
          sequence.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps|s)/g, function (x) {return obj[x]}) :
          isGcrsGapmers ?
            sequence.replace(/(moeA|\(5m\)moeC|moe5mC|moeG|moeT|A|T|\(5m\)C|5mC|G|ps|s)/g, function (x) {return obj[x]}) :
            sequence.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps|A|U|G|C)/g, function (x) {return obj[x]});
}

function getListOfPossibleRepresentationsByFirstMatchedCodes(sequence: string): string[] {
  let representations: string[] = [];
  Object.keys(map).forEach((representation: string) => {
    if (Object.keys(map[representation]).some((s) => s == sequence.slice(0, s.length)))
      representations.push(representation);
  });
  return representations;
}

function isValid(sequence: string) {
  //TODO: rewrite using switch case
  let possibleRepresentations = getListOfPossibleRepresentationsByFirstMatchedCodes(sequence);
  if (possibleRepresentations.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedRepresentation: null };

  let outputIndices = Array(possibleRepresentations.length).fill(0);

  const firstUniqueCharacters = ['r', 'd', 'f', 'm'],
    secondCharactersInNormalizedCodes = ["A", "U", "T", "C", "G"];

  possibleRepresentations.forEach((representation, representationIndex) => {
    while (outputIndices[representationIndex] < sequence.length) {

      let matchedCode = Object.keys(map[possibleRepresentations[representationIndex]])
        .find((s) => s == sequence.slice(outputIndices[representationIndex], outputIndices[representationIndex] + s.length));

      if (matchedCode == null)
        break;

      if (  // for mistake pattern 'rAA'
        outputIndices[representationIndex] > 1 &&
        secondCharactersInNormalizedCodes.includes(sequence[outputIndices[representationIndex]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[representationIndex] - 2])
      )
        break;

      if (  // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[representationIndex] + 1]) &&
        secondCharactersInNormalizedCodes.includes(sequence[outputIndices[representationIndex]])
      ) {
        outputIndices[representationIndex]++;
        break;
      }

      outputIndices[representationIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedRepresentation = Math.max.apply(Math, outputIndices);
  const indexOfFirstNotValidCharacter = (indexOfExpectedRepresentation == sequence.length) ? -1 : indexOfExpectedRepresentation;

  return {
    indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
    expectedRepresentation: possibleRepresentations[outputIndices.indexOf(indexOfFirstNotValidCharacter)]
  };
}

//name: Oligo Batch Calculator
//tags: app
export function OligoBatchCalculatorApp() {

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
      expectedRepresentations = Array(sequences.length);

    sequences.forEach((sequence, i) => {
      let output = isValid(sequence);
      indicesOfFirstNotValidCharacter[i] = output.indexOfFirstNotValidCharacter;
      expectedRepresentations[i] = output.expectedRepresentation;
      if (indicesOfFirstNotValidCharacter[i] < 0) {
        normalizedSequences[i] = normalizeSequence(sequence);
        if (normalizedSequences[i].length > 2) {
          try {
            molecularWeights[i] = molecularWeight(sequence, expectedRepresentations[i]);
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
      } else if (expectedRepresentations[i] == null)
        reasonsOfError[i] = "Not valid input";
      else
        reasonsOfError[i] = "Sequence is expected to be in representation '" +  expectedRepresentations[i] +
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

    let grid = DG.Viewer.grid(table);
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
  let clearSequences = ui.button('CLEAR', () => inputSequences.value = '');

  updateTable(defaultInput);

  let saveAsButton = ui.bigButton('SAVE AS CSV', () => {
    let link = document.createElement("a");
    link.setAttribute("href", "data:text/csv;charset=utf-8,\uFEFF" + encodeURI(table.toCsv()));
    link.setAttribute("download", "Oligo Properties.csv");
    link.click();
  });

  let tables = ui.divV([]);
  for (let representation of Object.keys(map)) {
    let tableRows = [];
    for (let [key, value] of Object.entries(map[representation]))
      tableRows.push({'name': value.name, 'code': key, 'weight': value['molecularWeight']});
    tables.append(
      DG.HtmlTable.create(
        tableRows,
        (item: {name: string; code: string; weight: number}) => [item['name'], item['code'], item['weight']],
        [representation, 'Code', 'Weight']
      ).root,
      ui.div([], {style: {height: '30px'}})
    );
  }

  let showCodesButton = ui.button('SHOW CODES', () => {
    ui.dialog('Codes')
      .add(tables)
      .show();
  });

  let title = ui.panel([ui.h1('Oligo Properties')], 'ui-panel ui-box');
  title.style.maxHeight = '45px';
  $(title).children('h1').css('margin', '0px');

  let view = grok.shell.newView('Oligo Batch Calculator', [
    ui.splitV([
      ui.box(
        ui.panel([
          ui.h1('Yield Amount & Units'),
          ui.divH([
            yieldAmount.root,
            units.root
          ]),
          ui.h1('Input Sequences'),
          ui.div([
            inputSequences.root
          ],'inputSequence'),
          clearSequences,
        ]), {style:{maxHeight:'270px'}}
      ),
      ui.splitV([
        title,
        ui.panel([gridDiv], 'ui-box')
      ]),
      ui.box(
        ui.panel([
          ui.divH([
            saveAsButton,
            showCodesButton
          ])
        ]), {style:{maxHeight:'60px'}}
      )
    ])
  ]);
  view.box = true;

  $('.inputSequence textarea')
    .css('resize','none')
    .css('min-height','70px')
    .css('width','100%')
    .css('font-family','monospace');
}