/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {map} from "./map";

export let _package = new DG.Package();

const STORAGE_NAME = 'oligo-batch-calculator-storage';

let weightsObj: {[code: string]: number} = {};
let normalizedObj: {[code: string]: string} = {};
for (let synthesizer of Object.keys(map))
  for (let technology of Object.keys(map[synthesizer]))
    for (let code of Object.keys(map[synthesizer][technology])) {
      weightsObj[code] = map[synthesizer][technology][code].weight;
      normalizedObj[code] = map[synthesizer][technology][code].normalized;
    }

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
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
  let clearSequences = ui.button('CLEAR', () => inputSequences.value = '');

  updateTable(defaultInput);

  let saveAsButton = ui.bigButton('SAVE AS CSV', () => {
    let link = document.createElement("a");
    link.setAttribute("href", "data:text/csv;charset=utf-8,\uFEFF" + encodeURI(table.toCsv()));
    link.setAttribute("download", "Oligo Properties.csv");
    link.click();
  });

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

  let showCodesButton = ui.button('SHOW CODES', () => {
    ui.dialog('Codes')
      .add(tables)
      .show();
  });

  async function loadAll() {
    let modifications: any[] = [];
    let entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, false);
    if (entries !== null && Object.keys(entries).length === 0)
      grok.shell.info('Storage is empty. Try to post something to the storage');
    else {
      Object.keys(entries).forEach((key) => {
        modifications.push(JSON.parse(entries[key]));
      });
    }
    return modifications;
  }

  const baseModificationItems = ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'];
  let addModificationButton = ui.button('ADD MODIFICATION',  async() => {
    let modifications = await loadAll();
    let longName = ui.stringInput('Long name', '');
    ui.tooltip.bind(longName.root, "Examples: 'Inverted Abasic', 'Cyanine 3 CPG', '5-Methyl dC'");
    let abbreviation = ui.stringInput('Abbreviation', '');
    ui.tooltip.bind(abbreviation.root, "Examples: 'invabasic', 'Cy3', '5MedC'");
    let molecularWeight = ui.floatInput('Molecular weight', 0);
    let baseModification = ui.choiceInput('Base modification', 'NO', baseModificationItems, (v: string) => {
      if (v != 'NO')
        extinctionCoefficient.value = 'Base';
      extinctionCoefficient.enabled = (v == 'NO');
    });
    let extinctionCoefficient = ui.stringInput('Extinction coefficient', '');
    ui.dialog('Add Modification')
      .add(ui.block([
        DG.Viewer.grid(
          DG.DataFrame.fromColumns([
            DG.Column.fromStrings(longName.caption, modifications.map((e) => e.longName)),
            DG.Column.fromStrings(abbreviation.caption, modifications.map((e) => e.abbreviation)),    // @ts-ignore
            DG.Column.fromFloat32Array(molecularWeight.caption, modifications.map((e) => (e.molecularWeight == undefined) ? 0 : e.molecularWeight)),
            DG.Column.fromStrings(baseModification.caption, modifications.map((e) => e.baseModification)),
            DG.Column.fromStrings(extinctionCoefficient.caption, modifications.map((e) => e.extinctionCoefficient))
          ])
        ).root,
        longName.root,
        abbreviation.root,
        molecularWeight.root,
        baseModification.root,
        extinctionCoefficient.root,
      ]))
      .onOK(() => {
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
      })
      .showModal(true);
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
            units.root,
          ]),
          ui.h1('Input Sequences'),
          ui.div([
            inputSequences.root
          ],'inputSequence'),
          ui.divH([
            clearSequences,
            showCodesButton,
            addModificationButton,
          ])
        ]), {style:{maxHeight:'270px'}}
      ),
      ui.splitV([
        title,
        ui.panel([gridDiv], 'ui-box')
      ]),
      ui.box(
        ui.panel([
          saveAsButton,
        ]), {style:{maxHeight:'60px'}}
      )
    ])
  ]);
  view.box = true;

  $('.inputSequence textarea')
    .css('resize','none')
    .css('min-height','70px')
    .css('width','100%')
    .css('font-family','monospace')
    .attr("spellcheck", "false");
}