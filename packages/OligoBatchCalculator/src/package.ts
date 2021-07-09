/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: opticalDensity
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['NMole', 'Milligrams', 'Micrograms']}
//output: double opticalDensity
export function od(sequence: string, amount: number, outputUnits: string) {
  if (outputUnits == 'Milligrams' || outputUnits == 'Micrograms') {
    const coefficient = outputUnits == 'Milligrams' ? 1 : 0.001;
    return coefficient * amount * extinctionCoefficient(sequence) / molecularWeight(sequence);
  }
  let coefficient = (outputUnits == 'Milligrams') ? 1 : 1000;
  if (outputUnits == 'NMole') coefficient = 1000000;
  return amount * extinctionCoefficient(sequence) / coefficient;
}

//name: nMole
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micrograms']}
//output: double nMole
export function nMole(sequence: string, amount: number, outputUnits: string): number {
  return (outputUnits == 'Optical Density') ? amount * 1000000 / extinctionCoefficient(sequence) : amount * 1000 / molecularWeight(sequence);
}

//name: molecularMass
//input: string sequence
//input: double amount
//input: string outputUnits {choices: ['Optical Density', 'Milligrams', 'Micromoles', 'Millimoles']}
//output: double nMole
export function molecularMass(sequence: string, amount: number, outputUnits: string): number {
  if (outputUnits == 'Optical Density')
    return 1000 * amount / extinctionCoefficient(sequence) * molecularWeight(sequence);
  const coefficient = (outputUnits == 'Milligrams' || outputUnits == 'Micromoles') ? 1 : 1000;
  return amount / extinctionCoefficient(sequence) * molecularWeight(sequence) * coefficient * od(sequence, amount, outputUnits) / nMole(sequence, amount, outputUnits);
}

//name: molecularWeight
//input: string sequence
//output: double molecularWeight
export function molecularWeight(sequence: string): number {
  const weights: {[index: string]: number} = {
    "ps":	16.07,
    "fA":	331.2, "fU":	308.16, "fC":	307.18, "fG":	347.19,
    "mA":	343.24, "mU":	320.2, "mC":	319.21, "mG":	359.24,
    "A": 313.21, "U": 306.17, "C": 289.18, "G": 329.21, "T": 304.2,
    "dA": 313.21, "dU": 306.17, "dC": 289.18, "dG": 329.21, "dT": 304.2,
    "rA": 329.21, "rU": 306.17, "rC": 305.18, "rG": 345.21
  };
  const slicingStep = /^[AUGCT]/g.test(sequence) ? 1 : 2;
  let molecularWeight = 0;
  for (let i = 0; i < sequence.length; i += slicingStep)
    molecularWeight += weights[sequence.slice(i, i + slicingStep)];
  return (sequence.length > 0) ? molecularWeight - 61.97 : 0;
}

//name: extinctionCoefficient
//input: string sequence
//output: double ec
export function extinctionCoefficient(sequence: string) {
  let sequences = normalizeSequences([sequence]);
  sequence = sequences[0];
  const individualBases: {[index: string]: number} = {
    'dA': 15400, 'dC': 7400, 'dG': 11500, 'dT': 8700,
    'rA': 15400, 'rC': 7200, 'rG': 11500, 'rU': 9900
    },
    nearestNeighbour: any = {
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
        nearestNeighbour['r' + ((sequence[i + 1] == 'T') ? 'U' : sequence[i + 1])]['r' + ((sequence[i + 3] == 'T') ? 'U' : sequence[i + 3])]
        +
        nearestNeighbour['d' + ((sequence[i + 1] == 'U') ? 'T' : sequence[i + 1])]['d' + ((sequence[i + 3] == 'U') ? 'T' : sequence[i + 3])]
      ) / 2;
  for (let i = 2; i < sequence.length - 2; i += 2)
    ec2 += individualBases[sequence.slice(i, i + 2)];
  return ec1 - ec2;
}

function calculateNMole(molecularWeights: number[], extinctionCoefficients: number[], amount: number, units: string): number[] {
  let nmoles = Array(molecularWeights.length);
  if (units == 'nmole' || units == 'µmole') return nmoles.fill(amount);
  if (units == 'OD') {
    for (let i = 0; i < molecularWeights.length; i++)
      nmoles[i] = amount * 1000000 / extinctionCoefficients[i];
    return (molecularWeights[0] > 0) ? nmoles : Array(molecularWeights.length).fill(0);
  }
  for (let i = 0; i < molecularWeights.length; i++)
    nmoles[i] = amount * 1000 / molecularWeights[i];
  return (molecularWeights[0] > 0) ? nmoles : Array(molecularWeights.length).fill(0);
}

function calculateMass(extinctionCoefficients: number[], molecularWeights: number[], nmoles: number[], od260: number[], amount: number, units: string) {
  let mass = Array(molecularWeights.length);
  if (units == 'mg' || units == 'µg') return mass.fill(amount);
  if (units == 'OD') {
    for (let i = 0; i < molecularWeights.length; i++)
      mass[i] = 1000 * amount / extinctionCoefficients[i] * molecularWeights[i];
    return (molecularWeights[0] > 0) ? mass : Array(molecularWeights.length).fill(0);
  }
  const coefficient = (units == 'mg' || units == 'µmole') ? 1 : 1000;
  for (let i = 0; i < extinctionCoefficients.length; i++)
    mass[i] = amount / extinctionCoefficients[i] * molecularWeights[i] * coefficient * od260[i] / nmoles[i];
  return mass;
}

function opticalDensity(extinctionCoefficients: number[], molecularWeights: number[], nmoles: number[], amount: number, units: string) {
  let od = Array(molecularWeights.length);
  if (units == 'OD') return od.fill(amount);
  if (units == 'mg' || units == 'µg') {
    const coefficient = units == 'mg' ? 1 : 0.001;
    for (let i = 0; i < extinctionCoefficients.length; i++)
      od[i] = coefficient * amount * extinctionCoefficients[i] / molecularWeights[i];
    return od;
  }
  let coefficient = (units == 'mg') ? 1 : 1000;
  if (units == 'nmole') coefficient = 1000000;
  for (let i = 0; i < extinctionCoefficients.length; i++)
    od[i] = amount * extinctionCoefficients[i] / coefficient;
  return od;
}

function normalizeSequences(sequences: string[]): string[] {
  let normalizedSequences = Array(sequences.length);
  for (let i = 0; i < sequences.length; i++) {
    const isRna = (/^[AUGC]+$/.test(sequences[i]));
    const isDna = (/^[ATGC]+$/.test(sequences[i]));
    const obj: {[index: string]: string} = isRna ?
      {"A": "rA", "U": "rU", "G": "rG", "C": "rC"} :
      isDna ?
        {"A": "dA", "T": "dT", "G": "dG", "C": "dC"} :
        {"fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG", "mU": "rU", "mA": "rA", "mC": "rC", "mG": "rG", "ps": ""};
    normalizedSequences[i] = isRna ?
      sequences[i].replace(/[AUGC]/g, function (x) {return obj[x]}) :
      isDna ?
        sequences[i].replace(/[ATGC]/g, function (x) {return obj[x]}) :
        sequences[i].replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x) {return obj[x]});
  }
  return normalizedSequences;
}

function prepareInputTextField(text: string) {
  return text.split('\n').map((s) => s.replace(/\s/g, '')).filter(item => item);
}

//name: Oligo Batch Calculator
//tags: app
export function OligoBatchCalculator() {

  const defaultInput = 'fAmCmGmAmCpsmU\nmApsmApsfGmAmUmCfGfAfC\nmAmUfGmGmUmCmAfAmGmA';

  function updateTable(text: string) {
    let sequences = prepareInputTextField(text);
    let normalizedSequences = normalizeSequences(sequences);
    tableDiv.innerHTML = '';
    let molecularWeights = sequences.map((s) => molecularWeight(s));
    let extinctionCoefficients = normalizedSequences.map((s) => extinctionCoefficient(s));
    let nMole = calculateNMole(molecularWeights, extinctionCoefficients, yieldAmount.value, units.value);
    let od260 = opticalDensity(extinctionCoefficients, molecularWeights, nMole, yieldAmount.value, units.value);
    let mass = calculateMass(extinctionCoefficients, molecularWeights, nMole, od260, yieldAmount.value, units.value);

    let moleName1 = (units.value == 'µmole' || units.value == 'mg') ? 'µmole' : 'nmole';
    let moleName2 = (units.value == 'µmole') ? 'µmole' : 'nmole';
    let massName = (units.value == 'µmole') ? 'mg' : (units.value == 'mg') ? units.value : 'µg';
    const coefficient = (units.value == 'mg' || units.value == 'µmole') ? 1000 : 1;

    table = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'Item #', Array(...Array(sequences.length + 1).keys()).slice(1)),
      DG.Column.fromList('string', 'Sequence', sequences),
      DG.Column.fromList('int', 'Length', normalizedSequences.map((s) => s.length / 2)),
      DG.Column.fromList('double', 'OD 260', od260),
      DG.Column.fromList('double', moleName1, nMole),
      DG.Column.fromList('double', 'Mass (' + massName + ')', mass),
      DG.Column.fromList('double', moleName2 + '/OD', nMole.map(function(n, i) {return coefficient * n / od260[i]})),
      DG.Column.fromList('double', 'µg/OD', mass.map(function(n, i) {return coefficient * n / od260[i]})),
      DG.Column.fromList('double', 'MW', molecularWeights),
      DG.Column.fromList('int', 'Ext. Coefficient', extinctionCoefficients)
    ]);
    tableDiv.append(DG.Viewer.grid(table).root);
  }

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let text2 = ui.divText('Search Modifications');

  let yieldAmount = ui.floatInput('', 1, () => updateTable(inputSequenceField.value));
  let enter2 = ui.stringInput('', '');
  let units = ui.choiceInput('', 'OD', ['OD', 'µg', 'mg', 'µmole', 'nmole'], () => {
    updateTable(inputSequenceField.value);
  });
  let analyzeAsSingleStrand = ui.button('Analyze As Single Strand', () => grok.shell.info('Coming soon'));
  let analyzeAsDuplexStrand = ui.button('Analyze As Duplex Strand', () => grok.shell.info('Coming soon'));
  let clearSequences = ui.button('CLEAR', () => inputSequenceField.value = '');
  let addModifications = ui.button('Add Modifications', () => grok.shell.info('Coming soon'));
  let threeMod = ui.boolInput("3' MOD", false);
  let internal = ui.boolInput('INTERNAL', false);
  let fiveMod = ui.boolInput("5' MOD", false);
  let table = DG.DataFrame.create();
  let inputSequenceField = ui.textInput("", defaultInput, async (txt: string) => updateTable(txt));

  let tableDiv = ui.box();
  updateTable(defaultInput);

  let saveAsButton = ui.bigButton('SAVE AS CSV', () => {
    let csvContent = table.toCsv();
    let encodedUri = encodeURI(csvContent);
    let link = document.createElement("a");
    link.setAttribute("href", "data:text/csv;charset=utf-8,\uFEFF" + encodedUri);
    link.setAttribute("download", "Oligo Properties.csv");
    link.click();
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
            inputSequenceField.root
          ],'inputSequence'),
          clearSequences,
        ]), {style:{maxHeight:'245px'}}
      ),
      ui.splitV([
        title,
        ui.panel([tableDiv], 'ui-box')
      ]),
      ui.box(
        ui.panel([saveAsButton]), {style:{maxHeight:'60px'}}
      ),
    ])
  ]);
  view.box = true;

  $('.inputSequence textarea')
    .css('resize','none')
    .css('min-height','50px')
    .css('width','100%')
    .css('font-family','monospace');
}