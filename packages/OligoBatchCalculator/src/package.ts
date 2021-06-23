/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

function calculateMolecularWeight(sequence: string): number {
  const step = (/^[AUGC]/g.test(sequence)) ? 1 : 2;
  let molecularWeight = 0;
  for (let i = 0; i < sequence.length; i += step)
    molecularWeight += weights[sequence.slice(i, i + step)];
  return (sequence.length > 0) ? molecularWeight - 61.97 : 0;
}

function calculateNMole(molecularWeights: number[], extinctionCoefficients: number[], amount: number, units: string): number[] {
  let nmoles = Array(molecularWeights.length);
  if (units == 'nmole' || units == 'µmole') return nmoles.fill(amount);
  if (units == 'OD') {
    for (let i = 0; i < molecularWeights.length; i++)
      nmoles[i] = amount * 1000000 / extinctionCoefficients[i];
    return (molecularWeights[0] > 0) ? nmoles : Array(molecularWeights.length).fill(0);
  }
  const coefficient = units == 'mg' ? 1000000 : 1000000000;
  for (let i = 0; i < molecularWeights.length; i++)
    nmoles[i] = amount * coefficient / molecularWeights[i]
  return (molecularWeights[0] > 0) ? nmoles : Array(molecularWeights.length).fill(0);
}

function calculateMass(extinctionCoefficients: number[], molecularWeights: number[], nmoles: number[], amount: number, units: string) {
  let mass = Array(molecularWeights.length);
  if (units == 'mg' || units == 'µg') return mass.fill(amount);
  if (units == 'OD') {
    for (let i = 0; i < molecularWeights.length; i++)
      mass[i] = amount / extinctionCoefficients[i] * molecularWeights[i];
    return (molecularWeights[0] > 0) ? mass : Array(molecularWeights.length).fill(0);
  }
  const coefficient = units == 'mg' ? 1 : 1000;
  for (let i = 0; i < extinctionCoefficients.length; i++)
    mass[i] = amount / extinctionCoefficients[i] * molecularWeights[i] * coefficient;
  return mass;
}

function calculateOd(extinctionCoefficients: number[], molecularWeights: number[], nmoles: number[], amount: number, units: string) {
  let od = Array(molecularWeights.length);
  if (units == 'OD') return od.fill(amount);
  if (units == 'mg' || units == 'µg') {
    for (let i = 0; i < extinctionCoefficients.length; i++)
      od[i] = amount * extinctionCoefficients[i] / molecularWeights[i];
    return od;
  }
  const coefficient = units == 'mg' ? 1 : 1000;
  for (let i = 0; i < extinctionCoefficients.length; i++)
    od[i] = amount * extinctionCoefficients[i] / coefficient;
  return od;
}

function normalizeSequence(sequence: string): string {
  const isRna = (/^[AUGC]/g.test(sequence));
  const obj: {[index: string]: string} = isRna ?
    {"A": "rA", "U": "rU", "G": "rG", "C": "rC"} :
    {"fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG", "mU": "rU", "mA": "rA", "mC": "rC", "mG": "rG", "ps": ""};
  if (isRna) return sequence.replace(/[AUGC]/g, function (x) {return obj[x]});
  return sequence.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x) {return obj[x]});
}

function prepareInputTextField(text: string) {
  return text.split('\n');
}

const weights: {[index: string]: number} = {
  "fU":	308.16,
  "fA":	331.2,
  "fC":	307.18,
  "fG":	347.19,
  "mU":	320.2,
  "mA":	343.24,
  "mC":	319.21,
  "mG":	359.24,
  "ps":	16.07,
  "A": 313.21,
  "U": 306.17,
  "G": 329.21,
  "C": 289.18
};

//name: OligoBatchCalculator
//tags: app
export function OligoBatchCalculator() {

  const individualDnaBases = {
    'dA': 15400,
    'dC': 7400,
    'dG': 11500,
    'dT': 8700
  },
  individualRnaBases: any = {
    'rA': 15400,
    'rC': 7200,
    'rG': 11500,
    'rU': 9900
  },
  nearestNeighbourDna = {
    'dA': {'dA': 27400, 'dC': 21200, 'dG': 25000, 'dT': 22800},
    'dC': {'dA': 21200, 'dC': 14600, 'dG': 18000, 'dT': 15200},
    'dG': {'dA': 25200, 'dC': 17600, 'dG': 21600, 'dT': 20000},
    'dT': {'dA': 23400, 'dC': 16200, 'dG': 19000, 'dT': 16800}
  },
  nearestNeighbourRna: any = {
    'rA': {'rA': 27400, 'rC': 21000, 'rG': 25000, 'rU': 24000},
    'rC': {'rA': 21000, 'rC': 14200, 'rG': 17800, 'rU': 16200},
    'rG': {'rA': 25200, 'rC': 17400, 'rG': 21600, 'rU': 21200},
    'rU': {'rA': 24600, 'rC': 17200, 'rG': 20000, 'rU': 19600}
  },
  defaultInput = 'fAmCmGmAmCpsmU\nmApsmApsfGmAmUmCfGfAfC\nmAmUfGmGmUmCmAfAmGmA';

  function getExtinctionCoefficientUsingNearestNeighborMethod(sequence: string) {
    sequence = normalizeSequence(sequence);
    let ec1 = 0, ec2 = 0;
    for (let i = 0; i < sequence.length - 2; i += 2)
      ec1 += nearestNeighbourRna[sequence.slice(i, i + 2)][sequence.slice(i + 2, i + 4)];
    for (let i = 2; i < sequence.length - 4; i += 2)
      ec2 += individualRnaBases[sequence.slice(i, i + 2)];
    return ec1 - ec2;
  }

  function updateTable(text: string) {
    let sequences = prepareInputTextField(text);
    tableDiv.innerHTML = '';
    let molecularWeights = sequences.map((s) => calculateMolecularWeight(s));
    let extinctionCoefficients = sequences.map((s) => getExtinctionCoefficientUsingNearestNeighborMethod(s));
    let nMole = calculateNMole(molecularWeights, extinctionCoefficients, yieldAmount.value, units.value);
    let mass = calculateMass(extinctionCoefficients, molecularWeights, nMole, yieldAmount.value, units.value);
    let od260 = calculateOd(extinctionCoefficients, molecularWeights, nMole, yieldAmount.value, units.value);

    let moleName = (units.value == 'µmole') ? 'µmole' : 'nmole';
    let massName = (units.value == 'µg') ? units.value : 'mg';

    table = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'Item#', [...new Int16Array(sequences.length + 1).keys()].slice(1)),
      DG.Column.fromList('string', 'Sequence', sequences),
      DG.Column.fromList('int', 'Length', sequences.map((s) => s.length)),
      DG.Column.fromList('int', 'OD-260', od260),
      DG.Column.fromList('double', moleName, nMole),
      DG.Column.fromList('double', 'Mass (' + massName + ')', mass),
      DG.Column.fromList('double', moleName + '/OD', nMole.map(function(n, i) {return n / od260[i]})),
      DG.Column.fromList('double', massName + '/OD', mass.map(function(n, i) {return n / od260[i]})),
      DG.Column.fromList('double', 'MW', molecularWeights),
      DG.Column.fromList('int', 'Ext. Coefficient', extinctionCoefficients),
    ]);
    tableDiv.append(DG.Viewer.grid(table).root);
  }

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let text1 = ui.h1('Yield Amount & Units');
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

  let tableDiv = ui.block([]);
  updateTable(defaultInput);

  let loadFileButton = ui.button('LOAD AS CSV', () => {grok.shell.info('Coming soon')});

  grok.shell.newView('Oligo Batch Calculator', [
    ui.divH([
      ui.divV([
        text1,
        ui.divH([
          yieldAmount.root,
          units.root,
          clearSequences
        ])
      ]),
      // analyzeAsSingleStrand,
      // analyzeAsDuplexStrand,
      // addModifications,
      // ui.divV([
      //   text2,
      //   enter2.root
      // ]),
      // threeMod.root,
      // internal.root,
      // fiveMod.root
    ]),
    ui.block([
      ui.div([
        ui.h1('Input sequences'),
        ui.div([
          inputSequenceField.root
        ],'input-base')
      ], 'sequenceInput'),
    ]),
    loadFileButton,
    tableDiv
  ]);

  $('.sequence')
    .children().css('padding','5px 0');
  $('.sequenceInput .input-base').css('margin','0');
  $('.sequenceInput textarea')
    .css('resize','none')
    .css('min-height','50px')
    .css('width','100%');
  $('.sequenceInput select')
    .css('width','100%');
}