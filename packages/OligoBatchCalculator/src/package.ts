/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

function calculateMolecularWeight(sequence: string): number {
  let molecularWeight = 0;
  for (let i = 0; i < sequence.length; i += 2)
    molecularWeight += weights[sequence.slice(i, i + 2)];
  return molecularWeight - 61.97;
}

function normalizeSequence(sequence: string): string {
  const obj: {[index: string]: string} = {
    "fU": "rU", "fA": "rA", "fC": "rC", "fG": "rG", "mU": "rU", "mA": "rA", "mC": "rC", "mG": "rG", "ps": ""
  };
  return sequence.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
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
  "ps":	16.07
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
  };

  function getExtinctionCoefficientUsingNearestNeighborMethod(sequence: string) {
    sequence = normalizeSequence(sequence);
    let ec1 = 0, ec2 = 0;
    for (let i = 0; i < sequence.length - 2; i += 2)
      ec1 += nearestNeighbourRna[sequence.slice(i, i + 2)][sequence.slice(i + 2, i + 4)];
    for (let i = 2; i < sequence.length - 4; i += 2)
      ec2 += individualRnaBases[sequence.slice(i, i + 2)];
    return ec1 - ec2;
  }

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let text1 = ui.divText('Enter Yield Amount & Units');
  let text2 = ui.divText('Search Modifications');

  let enter1 = ui.stringInput('', '');
  let enter2 = ui.stringInput('', '');
  let choose = ui.choiceInput('', 'Select', []);
  let analyzeAsSingleStrand = ui.button('Analyze As Single Strand', () => grok.shell.info('Coming soon'));
  let analyzeAsDuplexStrand = ui.button('Analyze As Duplex Strand', () => grok.shell.info('Coming soon'));
  let clearSequences = ui.button('Clear Sequences', () => grok.shell.info('Coming soon'));
  let addModifications = ui.button('Add Modifications', () => grok.shell.info('Coming soon'));
  let threeMod = ui.boolInput("3' MOD", false);
  let internal = ui.boolInput('INTERNAL', false);
  let fiveMod = ui.boolInput("5' MOD", false);
  let table = DG.DataFrame.create();
  let inputSequenceField = ui.textInput("", "", async (seq: string) => {
    let sequences = seq.split('\n');
    tableDiv.innerHTML = '';
    let od260 = Array(sequences.length).fill(1);
    let nMole = Array(sequences.length).fill(19);
    let mass = Array(sequences.length).fill(30);
    table = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'Item#', [...Array(sequences.length + 1).keys()].slice(1)),
      DG.Column.fromList('string', 'Sequence', sequences),
      DG.Column.fromList('int', 'Length', sequences.map((s) => s.length)),
      DG.Column.fromList('int', 'OD-260', od260),
      DG.Column.fromList('double', 'nmole', nMole),
      DG.Column.fromList('double', 'Mass', mass),
      DG.Column.fromList('double', 'nmole/OD', nMole.map(function(n, i) {return n / od260[i];})),
      DG.Column.fromList('double', 'μg/OD', mass.map(function(n, i) {return n / od260[i];})),
      DG.Column.fromList('double', 'MW', sequences.map((s) => calculateMolecularWeight(s))),
      DG.Column.fromList('int', 'Ext. Coefficient', sequences.map((s) => getExtinctionCoefficientUsingNearestNeighborMethod(s))),
    ]);
    tableDiv.append(DG.Viewer.grid(table).root);
  });

  let tableDiv = ui.block([]);

  grok.shell.newView('Sequence Translator', [
    ui.divH([
      ui.divV([
        text1,
        ui.divH([
          enter1.root,
          choose.root
        ])
      ]),
      analyzeAsSingleStrand,
      analyzeAsDuplexStrand,
      clearSequences,
      addModifications,
      ui.divV([
        text2,
        enter2.root
      ]),
      threeMod.root,
      internal.root,
      fiveMod.root
    ]),
    ui.block([
      ui.div([
        ui.h1('Input sequences'),
        ui.div([
          inputSequenceField.root
        ],'input-base')
      ], 'sequenceInput'),

    ]),
    ui.button('Load Data To Excel', () => {}),
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