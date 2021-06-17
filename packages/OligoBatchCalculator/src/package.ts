/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

function normalizeSeque

function calculateExtinctionCoefficientUsingNearestNeighborMethod(sequence: string) {
  let extinctionCoefficient = -1;
  const searchValue = '/(' + Object.keys(weights) + ')/g';
  // @ts-ignore
  sequence.replace(searchValue, function (x: string) {

  });
}

function calculateMolecularWeight(sequence: string): number {
  let molecularWeight = 0;
  const searchValue = '/(' + Object.keys(weights).join('|') + ')/g';
  // @ts-ignore
  sequence.replace("/(1|2|3|4|5|6|7|8|9|moeT|moeA|moe5mC|moeG|5mC|A|C|G|T|Uf|fU|Af|fA|Cf|fC|Gf|fG|u|mU|a|mA|c|mC|g|mG|*|s|ps)/g",  (nucleotideSymbol: string) => molecularWeight += weights[nucleotideSymbol]);
  return molecularWeight;
}

const weights: {[index: string]: number} = {
  "5": 378.27,
  "moeT": 378.27,
  "6": 387.29,
  "moeA": 387.29,
  "7": 377.29,
  "moe5mC": 377.29,
  "8": 403.28,
  "moeG": 403.28,
  "9": 303.21,
  "5mC": 303.21,
  "A": 313.21,
  "C": 289.18,
  "G": 329.21,
  "T": 304.20,
  "1": 308.16,
  "Uf": 308.16,
  "fU": 308.16,
  "2": 331.20,
  "Af": 331.20,
  "fA": 331.20,
  "3": 307.18,
  "Cf": 307.18,
  "fC": 307.18,
  "4": 347.19,
  "Gf": 347.19,
  "fG": 347.19,
  // "5": 320.20,
  "u": 320.20,
  "mU": 320.20,
  // "6": 343.24,
  "a": 343.24,
  "mA": 343.24,
  // "7": 319.21,
  "c": 319.21,
  "mC": 319.21,
  // "8": 359.24,
  "g": 359.24,
  "mG": 359.24,
  "*": 16.07,
  "s": 16.07,
  "ps": 16.07
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
  individualRnaBases = {
    'rA': 15400,
    'rC': 7200,
    'rG': 11500,
    'rU': 9900
  },
  nearestNeighbourRna = {
    'dA': {'dA': 27400, 'dC': 21200, 'dG': 25000, 'dT': 22800},
    'dC': {'dA': 21200, 'dC': 14600, 'dG': 18000, 'dT': 15200},
    'dG': {'dA': 25200, 'dC': 17600, 'dG': 21600, 'dT': 20000},
    'dT': {'dA': 23400, 'dC': 16200, 'dG': 19000, 'dT': 16800}
  },
  nearestNeighbourDna = {
    'rA': {'rA': 27400, 'rC': 21000, 'rG': 25000, 'rU': 24000},
    'rC': {'rA': 21000, 'rC': 14200, 'rG': 17800, 'rU': 16200},
    'rG': {'rA': 25200, 'rC': 17400, 'rG': 21600, 'rU': 21200},
    'rU': {'rA': 24600, 'rC': 17200, 'rG': 20000, 'rU': 19600}
  };

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
  let inputSequenceField = ui.textInput("", "", async (seq: string) => {
    let sequences = seq.split('\n');
    let tableRows = [];
    for (let i = 0; i < sequences.length; i++) {
      // @ts-ignore
      tableRows.push({
        'n': i + 1,
        'seq': sequences[i],
        'len': sequences[i].length,
        'od': '1.00 OD',
        'nmole': 'ex',
        'mass': calculateMolecularWeight(sequences[i]),
        'nmoleOD': 'nm',
        'mgod': 'k',
        'mw': 'mw',
        'g': 'g',
        'gc': 'gc',
        'ec': 'ec'
      });
    }
    tableDiv.innerHTML = '';
    tableDiv.append(
      DG.HtmlTable.create(
        tableRows,
        (item: {
          n: string; seq: string; len: string; od: string; nmole: string; mass: string, nmoleOD: string; mgod: string; mw: string; g: string; gc: string; ec: string;
        }) => [item.n, item.seq, item.len, item.od, item.nmole, item.mass, item.nmoleOD, item.mgod, item.mw, item.g, item.gc, item.ec],
        ['Item#', 'Sequence', 'Length', 'OD-260', 'nmole', 'Mass', 'nmole/OD', 'μg/OD', 'MW', 'G%', 'GC%', 'Ext. Coefficient']
      ).root
    );
  });

  let tableRows: never[] = [];
  let tableDiv = ui.div([]);

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
        ui.h1('Input sequence'),
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
