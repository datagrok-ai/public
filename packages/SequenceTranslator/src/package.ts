/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import $ from "cash-dom";
import {defineAxolabsPattern} from "./defineAxolabsPattern";
import {map} from "./map";

export let _package = new DG.Package();

const defaultInput = "AGGTCCTCTTGACTTAGGCC";
const minimalValidNumberOfCharacters = 6;
const smallNumberOfCharacters = "Length of input sequence should be at least " + minimalValidNumberOfCharacters + " characters";
const undefinedInputSequence = "Type of input sequence is undefined";
const noTranslationTableAvailable = "No translation table available";
const sequenceWasCopied = 'Copied';
const tooltipSequence = 'Copy sequence';

//name: Sequence Translator
//tags: app
export function sequenceTranslator() {

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  function updateTableAndSVG(sequence: string) {
    moleculeSvgDiv.innerHTML = "";
    outputTableDiv.innerHTML = "";
    let outputSequenceObj = convertSequence(sequence);
    let tableRows = [];
    for (let key of Object.keys(outputSequenceObj).slice(1)) {
      //@ts-ignore
      tableRows.push({'key': key, 'value': ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key]).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')})
    }
    outputTableDiv.append(
      ui.div([
        DG.HtmlTable.create(
          tableRows, (item: { key: string; value: string; }) => [item.key, item.value], ['Code', 'Sequence']
        ).root
      ], 'table')
    );
    semTypeOfInputSequence.textContent = 'Detected input type: ' + outputSequenceObj.type;
    if (!(outputSequenceObj.type == undefinedInputSequence || outputSequenceObj.type == smallNumberOfCharacters)) {
      let pi = DG.TaskBarProgressIndicator.create('Rendering molecule...');
      try {
        let flavor: string = (outputSequenceObj.Nucleotides.includes('U')) ? "RNA_both_caps" : "DNA_both_caps";
        (async () => {
          let smiles = (isSiRnaGcrsCode(inputSequenceField.value.replace(/\s/g, ''))) ? 
            modifiedToSmiles(inputSequenceField.value.replace(/\s/g, '')) :
            await nucleotidesToSmiles(outputSequenceObj.Nucleotides, flavor); 
          smiles = smiles.replace(/@/g, ''); // Remove StereoChemistry on the Nucleic acid chain and remove the Chiral label
          moleculeSvgDiv.append(grok.chem.svgMol(smiles, 900, 300));
        })();
      } finally {
        pi.close();
      }
    }
  }

  const appMainDescription = ui.info([
      ui.divText('\n How to convert one sequence:',{style:{'font-weight':'bolder'}}),
      ui.divText("Paste sequence into the text field below"),
      ui.divText('\n How to convert many sequences:',{style:{'font-weight':'bolder'}}),
      ui.divText("1. Drag & drop an Excel or CSV file with sequences into Datagrok. The platform will automatically detect columns with sequences"),
      ui.divText('2. Right-click on the column header, then see the \'Convert\' menu'),
      ui.divText("This will add the result column to the right of the table"),
    ], 'Convert oligonucleotide sequences between Nucleotides, BioSpring, Axolabs, and GCRS representations.'
  );

  let inputSequenceField = ui.textInput("", defaultInput, (sequence: string) => updateTableAndSVG(sequence));
  let outputSequencesObj = convertSequence(defaultInput);
  let semTypeOfInputSequence = ui.divText('Detected input type: ' + outputSequencesObj.type);

  let outputTableDiv = ui.div([
    DG.HtmlTable.create([
      {key: 'Nucleotides', value: ui.link(defaultInput, () => navigator.clipboard.writeText(defaultInput).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')},
      {key: 'BioSpring', value: ui.link(asoGapmersNucleotidesToBioSpring(defaultInput), () => navigator.clipboard.writeText(asoGapmersNucleotidesToBioSpring(defaultInput)).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')},
      {key: 'Axolabs', value: ui.link(noTranslationTableAvailable, () => navigator.clipboard.writeText(defaultInput).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')},
      {key: 'GCRS', value: ui.link(asoGapmersNucleotidesToGcrs(defaultInput), () => navigator.clipboard.writeText(asoGapmersNucleotidesToGcrs(defaultInput)).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')}
    ], (item: {key: string; value: string;}) => [item.key, item.value], ['Code', 'Sequence']).root
  ], 'table');

  let tables = ui.divV([]);
  for (let synthesizer of Object.keys(map)) {
    for (let technology of Object.keys(map[synthesizer])) {
      let tableRows = [];
      for (let [key, value] of Object.entries(map[synthesizer][technology]))
        tableRows.push({'name': value.name, 'code': key});
      tables.append(
        DG.HtmlTable.create(
          tableRows,
          (item: {name: string; code: string;}) => [item['name'], item['code']],
          [synthesizer + ' ' + technology, 'Code']
        ).root,
        ui.div([], {style: {height: '30px'}})
      );
    }
  }

  let showCodesButton = ui.button('SHOW CODES', () => ui.dialog('Codes').add(tables).show());

  let moleculeSvgDiv = ui.block([]);

  let flavor: string = (defaultInput.includes('U')) ? "RNA_both_caps" : "DNA_both_caps";
  (async () => moleculeSvgDiv.append(grok.chem.svgMol(<string> await nucleotidesToSmiles(defaultInput, flavor), 900, 300)))();

  let saveMolFileButton = ui.bigButton('SAVE MOL FILE', async () => {
    let outputSequenceObj = convertSequence(inputSequenceField.value);
    flavor = outputSequenceObj.Nucleotides.includes('U') ? "RNA_both_caps" : "DNA_both_caps";
    let smiles = (isSiRnaGcrsCode(inputSequenceField.value.replace(/\s/g, ''))) ? 
      modifiedToSmiles(inputSequenceField.value.replace(/\s/g, '')) :
      await nucleotidesToSmiles(outputSequenceObj.Nucleotides, flavor); 
    smiles = smiles.replace(/@/g, ''); // Remove StereoChemistry on the Nucleic acid chain and remove the Chiral label
    let mol = OCL.Molecule.fromSmiles(smiles);
    let result = `${mol.toMolfile()}\n`;// + '$$$$';
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', 
      ((isSiRnaGcrsCode(inputSequenceField.value.replace(/\s/g, ''))) ? inputSequenceField.value.replace(/\s/g, '') : outputSequenceObj.Nucleotides) + '.mol'
    );
    element.click();
  });

  let v = grok.shell.newView('Sequence Translator', [
    ui.tabControl({
      'MAIN': ui.div([
        appMainDescription,
        ui.panel([
          ui.div([
            ui.h1('Input sequence'),
            ui.div([
              inputSequenceField.root
            ],'input-base')
          ], 'sequenceInput'),
          semTypeOfInputSequence,
          ui.block([
            ui.h1('Output'),
            outputTableDiv
          ]),
          moleculeSvgDiv,
          ui.divH([saveMolFileButton, showCodesButton])
        ], 'sequence')
      ]),
      'AXOLABS': defineAxolabsPattern()
    })
  ]);
  v.box = true;

  $('.sequence')
    .children().css('padding','5px 0');
  $('.sequenceInput .input-base')
    .css('margin','0');
  $('.sequenceInput textarea')
    .css('resize','none')
    .css('min-height','50px')
    .css('width','100%')
    .attr("spellcheck", "false");
  $('.sequenceInput select')
    .css('width','100%');
}

function modifiedToSmiles(sequence: string) {
  const obj: {[index: string]: string} = {  
    'rA': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](O)[C@@H]1O',
    'rC': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](O)[C@@H]1O',
    'rG': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](O)[C@@H]1O',
    'rT': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))[C@H](O)[C@@H]1O',
    'rU': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](O)[C@@H]1O',
    'fU': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',
    'fA': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',
    'fC': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',
    'fG': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O',
    'mU': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',
    'mA': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
    'mC': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',
    'mG': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',
    'moe5mC': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O',
    'moeG': 'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
    'moeA': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OCCOC)[C@@H]1O',
    'moeT': 'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))[C@H](OCCOC)[C@@H]1O',
    '5mC': 'OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O',
    'ps': 'OP(=O)(O)S'
  };
  const stadardPhosphateLink = 'OP(=O)(O)O';
  const codes = Object.keys(obj);
  let i = 0;
  let smiles = '';
  let codesList = [];
  while (i < sequence.length) {
    let code = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (i < codesList.length - 1 && codesList[i+1] == 'ps') 
      smiles += obj[codesList[i]];
    else if (codesList[i] == 'ps')
      smiles += obj[codesList[i]];
    else
      smiles += obj[codesList[i]] + stadardPhosphateLink;
  }
  smiles = smiles.replace(/OO/g, 'O').replace(/SO/g, 'S');
  return smiles.slice(0, smiles.length - stadardPhosphateLink.length);
}

export async function nucleotidesToSmiles(nucleotides: string, flavor: string) {
  return await grok.functions.call('SequenceTranslator:convertFastaToSmiles', {
    'sequence_in_fasta_format': nucleotides,
    'flavor': flavor
  });
}

export function isDnaNucleotidesCode(sequence: string): boolean {return /^[ATGC]{6,}$/.test(sequence);}
export function isRnaNucleotidesCode(sequence: string): boolean {return /^[AUGC]{6,}$/.test(sequence);}
export function isAsoGapmerBioSpringCode(sequence: string): boolean {return /^[*56789ATGC]{6,}$/.test(sequence);}
export function isAsoGapmerGcrsCode(sequence: string): boolean {return /^(?=.*moe)(?=.*5mC)(?=.*ps){6,}/.test(sequence);}
export function isSiRnaBioSpringCode(sequence: string): boolean {return /^[*1-8]{6,}$/.test(sequence);}
export function isSiRnaAxolabsCode(sequence: string): boolean {return /^[fsACGUacgu]{6,}$/.test(sequence);}
export function isSiRnaGcrsCode(sequence: string): boolean {return /^[fmpsACGU]{6,}$/.test(sequence);}
export function isGcrsCode(sequence: string): boolean {return /^[fmpsACGU]{6,}$/.test(sequence);}
export function isMM12Code(sequence: string): boolean {return /^[IiJjKkLlEeFfGgHhQq]{6,}$/.test(sequence);}

function convertSequence(seq: string) {
  seq = seq.replace(/\s/g, '');
  if (seq.length < minimalValidNumberOfCharacters)
    return {
      type: smallNumberOfCharacters,
      Nucleotides: smallNumberOfCharacters,
      BioSpring: smallNumberOfCharacters,
      Axolabs: smallNumberOfCharacters,
      GCRS: smallNumberOfCharacters
    };
  if (isDnaNucleotidesCode(seq))
    return {
      type: "DNA Nucleotides Code",
      Nucleotides: seq,
      BioSpring: asoGapmersNucleotidesToBioSpring(seq),
      Axolabs: noTranslationTableAvailable,
      GCRS: asoGapmersNucleotidesToGcrs(seq)
    };
  if (isAsoGapmerBioSpringCode(seq))
    return {
      type: "ASO Gapmers / BioSpring Code",
      Nucleotides: asoGapmersBioSpringToNucleotides(seq),
      BioSpring: seq,
      Axolabs: noTranslationTableAvailable,
      GCRS: asoGapmersBioSpringToGcrs(seq)
    };
  if (isAsoGapmerGcrsCode(seq))
    return {
      type: "ASO Gapmers / GCRS Code",
      Nucleotides: asoGapmersGcrsToNucleotides(seq),
      BioSpring: asoGapmersGcrsToBioSpring(seq),
      Axolabs: noTranslationTableAvailable,
      MM12: gcrsToMM12(seq),
      GCRS: seq
    };
  if (isRnaNucleotidesCode(seq))
    return {
      type: "RNA Nucleotides Code",
      Nucleotides: seq,
      BioSpring: siRnaNucleotideToBioSpringSenseStrand(seq),
      Axolabs: siRnaNucleotideToAxolabsSenseStrand(seq),
      GCRS: siRnaNucleotidesToGcrs(seq)
    };
  if (isSiRnaBioSpringCode(seq))
    return {
      type: "siRNA / bioSpring Code",
      Nucleotides: siRnaBioSpringToNucleotides(seq),
      BioSpring: seq,
      Axolabs: siRnaBioSpringToAxolabs(seq),
      GCRS: siRnaBioSpringToGcrs(seq)
    };
  if (isSiRnaAxolabsCode(seq))
    return {
      type: "siRNA / Axolabs Code",
      Nucleotides: siRnaAxolabsToNucleotides(seq),
      BioSpring: siRnaAxolabsToBioSpring(seq),
      Axolabs: seq,
      GCRS: siRnaAxolabsToGcrs(seq)
    };
  if (isSiRnaGcrsCode(seq))
    return {
      type: "siRNA / GCRS Code",
      Nucleotides: siRnaGcrsToNucleotides(seq),
      BioSpring: siRnaGcrsToBioSpring(seq),
      Axolabs: siRnaGcrsToAxolabs(seq),
      MM12: gcrsToMM12(seq), 
      GCRS: seq
    };
  if (isGcrsCode(seq))
    return {
      type: "GCRS Code",
      Nucleotides: gcrsToNucleotides(seq),
      GCRS: seq,
      MM12: gcrsToMM12(seq)
    }
  if (isMM12Code(seq))
    return {
      type: "MM12 Code",
      Nucleotides: noTranslationTableAvailable,
      GCRS: noTranslationTableAvailable,
      MM12: seq
    };
  return {
    type: undefinedInputSequence,
    Nucleotides: undefinedInputSequence
  };
}

//name: asoGapmersNucleotidesToBioSpring
//input: string nucleotides {semType: DNA nucleotides}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersNucleotidesToBioSpring(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "5*", "A": "6*", "C": "7*", "G": "8*"};
  const objForCenter: {[index: string]: string} = {"C": "9*", "A": "A*", "T": "T*", "G": "G*"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    return (count > 4 && count < 15) ? objForCenter[x] : objForEdges[x];
  }).slice(0, 2 * count + 1);
}

//name: asoGapmersNucleotidesToGcrs
//input: string nucleotides {semType: DNA nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersNucleotidesToGcrs(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"T": "moeUnps", "A": "moeAnps", "C": "moe5mCnps", "G": "moeGnps"};
  const objForCenter: {[index: string]: string} = {"C": "5mCps", "A": "Aps", "T": "Tps", "G": "Gps"};
  return nucleotides.replace(/[ATCG]/g, function (x: string) {
    count++;
    if (count < 5) return (count == 4) ? objForEdges[x].slice(0, -3) + 'ps' : objForEdges[x];
    if (count < 15) return (count == 14) ? objForCenter[x].slice(0, -2) + 'nps' : objForCenter[x];
    return objForEdges[x];
  }).slice(0, -3);
}

//name: asoGapmersBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersBioSpringToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"*": "", "5": "T", "6": "A", "7": "C", "8": "G", "9": "C"};
  return nucleotides.replace(/[*56789]/g, function (x: string) {return obj[x];});
}

//name: asoGapmersBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersBioSpringToGcrs(nucleotides: string) {
  let count: number = -1;
  const obj: {[index: string]: string} = {
    "5*": "moeUnps", "6*": "moeAnps", "7*": "moe5mCnps", "8*": "moeGnps", "9*": "5mCps", "A*": "Aps", "T*": "Tps",
    "G*": "Gps", "C*": "Cps", "5": "moeU", "6": "moeA", "7": "moe5mC", "8": "moeG"
  };
  return nucleotides.replace(/(5\*|6\*|7\*|8\*|9\*|A\*|T\*|G\*|C\*|5|6|7|8)/g, function (x: string) {
    count++;
    return (count == 4) ? obj[x].slice(0, -3) + 'ps' : (count == 14) ? obj[x].slice(0, -2) + 'nps' : obj[x];
  });
}

//name: asoGapmersGcrsToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersGcrsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "moeT": "5", "moeA": "6", "moe5mC": "7", "moeG": "8", "moeU": "5", "5mC": "9", "nps": "*", "ps": "*", "U": "T"
  };
  return nucleotides.replace(/(moeT|moeA|moe5mC|moeG|moeU|5mC|nps|ps|U)/g, function (x: string) {return obj[x];});
}

//name: asoGapmersGcrsToNucleotides
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"moe": "", "5m": "", "n": "", "ps": "", "U": "T"};
  return nucleotides.replace(/(moe|5m|n|ps|U)/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaBioSpringToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"1": "U", "2": "A", "3": "C", "4": "G", "5": "U", "6": "A", "7": "C", "8": "G", "*": ""};
  return nucleotides.replace(/[12345678*]/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToAxolabs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: Axolabs / siRNA}
export function siRnaBioSpringToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {"1": "Uf", "2": "Af", "3": "Cf", "4": "Gf", "5": "u", "6": "a", "7": "c", "8": "g", "*": "s"};
  return nucleotides.replace(/[12345678*]/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: GCRS}
export function siRnaBioSpringToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {"1": "fU", "2": "fA", "3": "fC", "4": "fG", "5": "mU", "6": "mA", "7": "mC", "8": "mG", "*": "ps"};
  return nucleotides.replace(/[12345678*]/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToGcrs
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: GCRS}
export function siRnaAxolabsToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "fU", "Af": "fA", "Cf": "fC", "Gf": "fG", "u": "mU", "a": "mA", "c": "mC", "g": "mG", "s": "ps"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToBioSpring
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: BioSpring / siRNA}
export function siRnaAxolabsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "1", "Af": "2", "Cf": "3", "Gf": "4", "u": "5", "a": "6", "c": "7", "g": "8", "s": "*"
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToNucleotides
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaAxolabsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "Uf": "U", "Af": "A", "Cf": "C", "Gf": "G", "u": "U", "a": "A", "c": "C", "g": "G", "s": ""
  };
  return nucleotides.replace(/(Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToNucleotides
//input: string nucleotides {semType: GCRS}
//output: string result {semType: RNA nucleotides}
export function siRnaGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "fU": "U", "fA": "A", "fC": "C", "fG": "G", "mU": "U", "mA": "A", "mC": "C", "mG": "G", "ps": ""
  };
  return nucleotides.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToBioSpring
//input: string nucleotides {semType: GCRS}
//output: string result {semType: BioSpring / siRNA}
export function siRnaGcrsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "fU": "1", "fA": "2", "fC": "3", "fG": "4", "mU": "5", "mA": "6", "mC": "7", "mG": "8", "ps": "*"
  };
  return nucleotides.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToAxolabs
//input: string nucleotides {semType: GCRS}
//output: string result {semType: Axolabs / siRNA}
export function siRnaGcrsToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "fU": "Uf", "fA": "Af", "fC": "Cf", "fG": "Gf", "mU": "u", "mA": "a", "mC": "c", "mG": "g", "ps": "s"
  };
  return nucleotides.replace(/(fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaNucleotideToBioSpringSenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: BioSpring / siRNA}
export function siRnaNucleotideToBioSpringSenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {"A": "6*", "U": "5*", "G": "8*", "C": "7*"};
  const objForRightEdge: {[index: string]: string} = {"A": "*6", "U": "*5", "G": "*8", "C": "*7"};
  const objForOddIndices: {[index: string]: string} = {"A": "6", "U": "5", "G": "8", "C": "7"};
  const objForEvenIndices: {[index: string]: string} = {"A": "2", "U": "1", "G": "4", "C": "3"};
  return nucleotides.replace(/[AUGC]/g, function (x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count > nucleotides.length - 3) return objForRightEdge[x];
    return (count % 2 == 0) ? objForEvenIndices[x] : objForOddIndices[x];
  });
}

//name: siRnaNucleotidesToGcrs
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: GCRS}
export function siRnaNucleotidesToGcrs(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {"A": "mAps", "U": "mUps", "G": "mGps", "C": "mCps"};
  const objForRightEdge: {[index: string]: string} = {"A": "psmA", "U": "psmU", "G": "psmG", "C": "psmC"};
  const objForEvenIndices: {[index: string]: string} = {"A": "fA", "U": "fU", "G": "fG", "C": "fC"};
  const objForOddIndices: {[index: string]: string} = {"A": "mA", "U": "mU", "G": "mG", "C": "mC"};
  return nucleotides.replace(/[AUGC]/g, function (x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count > nucleotides.length - 3) return objForRightEdge[x];
    return (count % 2 == 0) ? objForEvenIndices[x] : objForOddIndices[x];
  });
}

//name: siRnaNucleotideToAxolabsSenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function siRnaNucleotideToAxolabsSenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {"A": "as", "U": "us", "G": "gs", "C": "cs"};
  const objForSomeIndices: {[index: string]: string} = {"A": "Af", "U": "Uf", "G": "Gf", "C": "Cf"};
  const obj: {[index: string]: string} = {"A": "a", "U": "u", "G": "g", "C": "c"};
  return nucleotides.replace(/[AUGC]/g, function (x: string) {
    count++;
    if (count < 2) return objForLeftEdge[x];
    if (count == 6 || (count > 7 && count < 11)) return objForSomeIndices[x]
    if (count == nucleotides.length - 1) return 'a';
    return obj[x];
  });
}

//name: siRnaNucleotideToAxolabsAntisenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: Axolabs}
export function siRnaNucleotideToAxolabsAntisenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForSmallLinkages: {[index: string]: string} = {"A": "as", "U": "us", "G": "gs", "C": "cs"};
  const objForBigLinkages: {[index: string]: string} = {"A": "Afs", "U": "Ufs", "G": "Gfs", "C": "Cfs"};
  const objForSomeIndices: {[index: string]: string} = {"A": "Af", "U": "Uf", "G": "Gf", "C": "Cf"};
  const obj: {[index: string]: string} = {"A": "a", "U": "u", "G": "g", "C": "c"};
  return nucleotides.replace(/[AUGC]/g, function (x: string) {
    count++;
    if (count > 19 && count < 22) return objForSmallLinkages[x];
    if (count == 0) return 'us';
    if (count == 1) return objForBigLinkages[x];
    return (count == 5 || count == 7 || count == 8 || count == 13 || count == 15) ? objForSomeIndices[x] : obj[x];
  });
}

//name: gcrsToNucleotides
//input: string nucleotides {semType: GCRS}
//output: string result {semType: RNA nucleotides}
export function gcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "mAps": "A", "mUps": "U", "mGps": "G", "mCps": "C", "fAps": "A", "fUps": "U", "fGps": "G", "fCps": "C",
    "fU": "U", "fA": "A", "fC": "C", "fG": "G", "mU": "U", "mA": "A", "mC": "C", "mG": "G"
  };
  return nucleotides.replace(/(mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g, function (x: string) {return obj[x];});
}

//name: gcrsToOP100
//input: string nucleotides {semType: GCRS}
//output: string result {semType: OP100}
export function gcrsToOP100(nucleotides: string) {
  let count: number = -1;
  const objForEvenIndicesAtLeftEdge: {[index: string]: string} = {
    "mAps": "a", "mUps": "u", "mGps": "g", "mCps": "c", "fAps": "a", "fUps": "u", "fGps": "g", "fCps": "c"
  };
  const objForOddIndicesAtLeftEdge: {[index: string]: string} = {
    "mAps": "a*", "mUps": "u*", "mGps": "g*", "mCps": "c*", "fAps": "a*", "fUps": "u*", "fGps": "g*", "fCps": "c*"
  };
  const objForOddIndicesAtRightEdge: {[index: string]: string} = {
    "mAps": "a", "mUps": "u", "mGps": "g", "mCps": "c", "fAps": "a", "fUps": "u", "fGps": "g", "fCps": "c"
  };
  const objForEvenIndicesAtCenter: {[index: string]: string} = {
    "fU": "u*", "fA": "a*", "fC": "c*", "fG": "g*", "mU": "u*", "mA": "a*", "mC": "c*", "mG": "g*"
  };
  const objForOddIndicesAtCenter: {[index: string]: string} = {
    "fU": "u", "fA": "a", "fC": "c", "fG": "g", "mU": "u", "mA": "a", "mC": "c", "mG": "g"
  };
  return nucleotides.replace(/(mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g, function (x: string) {
    count++;
    if (count < 3) return (count % 2 == 0) ? objForEvenIndicesAtLeftEdge[x] : objForOddIndicesAtLeftEdge[x];
    if (count == 19) return objForOddIndicesAtRightEdge[x];
    return (count % 2 == 1) ? objForEvenIndicesAtCenter[x] : objForOddIndicesAtCenter[x];
  });
}

//name: gcrsToMM12
//input: string nucleotides {semType: GCRS}
//output: string result {semType: MM12}
export function gcrsToMM12(nucleotides: string) {
  const obj: {[index: string]: string} = {
    "mAps": "e", "mUps": "h", "mGps": "g", "mCps": "f", "fAps": "i", "fUps": "l", "fGps": "k", "fCps": "j", "fU": "L",
    "fA": "I", "fC": "J", "fG": "K", "mU": "H", "mA": "E", "mC": "F", "mG": "G"
  };
  return nucleotides.replace(/(mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g, function (x: string) {return obj[x]});
}