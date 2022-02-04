/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import $ from "cash-dom";
import {defineAxolabsPattern} from "./defineAxolabsPattern";
import {map, stadardPhosphateLinkSmiles, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS} from "./map";

export let _package = new DG.Package();

const defaultInput = "AGGTCCTCTTGACTTAGGCC";
const undefinedInputSequence = "Type of input sequence is undefined";
const noTranslationTableAvailable = "No translation table available";
const sequenceWasCopied = 'Copied';
const tooltipSequence = 'Copy sequence';

function getAllCodesOfSynthesizer(synthesizer: string) {
  let codes: string[] = [];
  for (let technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes.concat(Object.keys(MODIFICATIONS));
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string): string[] {
  let synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    //TODO: get first non-dropdown code when there are two modifications
    let start = 0;
    for (let i = 0; i < sequence.length; i++)
      if (sequence[i] == ')') {
        start = i + 1;
        break;
      }
    if (codes.some((s: string) => s == sequence.slice(start, start + s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

function getListOfPossibleTechnologiesByFirstMatchedCode(sequence: string, synthesizer: string): string[] {
  let technologies: string[] = [];
  Object.keys(map[synthesizer]).forEach((technology: string) => {
    const codes = Object.keys(map[synthesizer][technology]).concat(Object.keys(MODIFICATIONS));
    if (codes.some((s) => s == sequence.slice(0, s.length)))
      technologies.push(technology);
  });
  return technologies;
}

function isValidSequence(sequence: string) {
  let possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence);
  if (possibleSynthesizers.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedType: null };

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
  const expectedSynthesizer = possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)];
  if (indexOfFirstNotValidCharacter != -1)
    return {
      indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
      expectedType: expectedSynthesizer
    };

  let possibleTechnologies = getListOfPossibleTechnologiesByFirstMatchedCode(sequence, expectedSynthesizer);
  if (possibleTechnologies.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedRepresentation: null };

  outputIndices = Array(possibleTechnologies.length).fill(0);

  possibleTechnologies.forEach((technology, technologyIndex) => {
    let codes = Object.keys(map[expectedSynthesizer][technology]);
    while (outputIndices[technologyIndex] < sequence.length) {

      let matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[technologyIndex], outputIndices[technologyIndex] + c.length));

      if (matchedCode == null)
        break;

      if (  // for mistake pattern 'rAA'
        outputIndices[technologyIndex] > 1 &&
        nucleotides.includes(sequence[outputIndices[technologyIndex]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[technologyIndex] - 2])
      ) break;

      if (  // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[technologyIndex] + 1]) &&
        nucleotides.includes(sequence[outputIndices[technologyIndex]])
      ) {
        outputIndices[technologyIndex]++;
        break;
      }

      outputIndices[technologyIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedTechnology = Math.max.apply(Math, outputIndices);
  const expectedTechnology = possibleTechnologies[outputIndices.indexOf(indexOfExpectedTechnology)];

  return {
    indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
    expectedType: expectedSynthesizer + ' ' + expectedTechnology
  };
}

function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a: string, b: string) { return b.length - a.length; });
}

function getObjectWithCodesAndSmiles() {
  let obj: {[code: string]: string} = {};
  for (let synthesizer of Object.keys(map))
    for (let technology of Object.keys(map[synthesizer]))
      for (let code of Object.keys(map[synthesizer][technology]))
        obj[code] = map[synthesizer][technology][code].SMILES;
  return obj;
}

function sequenceToSmiles(sequence: string) {
  const obj = getObjectWithCodesAndSmiles();
  let codes = sortByStringLengthInDescendingOrder(Object.keys(obj));
  let i = 0, smiles = '', codesList = [];
  const links = ['s', 'ps', '*'];
  const includesStandardLinkAlready = ["e", "h", "g", "f", "i", "l", "k", "j"];
  const dropdowns = Object.keys(MODIFICATIONS);
  codes = codes.concat(dropdowns);
  while (i < sequence.length) {
    let code = codes.find((s: string) => s == sequence.slice(i, i + s.length))!;
    i += code.length;
    codesList.push(code);
  }
  for (let i = 0; i < codesList.length; i++) {
    if (dropdowns.includes(codesList[i])) {
      smiles += (i >= codesList.length / 2) ?
        MODIFICATIONS[codesList[i]].right :
        MODIFICATIONS[codesList[i]].left;
    } else {
      if (links.includes(codesList[i]) && i > 1 && !includesStandardLinkAlready.includes(codesList[i - 1]))
        smiles = smiles.slice(0, smiles.length - stadardPhosphateLinkSmiles.length + 1);
      else if (links.includes(codesList[i]) ||
        includesStandardLinkAlready.includes(codesList[i]) ||
        (i < codesList.length - 1 && (links.includes(codesList[i + 1]) || dropdowns.includes(codesList[i + 1])))
      )
        smiles += obj[codesList[i]];
      else
        smiles += obj[codesList[i]] + stadardPhosphateLinkSmiles;
    }
  }
  smiles = smiles.replace(/OO/g, 'O');
  return (
    (
      links.includes(codesList[codesList.length - 1]) &&
      codesList.length > 1 &&
      !includesStandardLinkAlready.includes(codesList[codesList.length - 2])
    ) ||
    dropdowns.includes(codesList[codesList.length - 1]) ||
    includesStandardLinkAlready.includes(codesList[codesList.length - 1])
  ) ?
    smiles :
    smiles.slice(0, smiles.length - stadardPhosphateLinkSmiles.length + 1);
}

//name: Sequence Translator
//tags: app
export function sequenceTranslator() {

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  function updateTableAndMolecule(sequence: string) {
    moleculeSvgDiv.innerHTML = "";
    outputTableDiv.innerHTML = "";
    let pi = DG.TaskBarProgressIndicator.create('Rendering table and molecule...');
    try {
      let outputSequenceObj = convertSequence(sequence);
      let tableRows = [];
      for (let key of Object.keys(outputSequenceObj).slice(1)) {
        tableRows.push({
          'key': key,
          'value': ("indexOfFirstNotValidCharacter" in outputSequenceObj) ?
            ui.divH([
              ui.divText(sequence.slice(0, JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).indexOfFirstNotValidCharacter), {style: {color: "grey"}}),
              ui.tooltip.bind(
                ui.divText(sequence.slice(JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).indexOfFirstNotValidCharacter), {style: {color: "red"}}),
                "Expected format: " + JSON.parse(outputSequenceObj.indexOfFirstNotValidCharacter!).expectedType + ". Press 'SHOW CODES' button to see tables with valid codes"
              )
            ]) : //@ts-ignore
            ui.link(outputSequenceObj[key], () => navigator.clipboard.writeText(outputSequenceObj[key]).then(() => grok.shell.info(sequenceWasCopied)), tooltipSequence, '')
        })
      }
      outputTableDiv.append(
        ui.div([DG.HtmlTable.create(tableRows, (item: { key: string; value: string; }) => [item.key, item.value], ['Code', 'Sequence']).root], 'table')
      );
      semTypeOfInputSequence.textContent = 'Detected input type: ' + outputSequenceObj.type;

      let width = $(window).width();
      const canvas = ui.canvas(width, Math.round(width / 2));
      let smiles = sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''));
      // @ts-ignore
      OCL.StructureView.drawMolecule(canvas, OCL.Molecule.fromSmiles(smiles), { suppressChiralText: true });
      if (outputSequenceObj.type != undefinedInputSequence)
        moleculeSvgDiv.append(canvas);
    } finally {
      pi.close();
    }
  }

  let semTypeOfInputSequence = ui.divText('');
  let moleculeSvgDiv = ui.block([]);
  let outputTableDiv = ui.div([], 'table');
  let inputSequenceField = ui.textInput("", defaultInput, (sequence: string) => updateTableAndMolecule(sequence));
  updateTableAndMolecule(defaultInput);

  let tablesWithCodes = ui.divV([
    DG.HtmlTable.create(Object.keys(MODIFICATIONS), (item: string) => [item], ['Overhang modification']).root,
    ui.div([], {style: {height: '30px'}})
  ]);
  for (let synthesizer of Object.keys(map)) {
    for (let technology of Object.keys(map[synthesizer])) {
      let tableRows = [];
      for (let [key, value] of Object.entries(map[synthesizer][technology]))
        tableRows.push({'name': value.name, 'code': key});
      tablesWithCodes.append(
        DG.HtmlTable.create(
          tableRows,
          (item: {name: string; code: string;}) => [item['name'], item['code']],
          [synthesizer + ' ' + technology, 'Code']
        ).root,
        ui.div([], {style: {height: '30px'}})
      );
    }
  }
  let showCodesButton = ui.button('SHOW CODES', () => ui.dialog('Codes').add(tablesWithCodes).show());
  let copySmiles = ui.button(
    'COPY SMILES',
    () => navigator.clipboard.writeText(sequenceToSmiles(inputSequenceField.value.replace(/\s/g, '')))
      .then(() => grok.shell.info(sequenceWasCopied))
  );
  let saveMolFileButton = ui.bigButton('SAVE MOL FILE', () => {
    let smiles = sequenceToSmiles(inputSequenceField.value.replace(/\s/g, ''));
    let result = `${OCL.Molecule.fromSmiles(smiles).toMolfile()}\n`;
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', inputSequenceField.value.replace(/\s/g, '') + '.mol');
    element.click();
  });

  const appMainDescription = ui.info([
      ui.divText('\n How to convert one sequence:',{style:{'font-weight':'bolder'}}),
      ui.divText("Paste sequence into the text field below"),
      ui.divText('\n How to convert many sequences:',{style:{'font-weight':'bolder'}}),
      ui.divText("1. Drag & drop an Excel or CSV file with sequences into Datagrok. The platform will automatically detect columns with sequences"),
      ui.divText('2. Right-click on the column header, then see the \'Convert\' menu'),
      ui.divText("This will add the result column to the right of the table"),
    ], 'Convert oligonucleotide sequences between Nucleotides, BioSpring, Axolabs, Mermade 12 and GCRS representations.'
  );

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
          ui.divH([saveMolFileButton, showCodesButton, copySmiles])
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

function convertSequence(text: string) {
  text = text.replace(/\s/g, '');
  let seq = text;
  let output = isValidSequence(seq);
  if (output.indexOfFirstNotValidCharacter != -1)
    return {
      indexOfFirstNotValidCharacter: JSON.stringify(output),
      Error: undefinedInputSequence
    };
  if (output.expectedType == SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' + TECHNOLOGIES.DNA)
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' +  TECHNOLOGIES.DNA,
      Nucleotides: seq,
      BioSpring: asoGapmersNucleotidesToBioSpring(seq),
      GCRS: asoGapmersNucleotidesToGcrs(seq)
    };
  if (output.expectedType == SYNTHESIZERS.BIOSPRING + ' ' +  TECHNOLOGIES.ASO_GAPMERS)
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' +  TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersBioSpringToNucleotides(seq),
      BioSpring: seq,
      GCRS: asoGapmersBioSpringToGcrs(seq)
    };
  if (output.expectedType == SYNTHESIZERS.GCRS + ' ' +  TECHNOLOGIES.ASO_GAPMERS)
    return {
      type: SYNTHESIZERS.GCRS + ' ' +  TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersGcrsToNucleotides(seq),
      BioSpring: asoGapmersGcrsToBioSpring(seq),
      Mermade12: gcrsToMermade12(seq),
      GCRS: seq
    };
  if (output.expectedType == SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' +  TECHNOLOGIES.RNA)
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' +  TECHNOLOGIES.RNA,
      Nucleotides: seq,
      BioSpring: siRnaNucleotideToBioSpringSenseStrand(seq),
      Axolabs: siRnaNucleotideToAxolabsSenseStrand(seq),
      GCRS: siRnaNucleotidesToGcrs(seq)
    };
  if (output.expectedType == SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.SI_RNA)
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaBioSpringToNucleotides(seq),
      BioSpring: seq,
      Axolabs: siRnaBioSpringToAxolabs(seq),
      GCRS: siRnaBioSpringToGcrs(seq)
    };
  if (output.expectedType == SYNTHESIZERS.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA)
    return {
      type: SYNTHESIZERS.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaAxolabsToNucleotides(seq),
      BioSpring: siRnaAxolabsToBioSpring(seq),
      Axolabs: seq,
      GCRS: siRnaAxolabsToGcrs(seq)
    };
  if (output.expectedType == SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.SI_RNA)
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaGcrsToNucleotides(seq),
      BioSpring: siRnaGcrsToBioSpring(seq),
      Axolabs: siRnaGcrsToAxolabs(seq),
      MM12: gcrsToMermade12(seq),
      GCRS: seq
    };
  if (output.expectedType == SYNTHESIZERS.GCRS)
    return {
      type: SYNTHESIZERS.GCRS,
      Nucleotides: gcrsToNucleotides(seq),
      GCRS: seq,
      Mermade12: gcrsToMermade12(seq)
    }
  if (output.expectedType == SYNTHESIZERS.MERMADE_12)
    return {
      type: SYNTHESIZERS.MERMADE_12,
      Nucleotides: noTranslationTableAvailable,
      GCRS: noTranslationTableAvailable,
      Mermade12: seq
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
  const objForEdges: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "T": "5*", "A": "6*", "C": "7*", "G": "8*"};
  const objForCenter: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "T": "T*", "A": "A*", "C": "9*", "G": "G*"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|T|C|G)/g, function (x: string) {
    count++;
    return (count > 4 && count < 15) ? objForCenter[x] : objForEdges[x];
  }).slice(0, (nucleotides.endsWith("(invabasic)") || nucleotides.endsWith("(GalNAc-2-JNJ)")) ? nucleotides.length : 2 * count + 1);
}

//name: asoGapmersNucleotidesToGcrs
//input: string nucleotides {semType: DNA nucleotides}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersNucleotidesToGcrs(nucleotides: string) {
  let count: number = -1;
  const objForEdges: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "T": "moeUnps", "A": "moeAnps", "C": "moe5mCnps", "G": "moeGnps"};
  const objForCenter: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "C": "5mCps", "A": "Aps", "T": "Tps", "G": "Gps"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|T|C|G)/g, function (x: string) {
    count++;
    if (count < 5) return (count == 4) ? objForEdges[x].slice(0, -3) + 'ps' : objForEdges[x];
    if (count < 15) return (count == 14) ? objForCenter[x].slice(0, -2) + 'nps' : objForCenter[x];
    return objForEdges[x];
  }).slice(0, (nucleotides.endsWith("(invabasic)") || nucleotides.endsWith("(GalNAc-2-JNJ)")) ? nucleotides.length : -3);
}

//name: asoGapmersBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersBioSpringToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "*": "", "5": "T", "6": "A", "7": "C", "8": "G", "9": "C"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|5|6|7|8|9)/g, function (x: string) {return obj[x];});
}

//name: asoGapmersBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / Gapmers}
//output: string result {semType: GCRS / Gapmers}
export function asoGapmersBioSpringToGcrs(nucleotides: string) {
  let count: number = -1;
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "5*": "moeUnps", "6*": "moeAnps", "7*": "moe5mCnps", "8*": "moeGnps", "9*": "5mCps", "A*": "Aps", "T*": "Tps",
    "G*": "Gps", "C*": "Cps", "5": "moeU", "6": "moeA", "7": "moe5mC", "8": "moeG"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|5\*|6\*|7\*|8\*|9\*|A\*|T\*|G\*|C\*|5|6|7|8)/g, function (x: string) {
    count++;
    return (count == 4) ? obj[x].slice(0, -3) + 'ps' : (count == 14) ? obj[x].slice(0, -2) + 'nps' : obj[x];
  });
}

//name: asoGapmersGcrsToBioSpring
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: BioSpring / Gapmers}
export function asoGapmersGcrsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "moeT": "5", "moeA": "6", "moe5mC": "7", "moeG": "8", "moeU": "5", "5mC": "9", "nps": "*", "ps": "*", "U": "T"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|moeT|moeA|moe5mC|moeG|moeU|5mC|nps|ps|U)/g, function (x: string) {return obj[x];});
}

//name: asoGapmersGcrsToNucleotides
//input: string nucleotides {semType: GCRS / Gapmers}
//output: string result {semType: DNA nucleotides}
export function asoGapmersGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "moe": "", "5m": "", "n": "", "ps": "", "U": "T"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|moe|5m|n|ps|U)/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToNucleotides
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaBioSpringToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "1": "U", "2": "A", "3": "C", "4": "G", "5": "U", "6": "A", "7": "C", "8": "G", "*": ""};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToAxolabs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: Axolabs / siRNA}
export function siRnaBioSpringToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "1": "Uf", "2": "Af", "3": "Cf", "4": "Gf", "5": "u", "6": "a", "7": "c", "8": "g", "*": "s"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g, function (x: string) {return obj[x];});
}

//name: siRnaBioSpringToGcrs
//input: string nucleotides {semType: BioSpring / siRNA}
//output: string result {semType: GCRS}
export function siRnaBioSpringToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "1": "fU", "2": "fA", "3": "fC", "4": "fG", "5": "mU", "6": "mA", "7": "mC", "8": "mG", "*": "ps"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|1|2|3|4|5|6|7|8|\*)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToGcrs
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: GCRS}
export function siRnaAxolabsToGcrs(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "Uf": "fU", "Af": "fA", "Cf": "fC", "Gf": "fG", "u": "mU", "a": "mA", "c": "mC", "g": "mG", "s": "ps"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToBioSpring
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: BioSpring / siRNA}
export function siRnaAxolabsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "Uf": "1", "Af": "2", "Cf": "3", "Gf": "4", "u": "5", "a": "6", "c": "7", "g": "8", "s": "*"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaAxolabsToNucleotides
//input: string nucleotides {semType: Axolabs / siRNA}
//output: string result {semType: RNA nucleotides}
export function siRnaAxolabsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "Uf": "U", "Af": "A", "Cf": "C", "Gf": "G", "u": "U", "a": "A", "c": "C", "g": "G", "s": ""
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|Uf|Af|Cf|Gf|u|a|c|g|s)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToNucleotides
//input: string nucleotides {semType: GCRS}
//output: string result {semType: RNA nucleotides}
export function siRnaGcrsToNucleotides(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "fU": "U", "fA": "A", "fC": "C", "fG": "G", "mU": "U", "mA": "A", "mC": "C", "mG": "G", "ps": ""
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToBioSpring
//input: string nucleotides {semType: GCRS}
//output: string result {semType: BioSpring / siRNA}
export function siRnaGcrsToBioSpring(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "fU": "1", "fA": "2", "fC": "3", "fG": "4", "mU": "5", "mA": "6", "mC": "7", "mG": "8", "ps": "*"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaGcrsToAxolabs
//input: string nucleotides {semType: GCRS}
//output: string result {semType: Axolabs / siRNA}
export function siRnaGcrsToAxolabs(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "fU": "Uf", "fA": "Af", "fC": "Cf", "fG": "Gf", "mU": "u", "mA": "a", "mC": "c", "mG": "g", "ps": "s"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|fU|fA|fC|fG|mU|mA|mC|mG|ps)/g, function (x: string) {return obj[x];});
}

//name: siRnaNucleotideToBioSpringSenseStrand
//input: string nucleotides {semType: RNA nucleotides}
//output: string result {semType: BioSpring / siRNA}
export function siRnaNucleotideToBioSpringSenseStrand(nucleotides: string) {
  let count: number = -1;
  const objForLeftEdge: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "6*", "U": "5*", "G": "8*", "C": "7*"};
  const objForRightEdge: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "*6", "U": "*5", "G": "*8", "C": "*7"};
  const objForOddIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "6", "U": "5", "G": "8", "C": "7"};
  const objForEvenIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "2", "U": "1", "G": "4", "C": "3"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function (x: string) {
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
  const objForLeftEdge: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "mAps", "U": "mUps", "G": "mGps", "C": "mCps"};
  const objForRightEdge: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "psmA", "U": "psmU", "G": "psmG", "C": "psmC"};
  const objForEvenIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "fA", "U": "fU", "G": "fG", "C": "fC"};
  const objForOddIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "mA", "U": "mU", "G": "mG", "C": "mC"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function (x: string) {
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
  const objForLeftEdge: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "as", "U": "us", "G": "gs", "C": "cs"};
  const objForSomeIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "Af", "U": "Uf", "G": "Gf", "C": "Cf"};
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "a", "U": "u", "G": "g", "C": "c"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function (x: string) {
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
  const objForSmallLinkages: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "as", "U": "us", "G": "gs", "C": "cs"};
  const objForBigLinkages: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "Afs", "U": "Ufs", "G": "Gfs", "C": "Cfs"};
  const objForSomeIndices: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "Af", "U": "Uf", "G": "Gf", "C": "Cf"};
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)", "A": "a", "U": "u", "G": "g", "C": "c"};
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C)/g, function (x: string) {
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
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "mAps": "A", "mUps": "U", "mGps": "G", "mCps": "C", "fAps": "A", "fUps": "U", "fGps": "G", "fCps": "C",
    "fU": "U", "fA": "A", "fC": "C", "fG": "G", "mU": "U", "mA": "A", "mC": "C", "mG": "G"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g, function (x: string) {return obj[x];});
}

//name: gcrsToMermade12
//input: string nucleotides {semType: GCRS}
//output: string result {semType: Mermade 12 / siRNA}
export function gcrsToMermade12(nucleotides: string) {
  const obj: {[index: string]: string} = {"(invabasic)": "(invabasic)", "(GalNAc-2-JNJ)": "(GalNAc-2-JNJ)",
    "mAps": "e", "mUps": "h", "mGps": "g", "mCps": "f", "fAps": "i", "fUps": "l", "fGps": "k", "fCps": "j", "fU": "L",
    "fA": "I", "fC": "J", "fG": "K", "mU": "H", "mA": "E", "mC": "F", "mG": "G"
  };
  return nucleotides.replace(/(\(invabasic\)|\(GalNAc-2-JNJ\)|mAps|mUps|mGps|mCps|fAps|fUps|fGps|fCps|fU|fA|fC|fG|mU|mA|mC|mG)/g, function (x: string) {return obj[x]});
}
