
import {map, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS} from './map';
import {asoGapmersNucleotidesToBioSpring, asoGapmersNucleotidesToGcrs,
  asoGapmersBioSpringToNucleotides, asoGapmersBioSpringToGcrs, asoGapmersGcrsToNucleotides,
  asoGapmersGcrsToBioSpring, gcrsToMermade12, siRnaNucleotideToBioSpringSenseStrand,
  siRnaNucleotideToAxolabsSenseStrand, siRnaNucleotidesToGcrs, siRnaBioSpringToNucleotides,
  siRnaBioSpringToAxolabs, siRnaBioSpringToGcrs, siRnaAxolabsToNucleotides,
  siRnaAxolabsToBioSpring, siRnaAxolabsToGcrs, siRnaGcrsToNucleotides,
  siRnaGcrsToBioSpring, siRnaGcrsToAxolabs, gcrsToNucleotides, gcrsToLcms} from './converters';

const noTranslationTableAvailable = 'No translation table available';
export const undefinedInputSequence = 'Type of input sequence is undefined';

export function getFormat(sequence: string): string | null {
  const possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence);

  if (possibleSynthesizers.length == 0)
    return null;

  let outputIndex = 0;

  const firstUniqueCharacters = ['r', 'd'];
  const nucleotides = ['A', 'U', 'T', 'C', 'G'];

  possibleSynthesizers.forEach((synthesizer) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    while (outputIndex < sequence.length) {
      const matchedCode = codes.find((c) => c == sequence.slice(outputIndex, outputIndex + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndex > 1 &&
        nucleotides.includes(sequence[outputIndex]) &&
        firstUniqueCharacters.includes(sequence[outputIndex - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndex + 1]) &&
        nucleotides.includes(sequence[outputIndex])
      ) {
        outputIndex++;
        break;
      }

      outputIndex += matchedCode.length;
    }
  });

  const indexOfFirstNotValidChar = (outputIndex == sequence.length) ? -1 : outputIndex;
  if (indexOfFirstNotValidChar != -1)
    return possibleSynthesizers[0];

  const possibleTechnologies = getListOfPossibleTechnologiesByFirstMatchedCode(sequence, possibleSynthesizers[0]);

  if (possibleTechnologies.length == 0)
    return null;

  outputIndex = 0;

  possibleTechnologies.forEach((technology: string) => {
    const codes = Object.keys(map[possibleSynthesizers[0]][technology]);
    while (outputIndex < sequence.length) {
      const matchedCode = codes.find((c) => c == sequence.slice(outputIndex, outputIndex + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndex > 1 &&
        nucleotides.includes(sequence[outputIndex]) &&
        firstUniqueCharacters.includes(sequence[outputIndex - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndex + 1]) &&
        nucleotides.includes(sequence[outputIndex])
      ) {
        outputIndex++;
        break;
      }

      outputIndex += matchedCode.length;
    }
  });

  return possibleSynthesizers[0];
}

export function isValidSequence(sequence: string, format: string | null): {
  indexOfFirstNotValidChar: number,
  synthesizer: string[] | null,
  technology: string[] | null
} {
  const possibleSynthesizers = format == null ?
    getListOfPossibleSynthesizersByFirstMatchedCode(sequence) :
    [format];

  // if (possibleSynthesizers.length > 1) {
  //   const synthesizer = ui.choiceInput('Choose synthesizer from list: ', possibleSynthesizers[0],
  //  possibleSynthesizers);
  //   ui.dialog('Choose Synthesizer')
  //     .add(ui.panel([synthesizer.root], {style: {fontWeight: 'bold'}}))
  //     .onOK(() => possibleSynthesizers = [synthesizer.value])
  //     .onCancel(() => {
  //       possibleSynthesizers = [possibleSynthesizers[0]];
  //       grok.shell.warning('Input sequence is expected to be in format ' + possibleSynthesizers[0]);
  //     })
  //     .show();
  // } else if (possibleSynthesizers.length == 0)
  if (possibleSynthesizers.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: null, technology: null};

  let outputIndex = 0;

  const firstUniqueCharacters = ['r', 'd'];
  const nucleotides = ['A', 'U', 'T', 'C', 'G'];

  possibleSynthesizers.forEach((synthesizer) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    while (outputIndex < sequence.length) {
      const matchedCode = codes.find((c) => c == sequence.slice(outputIndex, outputIndex + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndex > 1 &&
        nucleotides.includes(sequence[outputIndex]) &&
        firstUniqueCharacters.includes(sequence[outputIndex - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndex + 1]) &&
        nucleotides.includes(sequence[outputIndex])
      ) {
        outputIndex++;
        break;
      }

      outputIndex += matchedCode.length;
    }
  });

  const indexOfFirstNotValidChar = (outputIndex == sequence.length) ? -1 : outputIndex;
  if (indexOfFirstNotValidChar != -1) {
    return {
      indexOfFirstNotValidChar: indexOfFirstNotValidChar,
      synthesizer: possibleSynthesizers,
      technology: null,
    };
  }

  const possibleTechnologies = getListOfPossibleTechnologiesByFirstMatchedCode(sequence, possibleSynthesizers[0]);

  // if (possibleTechnologies.length > 1) {
  //   const technology = ui.choiceInput('Choose technology from list: ', possibleTechnologies[0],
  // possibleTechnologies);
  //   ui.dialog('Choose Technology')
  //     .add(ui.panel([technology.root], {style: {fontWeight: 'bold'}}))
  //     .onOK(() => possibleTechnologies = [technology.value])
  //     .onCancel(() => {
  //       possibleTechnologies = [possibleTechnologies[0]];
  //       grok.shell.warning('Input sequence is expected to be in format ' + possibleTechnologies[0]);
  //     })
  //     .show();
  // } else if (possibleTechnologies.length == 0)
  if (possibleTechnologies.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: [possibleSynthesizers[3]], technology: null};

  outputIndex = 0;

  // possibleTechnologies.forEach((technology: string) => {
  //   const codes = Object.keys(map[possibleSynthesizers[0]][technology]);
  //   while (outputIndex < sequence.length) {
  //     const matchedCode = codes.find((c) => c == sequence.slice(outputIndex, outputIndex + c.length));

  //     if (matchedCode == null)
  //       break;

  //     if ( // for mistake pattern 'rAA'
  //       outputIndex > 1 &&
  //       nucleotides.includes(sequence[outputIndex]) &&
  //       firstUniqueCharacters.includes(sequence[outputIndex - 2])
  //     ) break;

  //     if ( // for mistake pattern 'ArA'
  //       firstUniqueCharacters.includes(sequence[outputIndex + 1]) &&
  //       nucleotides.includes(sequence[outputIndex])
  //     ) {
  //       outputIndex++;
  //       break;
  //     }

  //     outputIndex += matchedCode.length;
  //   }
  // });

  return {
    indexOfFirstNotValidChar: indexOfFirstNotValidChar,
    synthesizer: possibleSynthesizers,
    technology: [possibleTechnologies[outputIndex]],
  };
}

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (const technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes.concat(Object.keys(MODIFICATIONS)).concat(',');
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string): string[] {
  const synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    //TODO: get first non-dropdown code when there are two modifications
    let start = 0;
    for (let i = 0; i < sequence.length; i++) {
      if (sequence[i] == ')' && i != sequence.length - 1) {
        start = i + 1;
        break;
      }
    }
    if (codes.some((s: string) => s == sequence.slice(start, start + s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

function getListOfPossibleTechnologiesByFirstMatchedCode(sequence: string, synthesizer: string): string[] {
  const technologies: string[] = [];
  Object.keys(map[synthesizer]).forEach((technology: string) => {
    const codes = Object.keys(map[synthesizer][technology]).concat(Object.keys(MODIFICATIONS));
    if (codes.some((s) => s == sequence.slice(0, s.length)))
      technologies.push(technology);
  });
  return technologies;
}

export function convertSequence(sequence: string, output: {
  indexOfFirstNotValidChar: number, synthesizer: string[] | null, technology: string[] | null}) {
  if (output.indexOfFirstNotValidChar != -1) {
    return {
      // type: '',
      indexOfFirstNotValidChar: JSON.stringify(output),
      Error: undefinedInputSequence,
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.RAW_NUCLEOTIDES)) {//&& output.technology!.includes(TECHNOLOGIES.DNA)) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES, // + ' ' + TECHNOLOGIES.DNA,
      Nucleotides: sequence,
      BioSpring: asoGapmersNucleotidesToBioSpring(sequence),
      GCRS: asoGapmersNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.BIOSPRING) && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersBioSpringToNucleotides(sequence),
      BioSpring: sequence,
      GCRS: asoGapmersBioSpringToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.GCRS) && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersGcrsToNucleotides(sequence),
      BioSpring: asoGapmersGcrsToBioSpring(sequence),
      Mermade12: gcrsToMermade12(sequence),
      GCRS: sequence,
      LCMS: gcrsToLcms(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.RAW_NUCLEOTIDES) && output.technology!.includes(TECHNOLOGIES.RNA)) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' + TECHNOLOGIES.RNA,
      Nucleotides: sequence,
      BioSpring: siRnaNucleotideToBioSpringSenseStrand(sequence),
      Axolabs: siRnaNucleotideToAxolabsSenseStrand(sequence),
      GCRS: siRnaNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.BIOSPRING) && output.technology!.includes(TECHNOLOGIES.SI_RNA)) {
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaBioSpringToNucleotides(sequence),
      BioSpring: sequence,
      Axolabs: siRnaBioSpringToAxolabs(sequence),
      GCRS: siRnaBioSpringToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS)) {
    return {
      type: SYNTHESIZERS.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaAxolabsToNucleotides(sequence),
      BioSpring: siRnaAxolabsToBioSpring(sequence),
      Axolabs: sequence,
      GCRS: siRnaAxolabsToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.GCRS) && output.technology!.includes(TECHNOLOGIES.SI_RNA)) {
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaGcrsToNucleotides(sequence),
      BioSpring: siRnaGcrsToBioSpring(sequence),
      Axolabs: siRnaGcrsToAxolabs(sequence),
      MM12: gcrsToMermade12(sequence),
      GCRS: sequence,
      LCMS: gcrsToLcms(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.GCRS)) {
    return {
      type: SYNTHESIZERS.GCRS,
      Nucleotides: gcrsToNucleotides(sequence),
      GCRS: sequence,
      Mermade12: gcrsToMermade12(sequence),
      LCMS: gcrsToLcms(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12)) {
    return {
      type: SYNTHESIZERS.MERMADE_12,
      Nucleotides: noTranslationTableAvailable,
      GCRS: noTranslationTableAvailable,
      Mermade12: sequence,
    };
  }
  return {
    type: undefinedInputSequence,
    Nucleotides: undefinedInputSequence,
  };
}
