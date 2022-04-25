import {map, SYNTHESIZERS, TECHNOLOGIES, MODIFICATIONS} from './map';
import {asoGapmersNucleotidesToBioSpring, asoGapmersNucleotidesToGcrs,
  asoGapmersBioSpringToNucleotides, asoGapmersBioSpringToGcrs, asoGapmersGcrsToNucleotides,
  asoGapmersGcrsToBioSpring, gcrsToMermade12, siRnaNucleotideToBioSpringSenseStrand,
  siRnaNucleotideToAxolabsSenseStrand, siRnaNucleotidesToGcrs, siRnaBioSpringToNucleotides,
  siRnaBioSpringToAxolabs, siRnaBioSpringToGcrs, siRnaAxolabsToNucleotides,
  siRnaAxolabsToBioSpring, siRnaAxolabsToGcrs, siRnaGcrsToNucleotides,
  siRnaGcrsToBioSpring, siRnaGcrsToAxolabs, gcrsToNucleotides} from './converters';

const noTranslationTableAvailable = 'No translation table available';
export const undefinedInputSequence = 'Type of input sequence is undefined';

export function isValidSequence(sequence: string): {
  indexOfFirstNotValidChar: number,
  synthesizer: string | null,
  technology: string | null
} {
  const possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence);
  if (possibleSynthesizers.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: null, technology: null};

  let outputIndices = Array(possibleSynthesizers.length).fill(0);

  const firstUniqueCharacters = ['r', 'd'];
  const nucleotides = ['A', 'U', 'T', 'C', 'G'];

  possibleSynthesizers.forEach((synthesizer, synthesizerIndex) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    while (outputIndices[synthesizerIndex] < sequence.length) {
      const matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[synthesizerIndex], outputIndices[synthesizerIndex] + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndices[synthesizerIndex] > 1 &&
        nucleotides.includes(sequence[outputIndices[synthesizerIndex]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[synthesizerIndex] - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[synthesizerIndex] + 1]) &&
        nucleotides.includes(sequence[outputIndices[synthesizerIndex]])
      ) {
        outputIndices[synthesizerIndex]++;
        break;
      }

      outputIndices[synthesizerIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedSythesizer = Math.max(...outputIndices);
  const indexOfFirstNotValidChar = (indexOfExpectedSythesizer == sequence.length) ? -1 : indexOfExpectedSythesizer;
  const synthesizer = possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)];
  if (indexOfFirstNotValidChar != -1) {
    return {
      indexOfFirstNotValidChar: indexOfFirstNotValidChar,
      synthesizer: synthesizer,
      technology: null,
    };
  }

  const possibleTechnologies = getListOfPossibleTechnologiesByFirstMatchedCode(sequence, synthesizer);
  if (possibleTechnologies.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: null, technology: null};

  outputIndices = Array(possibleTechnologies.length).fill(0);

  possibleTechnologies.forEach((technology: string, technologyIndex: number) => {
    const codes = Object.keys(map[synthesizer][technology]);
    while (outputIndices[technologyIndex] < sequence.length) {
      const matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[technologyIndex], outputIndices[technologyIndex] + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndices[technologyIndex] > 1 &&
        nucleotides.includes(sequence[outputIndices[technologyIndex]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[technologyIndex] - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[technologyIndex] + 1]) &&
        nucleotides.includes(sequence[outputIndices[technologyIndex]])
      ) {
        outputIndices[technologyIndex]++;
        break;
      }

      outputIndices[technologyIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedTechnology = Math.max(...outputIndices);
  const expectedTechnology = possibleTechnologies[outputIndices.indexOf(indexOfExpectedTechnology)];

  return {
    indexOfFirstNotValidChar: indexOfFirstNotValidChar,
    synthesizer: synthesizer,
    technology: expectedTechnology,
  };
}

function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (const technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes.concat(Object.keys(MODIFICATIONS));
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
  indexOfFirstNotValidChar: number, synthesizer: string | null, technology: string | null}) {
  if (output.indexOfFirstNotValidChar != -1) {
    return {
      // type: '',
      indexOfFirstNotValidChar: JSON.stringify(output),
      Error: undefinedInputSequence,
    };
  }
  if (output.synthesizer == SYNTHESIZERS.RAW_NUCLEOTIDES && output.technology == TECHNOLOGIES.DNA) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' + TECHNOLOGIES.DNA,
      Nucleotides: sequence,
      BioSpring: asoGapmersNucleotidesToBioSpring(sequence),
      GCRS: asoGapmersNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.BIOSPRING && output.technology == TECHNOLOGIES.ASO_GAPMERS) {
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersBioSpringToNucleotides(sequence),
      BioSpring: sequence,
      GCRS: asoGapmersBioSpringToGcrs(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.GCRS && output.technology == TECHNOLOGIES.ASO_GAPMERS) {
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersGcrsToNucleotides(sequence),
      BioSpring: asoGapmersGcrsToBioSpring(sequence),
      Mermade12: gcrsToMermade12(sequence),
      GCRS: sequence,
    };
  }
  if (output.synthesizer == SYNTHESIZERS.RAW_NUCLEOTIDES && output.technology == TECHNOLOGIES.RNA) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' + TECHNOLOGIES.RNA,
      Nucleotides: sequence,
      BioSpring: siRnaNucleotideToBioSpringSenseStrand(sequence),
      Axolabs: siRnaNucleotideToAxolabsSenseStrand(sequence),
      GCRS: siRnaNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.BIOSPRING && output.technology == TECHNOLOGIES.SI_RNA) {
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaBioSpringToNucleotides(sequence),
      BioSpring: sequence,
      Axolabs: siRnaBioSpringToAxolabs(sequence),
      GCRS: siRnaBioSpringToGcrs(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.AXOLABS && output.technology == TECHNOLOGIES.SI_RNA) {
    return {
      type: SYNTHESIZERS.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaAxolabsToNucleotides(sequence),
      BioSpring: siRnaAxolabsToBioSpring(sequence),
      Axolabs: sequence,
      GCRS: siRnaAxolabsToGcrs(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.GCRS && output.technology == TECHNOLOGIES.SI_RNA) {
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: siRnaGcrsToNucleotides(sequence),
      BioSpring: siRnaGcrsToBioSpring(sequence),
      Axolabs: siRnaGcrsToAxolabs(sequence),
      MM12: gcrsToMermade12(sequence),
      GCRS: sequence,
    };
  }
  if (output.synthesizer == SYNTHESIZERS.GCRS) {
    return {
      type: SYNTHESIZERS.GCRS,
      Nucleotides: gcrsToNucleotides(sequence),
      GCRS: sequence,
      Mermade12: gcrsToMermade12(sequence),
    };
  }
  if (output.synthesizer == SYNTHESIZERS.MERMADE_12) {
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
