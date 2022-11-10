import {map, SYNTHESIZERS, NUCLEOTIDES, FIRST_UNIQUE_CHARACTERS} from './constants';
import {sortByStringLengthInDescOrder, getAllCodesOfSynthesizer} from './helpers';

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string, additionalCodes: string[]): string[] {
  const synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    const matched = additionalCodes.find((s) => s == sequence.slice(0, s.length));
    const cutoffIndex = (matched == null) ? 0 : matched.length;
    const codes = getAllCodesOfSynthesizer(synthesizer);
    if (codes.some((s) => s == sequence.slice(cutoffIndex, cutoffIndex + s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

function possibleTechnologiesByFirstMatchedCode(sequence: string, synthesizer: string, additCodes: string[]): string[] {
  const technologies: string[] = [];
  Object.keys(map[synthesizer]).forEach((technology: string) => {
    const matched = additCodes.find((s) => s == sequence.slice(0, s.length));
    const cutoffIndex = (matched == null) ? 0 : matched.length;
    const codes = Object.keys(map[synthesizer][technology]);
    if (codes.some((s) => s == sequence.slice(cutoffIndex, cutoffIndex + s.length)))
      technologies.push(technology);
  });
  return technologies;
}

export function sequenceValidation(sequence: string, additionalCodes: string[]): {
  indexOfFirstNotValidChar: number,
  synthesizer: string[] | null,
  technology: string[] | null
} {
  const sortedAdditionalCodes = sortByStringLengthInDescOrder(additionalCodes);

  const possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence, sortedAdditionalCodes);
  if (possibleSynthesizers.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: null, technology: null};

  let outputIndices = Array(possibleSynthesizers.length).fill(0);

  possibleSynthesizers.forEach((synthesizer, synthesizerIndex) => {
    const codes = sortByStringLengthInDescOrder(getAllCodesOfSynthesizer(synthesizer)
      .concat(sortedAdditionalCodes));
    while (outputIndices[synthesizerIndex] < sequence.length) {
      const matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[synthesizerIndex], outputIndices[synthesizerIndex] + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndices[synthesizerIndex] > 1 &&
        NUCLEOTIDES.includes(sequence[outputIndices[synthesizerIndex]]) &&
        FIRST_UNIQUE_CHARACTERS.includes(sequence[outputIndices[synthesizerIndex] - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        FIRST_UNIQUE_CHARACTERS.includes(sequence[outputIndices[synthesizerIndex] + 1]) &&
        NUCLEOTIDES.includes(sequence[outputIndices[synthesizerIndex]])
      ) {
        outputIndices[synthesizerIndex]++;
        break;
      }

      outputIndices[synthesizerIndex] += matchedCode.length;
    }
  });

  const indexOfExpectedSythesizer = Math.max(...outputIndices); //Math.max.apply(Math, outputIndices);
  const indexOfFirstNotValidChar = (indexOfExpectedSythesizer == sequence.length) ? -1 : indexOfExpectedSythesizer;
  const synthesizer = possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)];
  if (indexOfFirstNotValidChar != -1) {
    return {
      indexOfFirstNotValidChar: indexOfFirstNotValidChar,
      synthesizer: [synthesizer],
      technology: null,
    };
  }

  const possibleTechnologies = possibleTechnologiesByFirstMatchedCode(
    sequence, synthesizer, sortedAdditionalCodes);
  if (possibleTechnologies.length == 0)
    return {indexOfFirstNotValidChar: 0, synthesizer: null, technology: null};

  outputIndices = Array(possibleTechnologies.length).fill(0);

  possibleTechnologies.forEach((technology: string, technologyIndex: number) => {
    const codes = sortByStringLengthInDescOrder(
      Object.keys(map[synthesizer][technology]).concat(sortedAdditionalCodes));
    while (outputIndices[technologyIndex] < sequence.length) {
      const matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[technologyIndex], outputIndices[technologyIndex] + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
        outputIndices[technologyIndex] > 1 &&
        NUCLEOTIDES.includes(sequence[outputIndices[technologyIndex]]) &&
        FIRST_UNIQUE_CHARACTERS.includes(sequence[outputIndices[technologyIndex] - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        FIRST_UNIQUE_CHARACTERS.includes(sequence[outputIndices[technologyIndex] + 1]) &&
        NUCLEOTIDES.includes(sequence[outputIndices[technologyIndex]])
      ) {
        outputIndices[technologyIndex]++;
        break;
      }

      outputIndices[technologyIndex] += matchedCode.length;
    }
  });

  const indexOfTechnology = Math.max(...outputIndices); //Math.max.apply(Math, outputIndices);
  const technology = possibleTechnologies[outputIndices.indexOf(indexOfTechnology)];

  return {
    indexOfFirstNotValidChar: indexOfFirstNotValidChar,
    synthesizer: [synthesizer],
    technology: [technology],
  };
}

export function validate(sequence: string, additionalCodes: string[]): number {
  const codes = sortByStringLengthInDescOrder(getAllCodesOfSynthesizer(SYNTHESIZERS.GCRS).concat(additionalCodes));
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((c) => c == sequence.slice(i, i + c.length));

    if (matchedCode == null)
      break;

    if (i > 1 && NUCLEOTIDES.includes(sequence[i]) && FIRST_UNIQUE_CHARACTERS.includes(sequence[i - 2]))
      break;

    // if (firstUniqueCharacters.includes(sequence[i + 1]) && nucleotides.includes(sequence[i])) {
    //   i++;
    //   break;
    // }

    i += matchedCode.length;
  }
  return (i == sequence.length) ? -1 : i;
}

export function isValidSequence(indicesOfFirstNotValidCharacter: number): boolean {
  return indicesOfFirstNotValidCharacter == -1;
}
