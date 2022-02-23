import {map} from "./map";
import { sortByStringLengthInDescendingOrder } from "./package";

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (let technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes;
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string, overhangCodes: string[]): string[] {
  let synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    let matched = overhangCodes.find((s) => s == sequence.slice(0, s.length));
    let cutoffIndex = (matched == null) ? 0 : matched.length;
    const codes = getAllCodesOfSynthesizer(synthesizer);
    if (codes.some((s) => s == sequence.slice(cutoffIndex, cutoffIndex + s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

function getListOfPossibleTechnologiesByFirstMatchedCode(sequence: string, synthesizer: string, overhangCodes: string[]): string[] {
  const technologies: string[] = [];
  Object.keys(map[synthesizer]).forEach((technology: string) => {
    let matched = overhangCodes.find((s) => s == sequence.slice(0, s.length));
    let cutoffIndex = (matched == null) ? 0 : matched.length;
    const codes = Object.keys(map[synthesizer][technology]);  // .concat(Object.keys(MODIFICATIONS));
    if (codes.some((s) => s == sequence.slice(cutoffIndex, cutoffIndex + s.length)))
      technologies.push(technology);
  });
  return technologies;
}

export function isValidSequence(sequence: string, overhangCodes: string[]): {
  indexOfFirstNotValidCharacter: number,
  expectedSynthesizer: string | null,
  expectedTechnology: string | null
} {
  const sortedOverhangCodes = sortByStringLengthInDescendingOrder(overhangCodes);

  const possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence, sortedOverhangCodes);
  if (possibleSynthesizers.length == 0)
    return {indexOfFirstNotValidCharacter: 0, expectedSynthesizer: null, expectedTechnology: null};

  let outputIndices = Array(possibleSynthesizers.length).fill(0);

  const firstUniqueCharacters = ['r', 'd'];
  const nucleotides = ['A', 'U', 'T', 'C', 'G'];

  possibleSynthesizers.forEach((synthesizer, synthesizerIndex) => {
    const codes = sortByStringLengthInDescendingOrder(getAllCodesOfSynthesizer(synthesizer).concat(sortedOverhangCodes));//.concat(sorte);
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

  const indexOfExpectedSythesizer = Math.max.apply(Math, outputIndices);
  const indexOfFirstNotValidCharacter = (indexOfExpectedSythesizer == sequence.length) ? -1 : indexOfExpectedSythesizer;
  const expectedSynthesizer = possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)];
  if (indexOfFirstNotValidCharacter != -1)
    return {
      indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
      expectedSynthesizer: expectedSynthesizer,
      expectedTechnology: null
    };

  let possibleTechnologies = getListOfPossibleTechnologiesByFirstMatchedCode(sequence, expectedSynthesizer, sortedOverhangCodes);
  if (possibleTechnologies.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedSynthesizer: null, expectedTechnology: null };

  outputIndices = Array(possibleTechnologies.length).fill(0);

  possibleTechnologies.forEach((technology: string, technologyIndex: number) => {
    let codes = sortByStringLengthInDescendingOrder(Object.keys(map[expectedSynthesizer][technology]).concat(sortedOverhangCodes));
    while (outputIndices[technologyIndex] < sequence.length) {

      let matchedCode = codes
        .find((c) => c == sequence.slice(outputIndices[technologyIndex], outputIndices[technologyIndex] + c.length));

      if (matchedCode == null)
        break;

      if ( // for mistake pattern 'rAA'
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
    expectedSynthesizer: expectedSynthesizer,
    expectedTechnology: expectedTechnology
  };
}
