import {map} from "./map";

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (let technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes;
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string): string[] {
  let synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    const codes = getAllCodesOfSynthesizer(synthesizer);
    if (codes.some((s) => s == sequence.slice(0, s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

export function isValid(sequence: string): { indexOfFirstNotValidCharacter: number, expectedSynthesizer: string | null } {
  let possibleSynthesizers = getListOfPossibleSynthesizersByFirstMatchedCode(sequence);
  if (possibleSynthesizers.length == 0)
    return { indexOfFirstNotValidCharacter: 0, expectedSynthesizer: null };

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

  return {
    indexOfFirstNotValidCharacter: indexOfFirstNotValidCharacter,
    expectedSynthesizer: possibleSynthesizers[outputIndices.indexOf(indexOfExpectedSythesizer)]
  };
}
