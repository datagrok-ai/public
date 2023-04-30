import {
  map,
} from '../../hardcode-to-be-eliminated/map';
import {SYNTHESIZERS, TECHNOLOGIES, DELIMITER, NUCLEOTIDES} from '../const';
import {sortByStringLengthInDescendingOrder} from '../helpers';
import {
  asoGapmersNucleotidesToBioSpring, asoGapmersNucleotidesToGcrs,
  asoGapmersBioSpringToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12, siRnaNucleotideToBioSpringSenseStrand,
  siRnaNucleotideToAxolabsSenseStrand, siRnaNucleotidesToGcrs, siRnaBioSpringToNucleotides,
  siRnaBioSpringToAxolabs, siRnaBioSpringToGcrs, siRnaAxolabsToNucleotides,
  siRnaAxolabsToBioSpring, siRnaAxolabsToGcrs, siRnaGcrsToNucleotides,
  siRnaGcrsToBioSpring, siRnaGcrsToAxolabs, gcrsToNucleotides, gcrsToLcms
} from '../../hardcode-to-be-eliminated/converters';

import {MonomerLibWrapper} from '../monomer-lib-utils/lib-wrapper';

const noTranslationTableAvailable = 'No translation table available';
export const undefinedInputSequence = 'Type of input sequence is undefined';

export function isValidSequence(sequence: string, format: string | null): {
  indexOfFirstInvalidChar: number,
  synthesizer: string[] | null,
} {
  const possibleSynthesizers = format === null ?
    getListOfPossibleSynthesizersByFirstMatchedCode(sequence) :
    [format];
  if (possibleSynthesizers.length === 0)
    return {indexOfFirstInvalidChar: 0, synthesizer: null};

  const outputIndices = Array(possibleSynthesizers.length).fill(0);

  const firstUniqueCharacters = ['r', 'd'];
  possibleSynthesizers.forEach(function(synthesizer, i) {
    const codes = sortByStringLengthInDescendingOrder(getAllCodesOfSynthesizer(synthesizer));
    while (outputIndices[i] < sequence.length) {
      const matchedCode = codes.find((c) => c === sequence.slice(outputIndices[i], outputIndices[i] + c.length));

      if (!matchedCode) break;

      if ( // for mistake pattern 'rAA'
        outputIndices[i] > 1 &&
        NUCLEOTIDES.includes(sequence[outputIndices[i]]) &&
        firstUniqueCharacters.includes(sequence[outputIndices[i] - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(sequence[outputIndices[i] + 1]) &&
        NUCLEOTIDES.includes(sequence[outputIndices[i]])
      ) {
        outputIndices[i]++;
        break;
      }

      outputIndices[i] += matchedCode.length;
    }
  });

  const outputIndex = Math.max(...outputIndices);
  const synthesizer = possibleSynthesizers[outputIndices.indexOf(outputIndex)];
  const indexOfFirstInvalidChar = (outputIndex === sequence.length) ? -1 : outputIndex;
  if (indexOfFirstInvalidChar !== -1) {
    return {
      indexOfFirstInvalidChar: indexOfFirstInvalidChar,
      synthesizer: [synthesizer],
    };
  }

  return {
    indexOfFirstInvalidChar: indexOfFirstInvalidChar,
    synthesizer: [synthesizer],
  };
}

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  console.log('got to get All codes')
  return MonomerLibWrapper.getInstance().getCodesByFromat(synthesizer);
  // let codes: string[] = [];
  // for (const technology of Object.keys(map[synthesizer]))
  //   codes = codes.concat(Object.keys(map[synthesizer][technology]));
  // return codes.concat(Object.keys(MODIFICATIONS)).concat(DELIMITER);
}

function getListOfPossibleSynthesizersByFirstMatchedCode(sequence: string): string[] {
  let synthesizers: string[] = [];
  Object.keys(map).forEach((synthesizer: string) => {
    let codes = sortByStringLengthInDescendingOrder(getAllCodesOfSynthesizer(synthesizer));
    // if (synthesizer === SYNTHESIZERS.GCRS)
    //   codes = codes.concat(gcrsCodesWithoutSmiles);
    //TODO: get first non-dropdown code when there are two modifications
    let start = 0;
    for (let i = 0; i < sequence.length; i++) {
      if (sequence[i] === ')' && i !== sequence.length - 1) {
        start = i + 1;
        break;
      }
    }
    // if (gcrsCodesWithoutSmiles.some((s: string) => s === sequence.slice(start, start + s.length)))
    //   synthesizers = [SYNTHESIZERS.GCRS];
    if (codes.some((s: string) => s === sequence.slice(start, start + s.length)))
      synthesizers.push(synthesizer);
  });
  return synthesizers;
}

export function convertSequence(sequence: string, output: {
  indexOfFirstInvalidChar: number, synthesizer: string[] | null
}) {
  if (output.indexOfFirstInvalidChar !== -1) {
    return {
      // type: '',
      indexOfFirstInvalidChar: JSON.stringify(output),
      Error: undefinedInputSequence,
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.RAW_NUCLEOTIDES) /*&& output.technology!.includes(TECHNOLOGIES.DNA)*/) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES, // + ' ' + TECHNOLOGIES.DNA,
      Nucleotides: sequence,
      BioSpring: asoGapmersNucleotidesToBioSpring(sequence),
      GCRS: asoGapmersNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.BIOSPRING)) {
    // && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    return {
      type: SYNTHESIZERS.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: asoGapmersBioSpringToNucleotides(sequence),
      BioSpring: sequence,
      GCRS: asoGapmersBioSpringToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.GCRS)) { // && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    return {
      type: SYNTHESIZERS.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: gcrsToNucleotides(sequence),
      BioSpring: siRnaGcrsToBioSpring(sequence),
      Axolabs: siRnaGcrsToAxolabs(sequence),
      Mermade12: gcrsToMermade12(sequence),
      GCRS: sequence,
      LCMS: gcrsToLcms(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.RAW_NUCLEOTIDES)) {
    // && output.technology!.includes(TECHNOLOGIES.RNA)) {
    return {
      type: SYNTHESIZERS.RAW_NUCLEOTIDES + ' ' + TECHNOLOGIES.RNA,
      Nucleotides: sequence,
      BioSpring: siRnaNucleotideToBioSpringSenseStrand(sequence),
      Axolabs: siRnaNucleotideToAxolabsSenseStrand(sequence),
      GCRS: siRnaNucleotidesToGcrs(sequence),
    };
  }
  if (output.synthesizer!.includes(SYNTHESIZERS.BIOSPRING)) { // && output.technology!.includes(TECHNOLOGIES.SI_RNA)) {
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
  if (output.synthesizer!.includes(SYNTHESIZERS.GCRS)) { // && output.technology!.includes(TECHNOLOGIES.SI_RNA)) {
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
