import {SYNTHESIZERS, TECHNOLOGIES, DELIMITER, NUCLEOTIDES} from '../const';
import {
  asoGapmersNucleotidesToBioSpring, asoGapmersNucleotidesToGcrs,
  asoGapmersBioSpringToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12, siRnaNucleotideToBioSpringSenseStrand,
  siRnaNucleotideToAxolabsSenseStrand, siRnaNucleotidesToGcrs, siRnaBioSpringToNucleotides,
  siRnaBioSpringToAxolabs, siRnaBioSpringToGcrs, siRnaAxolabsToNucleotides,
  siRnaAxolabsToBioSpring, siRnaAxolabsToGcrs, siRnaGcrsToNucleotides,
  siRnaGcrsToBioSpring, siRnaGcrsToAxolabs, gcrsToNucleotides, gcrsToLcms
} from '../../hardcode-to-be-eliminated/converters';

import {FormatDetector} from '../parsing-validation-utils/format-detector';

const noTranslationTableAvailable = 'No translation table available';
export const undefinedInputSequence = 'Type of input sequence is undefined';

export function isValidSequence(sequence: string, format: string | null): {
  indexOfFirstInvalidChar: number,
  synthesizer: string[] | null,
} {
  const formatDetector = new FormatDetector(sequence);
  const synthesizer = format ? format : formatDetector.getFormat();
  if (!synthesizer)
    return {indexOfFirstInvalidChar: 0, synthesizer: null};
  const indexOfFirstInvalidChar = formatDetector.getInvalidCodeIndex();
  return {
    indexOfFirstInvalidChar: indexOfFirstInvalidChar,
    synthesizer: [synthesizer],
  };
}

export function convertSequence(sequence: string, output: {
  indexOfFirstInvalidChar: number, synthesizer: string[] | null
}) {
  if (output.indexOfFirstInvalidChar !== -1) {
    return {
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
