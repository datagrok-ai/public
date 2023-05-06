import {SYNTHESIZERS as FORMAT, TECHNOLOGIES} from '../const';
import {FormatConverter} from './format-converter';

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
  if (output.synthesizer!.includes(FORMAT.NUCLEOTIDES) /*&& output.technology!.includes(TECHNOLOGIES.DNA)*/) {
    const converter = new FormatConverter(sequence, FORMAT.NUCLEOTIDES);
    return {
      type: FORMAT.NUCLEOTIDES, // + ' ' + TECHNOLOGIES.DNA,
      Nucleotides: sequence,
      BioSpring: converter.convert(FORMAT.BIOSPRING),
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (output.synthesizer!.includes(FORMAT.BIOSPRING)) {
    const converter = new FormatConverter(sequence, FORMAT.BIOSPRING);
    // && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    return {
      type: FORMAT.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: converter.convert(FORMAT.NUCLEOTIDES),
      BioSpring: sequence,
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (output.synthesizer!.includes(FORMAT.GCRS)) { // && output.technology!.includes(TECHNOLOGIES.ASO_GAPMERS)) {
    const converter = new FormatConverter(sequence, FORMAT.GCRS);
    return {
      type: FORMAT.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: converter.convert(FORMAT.NUCLEOTIDES),
      BioSpring: converter.convert(FORMAT.BIOSPRING),
      Axolabs: converter.convert(FORMAT.AXOLABS),
      Mermade12: converter.convert(FORMAT.MERMADE_12),
      GCRS: sequence,
      LCMS: converter.convert(FORMAT.LCMS),
    };
  }
  if (output.synthesizer!.includes(FORMAT.AXOLABS)) {
    const converter = new FormatConverter(sequence, FORMAT.AXOLABS);
    return {
      type: FORMAT.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: converter.convert(FORMAT.NUCLEOTIDES),
      BioSpring: converter.convert(FORMAT.BIOSPRING),
      Axolabs: sequence,
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (output.synthesizer!.includes(FORMAT.MERMADE_12)) {
    return {
      type: FORMAT.MERMADE_12,
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
