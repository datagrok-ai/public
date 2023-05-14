import {SYNTHESIZERS as FORMAT, TECHNOLOGIES} from '../const';
import {FormatConverter} from './format-converter';

const NO_TRANSLATION_MSG = 'No translation table available';
export const UNDEFINED_SEQ_MSG = 'Type of input sequence is undefined';

export function convertSequence(sequence: string, indexOfFirstInvalidChar: number, format: FORMAT | null): {[key: string]: string} {
  if (indexOfFirstInvalidChar !== -1) {
    return {
      indexOfFirstInvalidChar: indexOfFirstInvalidChar.toString(),
      Error: UNDEFINED_SEQ_MSG,
    };
  }
  if (format === FORMAT.NUCLEOTIDES ) {
    const converter = new FormatConverter(sequence, FORMAT.NUCLEOTIDES);
    return {
      type: FORMAT.NUCLEOTIDES,
      Nucleotides: sequence,
      BioSpring: converter.convert(FORMAT.BIOSPRING),
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (format === FORMAT.BIOSPRING) {
    const converter = new FormatConverter(sequence, FORMAT.BIOSPRING);
    return {
      type: FORMAT.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: converter.convert(FORMAT.NUCLEOTIDES),
      BioSpring: sequence,
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (format === FORMAT.GCRS) {
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
  if (format === FORMAT.AXOLABS) {
    const converter = new FormatConverter(sequence, FORMAT.AXOLABS);
    return {
      type: FORMAT.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: converter.convert(FORMAT.NUCLEOTIDES),
      BioSpring: converter.convert(FORMAT.BIOSPRING),
      Axolabs: sequence,
      GCRS: converter.convert(FORMAT.GCRS),
    };
  }
  if (format === FORMAT.MERMADE_12) {
    return {
      type: FORMAT.MERMADE_12,
      Nucleotides: NO_TRANSLATION_MSG,
      GCRS: NO_TRANSLATION_MSG,
      Mermade12: sequence,
    };
  }
  return {
    type: UNDEFINED_SEQ_MSG,
    Nucleotides: UNDEFINED_SEQ_MSG,
  };
}
