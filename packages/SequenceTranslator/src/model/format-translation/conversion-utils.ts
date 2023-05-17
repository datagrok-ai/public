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
      BioSpring: converter.convertTo(FORMAT.BIOSPRING),
      GCRS: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (format === FORMAT.BIOSPRING) {
    const converter = new FormatConverter(sequence, FORMAT.BIOSPRING);
    return {
      type: FORMAT.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: converter.convertTo(FORMAT.NUCLEOTIDES),
      BioSpring: sequence,
      GCRS: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (format === FORMAT.GCRS) {
    const converter = new FormatConverter(sequence, FORMAT.GCRS);
    return {
      type: FORMAT.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      Nucleotides: converter.convertTo(FORMAT.NUCLEOTIDES),
      BioSpring: converter.convertTo(FORMAT.BIOSPRING),
      Axolabs: converter.convertTo(FORMAT.AXOLABS),
      Mermade12: converter.convertTo(FORMAT.MERMADE_12),
      GCRS: sequence,
      LCMS: converter.convertTo(FORMAT.LCMS),
    };
  }
  if (format === FORMAT.AXOLABS) {
    const converter = new FormatConverter(sequence, FORMAT.AXOLABS);
    return {
      type: FORMAT.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      Nucleotides: converter.convertTo(FORMAT.NUCLEOTIDES),
      BioSpring: converter.convertTo(FORMAT.BIOSPRING),
      Axolabs: sequence,
      GCRS: converter.convertTo(FORMAT.GCRS),
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
