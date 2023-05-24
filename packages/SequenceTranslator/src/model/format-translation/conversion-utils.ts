import {SYNTHESIZERS as FORMAT, TECHNOLOGIES} from '../const';
import {FormatConverter} from './format-converter';

const NO_TRANSLATION_MSG = 'No translation table available';
export const UNDEFINED_SEQ_MSG = 'Type of input sequence is undefined';

export function convertSequence(sequence: string, indexOfFirstInvalidChar: number, sourceFormat: FORMAT | null): {[key: string]: string} {
  if (indexOfFirstInvalidChar !== -1 && sourceFormat !== FORMAT.HELM) {
    return {
      indexOfFirstInvalidChar: indexOfFirstInvalidChar.toString(),
      Error: UNDEFINED_SEQ_MSG,
    };
  }
  if (sourceFormat === FORMAT.NUCLEOTIDES) {
    const converter = new FormatConverter(sequence, FORMAT.NUCLEOTIDES);
    return {
      type: FORMAT.NUCLEOTIDES,
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.BIOSPRING) {
    const converter = new FormatConverter(sequence, FORMAT.BIOSPRING);
    return {
      type: FORMAT.BIOSPRING + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.GCRS) {
    const converter = new FormatConverter(sequence, FORMAT.GCRS);
    return {
      type: FORMAT.GCRS + ' ' + TECHNOLOGIES.ASO_GAPMERS,
      [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.AXOLABS]: converter.convertTo(FORMAT.AXOLABS),
      [FORMAT.MERMADE_12]: converter.convertTo(FORMAT.MERMADE_12),
      [FORMAT.LCMS]: converter.convertTo(FORMAT.LCMS),
      [FORMAT.HELM]: converter.convertTo(FORMAT.HELM)
    };
  }
  if (sourceFormat === FORMAT.AXOLABS) {
    const converter = new FormatConverter(sequence, FORMAT.AXOLABS);
    return {
      type: FORMAT.AXOLABS + ' ' + TECHNOLOGIES.SI_RNA,
      [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.MERMADE_12) {
    return {
      type: FORMAT.MERMADE_12,
      [FORMAT.NUCLEOTIDES]: NO_TRANSLATION_MSG,
      [FORMAT.GCRS]: NO_TRANSLATION_MSG,
    };
  }
  if (sourceFormat === FORMAT.HELM) {
    const gcrsSequence = (new FormatConverter(sequence, sourceFormat)).convertTo(FORMAT.GCRS);
    const converter = new FormatConverter(gcrsSequence, FORMAT.GCRS);
    return {
      type: FORMAT.HELM,
      [FORMAT.GCRS]: gcrsSequence,
      [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.AXOLABS]: converter.convertTo(FORMAT.AXOLABS),
      [FORMAT.MERMADE_12]: converter.convertTo(FORMAT.MERMADE_12),
      [FORMAT.LCMS]: converter.convertTo(FORMAT.LCMS),
    }
  }
  return {
    type: UNDEFINED_SEQ_MSG,
    [FORMAT.NUCLEOTIDES]: UNDEFINED_SEQ_MSG,
  };
}
