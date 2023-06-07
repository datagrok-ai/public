import {FORMAT} from '../const';
import {FormatConverter} from './format-converter';

const NO_TRANSLATION_MSG = 'No translation table available';
export const UNDEFINED_SEQ_MSG = 'Type of input sequence is undefined';

export function getTranslatedSequences(sequence: string, indexOfFirstInvalidChar: number, sourceFormat: FORMAT | null): {[key: string]: string} {
  if (indexOfFirstInvalidChar !== -1 && sourceFormat !== FORMAT.HELM) {
    return {
      indexOfFirstInvalidChar: indexOfFirstInvalidChar.toString(),
      Error: UNDEFINED_SEQ_MSG,
    };
  }
  if (sourceFormat === FORMAT.NUCLEOTIDES) {
    const converter = new FormatConverter(sequence, FORMAT.NUCLEOTIDES);
    return {
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.BIOSPRING) {
    const converter = new FormatConverter(sequence, FORMAT.BIOSPRING);
    return {
      // [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.HELM]: converter.convertTo(FORMAT.HELM),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.GCRS) {
    const converter = new FormatConverter(sequence, FORMAT.GCRS);
    return {
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
      [FORMAT.NUCLEOTIDES]: converter.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.HELM]: converter.convertTo(FORMAT.HELM),
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.GCRS]: converter.convertTo(FORMAT.GCRS),
    };
  }
  if (sourceFormat === FORMAT.MERMADE_12) {
    return {
      [FORMAT.NUCLEOTIDES]: NO_TRANSLATION_MSG,
      [FORMAT.GCRS]: NO_TRANSLATION_MSG,
    };
  }
  if (sourceFormat === FORMAT.HELM) {
    const converter = new FormatConverter(sequence, sourceFormat);
    const gcrsSequence = converter.convertTo(FORMAT.GCRS);
    const fromGcrs = new FormatConverter(gcrsSequence, FORMAT.GCRS);
    return {
      [FORMAT.GCRS]: gcrsSequence,
      [FORMAT.NUCLEOTIDES]: fromGcrs.convertTo(FORMAT.NUCLEOTIDES),
      [FORMAT.BIOSPRING]: converter.convertTo(FORMAT.BIOSPRING),
      [FORMAT.AXOLABS]: converter.convertTo(FORMAT.AXOLABS),
      [FORMAT.MERMADE_12]: fromGcrs.convertTo(FORMAT.MERMADE_12),
      [FORMAT.LCMS]: fromGcrs.convertTo(FORMAT.LCMS),
    }
  }
  return {
    [FORMAT.NUCLEOTIDES]: UNDEFINED_SEQ_MSG,
  };
}
