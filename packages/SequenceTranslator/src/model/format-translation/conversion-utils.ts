import {FORMAT} from '../const';
import {FormatConverter} from './format-converter';
import {codesToHelmDictionary} from '../data-loading-utils/json-loader';

export function getTranslatedSequences(sequence: string, indexOfFirstInvalidChar: number, sourceFormat: FORMAT): {[key: string]: string} {
  const supportedFormats = Object.keys(codesToHelmDictionary).concat([FORMAT.NUCLEOTIDES, FORMAT.HELM]) as FORMAT[];

  if (!sequence)
    return {};

  if (indexOfFirstInvalidChar !== -1 && sourceFormat !== FORMAT.HELM) {
    return {
      'Error: invalid sequence': '',
    };
  }
  if (!supportedFormats.includes(sourceFormat))
    throw new Error(`${sourceFormat} format is not supported by SequenceTranslator`)

  const outputFormats = supportedFormats
    .filter((el) => el != sourceFormat)
    .sort((a, b) => a.localeCompare(b));
  const converter = new FormatConverter(sequence, sourceFormat);
  const result = Object.fromEntries(
    outputFormats.map((format) => {
      let translation;
      try {
        translation = converter.convertTo(format);
      } catch {
        translation = null;
      }
      return [format, translation ? translation : 'No translation available'];
    })
  )
  console.log(`${sourceFormat}:`, result);
  return result;
}
