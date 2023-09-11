import {DEFAULT_FORMATS, NUCLEOTIDES} from '../const';
import {UNKNOWN_SYMBOL} from './const';
import {FormatConverter} from './format-converter';
import {codesToHelmDictionary} from '../data-loading-utils/json-loader';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';

export function getTranslatedSequences(sequence: string, indexOfFirstInvalidChar: number, sourceFormat: string): {[key: string]: string} {
  const supportedFormats = Object.keys(codesToHelmDictionary).concat([DEFAULT_FORMATS.HELM]) as string[];

  if (!sequence || (indexOfFirstInvalidChar !== -1 && sourceFormat !== DEFAULT_FORMATS.HELM))
    return {};

  if (!supportedFormats.includes(sourceFormat))
    throw new Error(`${sourceFormat} format is not supported by SequenceTranslator`)

  const outputFormats = supportedFormats.filter((el) => el != sourceFormat)
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
      return [format, translation];
    }).filter(([format, translation]) => translation)
  )
  const nucleotides = getNucleotidesSequence(result[DEFAULT_FORMATS.HELM], MonomerLibWrapper.getInstance());
  if (nucleotides)
    result['Nucleotides'] = nucleotides;
  return result;
}

export function getNucleotidesSequence(helmString: string, monomerLib: MonomerLibWrapper): string | null {
  const re = new RegExp('\\([^()]*\\)', 'g');
  const branches = helmString.match(re);
  if (!branches)
    return null;
  const nucleotides = branches!.map((branch) => {
    const stripped = branch.replace(/[\[\]()]/g, '');
    if (NUCLEOTIDES.includes(stripped))
      return stripped;
    return monomerLib.getNaturalAnalogBySymbol(stripped);
  }).map((el) => el ? el : UNKNOWN_SYMBOL).join('');
  return nucleotides;
}
