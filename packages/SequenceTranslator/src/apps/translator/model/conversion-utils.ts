import {DEFAULT_FORMATS, NUCLEOTIDES} from '../../common/model/const';
import {NUCLEOTIDES_FORMAT} from '../view/const';
import {UNKNOWN_SYMBOL} from './const';
import {FormatConverter} from './format-converter';
import {CODES_TO_HELM_DICT} from '../../common/data-loader/json-loader';
import {MonomerLibWrapper} from '../../common/monomer-lib/lib-wrapper';

export function getTranslatedSequences(sequence: string, indexOfFirstInvalidChar: number, sourceFormat: string): {[key: string]: string} {
  const supportedFormats = Object.keys(CODES_TO_HELM_DICT).concat([DEFAULT_FORMATS.HELM]) as string[];

  if (!sequence || (indexOfFirstInvalidChar !== -1 && sourceFormat !== DEFAULT_FORMATS.HELM))
    return {};

  if (!supportedFormats.includes(sourceFormat))
    throw new Error(`${sourceFormat} format is not supported by SequenceTranslator`);

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
    }).filter(([_, translation]) => translation)
  );
  const helm = (sourceFormat === DEFAULT_FORMATS.HELM) ? sequence : result[DEFAULT_FORMATS.HELM];
  const nucleotides = getNucleotidesSequence(helm, MonomerLibWrapper.getInstance());
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

// todo: remove after refactoring as a workaround
export function convert(sequence: string, sourceFormat: string, targetFormat: string): string | null {
  const converter = new FormatConverter(sequence, sourceFormat);
  if (targetFormat === NUCLEOTIDES_FORMAT) {
    const helm = converter.convertTo(DEFAULT_FORMATS.HELM);
    const nucleotides = getNucleotidesSequence(helm, MonomerLibWrapper.getInstance());
    return nucleotides;
  }

  return converter.convertTo(targetFormat);
}

export function getSupportedTargetFormats(): string[] {
  const supportedTargetFormats = Object.keys(CODES_TO_HELM_DICT).concat([DEFAULT_FORMATS.HELM, NUCLEOTIDES_FORMAT]).sort() as string[];
  return supportedTargetFormats;
}
