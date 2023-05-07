import {FormatDetector} from '../parsing-validation-utils/format-detector';

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
