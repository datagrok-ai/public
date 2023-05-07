import {FormatDetector} from '../parsing-validation-utils/format-detector';

export function isValidSequence(sequence: string, format: string | null): {
  indexOfFirstInvalidChar: number,
  synthesizer: string[] | null,
} {
  const formatDetector = new FormatDetector(sequence);
  // we get format if not specified
  const synthesizer = format ? format : formatDetector.getFormat();
  // return invalid status if format detection failed
  if (!synthesizer)
    return {indexOfFirstInvalidChar: 0, synthesizer: null};
  // return (possibly) detected format AND get the invalid code idx
  const indexOfFirstInvalidChar = formatDetector.getInvalidCodeIndex();
  return {
    indexOfFirstInvalidChar: indexOfFirstInvalidChar,
    synthesizer: [synthesizer],
  };
}
