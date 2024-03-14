import {NUCLEOTIDES} from '../const';
import {MonomerLibWrapper} from '../../monomer-lib/lib-wrapper';
import {sortByReverseLength} from '../helpers';
import {DEFAULT_FORMATS} from '../const';

export class SequenceValidator {
  constructor(private sequence: string) {
    this.libWrapper = MonomerLibWrapper.getInstance();
  };
  private libWrapper: MonomerLibWrapper;

  getInvalidCodeIndex(format: string): number {
    if (format === DEFAULT_FORMATS.HELM)
      return this.sequence.length;
    const firstUniqueCharacters = ['r', 'd']; // what for?
    const codes = sortByReverseLength(
      this.libWrapper.getCodesByFormat(format)
    );
    let indexOfFirstInvalidChar = 0;
    while (indexOfFirstInvalidChar < this.sequence.length) {
      const matchedCode = codes.find((code) => {
        const subSequence = this.sequence.substring(indexOfFirstInvalidChar, indexOfFirstInvalidChar + code.length);
        return code === subSequence;
      });

      if (!matchedCode) break;

      // todo: refactor the vague condition
      if ( // for mistake pattern 'rAA'
        indexOfFirstInvalidChar > 1 &&
        NUCLEOTIDES.includes(this.sequence[indexOfFirstInvalidChar]) &&
        firstUniqueCharacters.includes(this.sequence[indexOfFirstInvalidChar - 2])
      ) break;

      if ( // for mistake pattern 'ArA'
        firstUniqueCharacters.includes(this.sequence[indexOfFirstInvalidChar + 1]) &&
        NUCLEOTIDES.includes(this.sequence[indexOfFirstInvalidChar])
      ) {
        indexOfFirstInvalidChar++;
        break;
      }
      indexOfFirstInvalidChar += matchedCode.length;
    }
    if (indexOfFirstInvalidChar === this.sequence.length)
      indexOfFirstInvalidChar = -1;
    return indexOfFirstInvalidChar;
  }

  isValidSequence(format: string): boolean {
    return this.getInvalidCodeIndex(format) === -1;
  }
}
