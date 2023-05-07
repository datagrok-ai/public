/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {sortByReverseLength} from '../helpers';
import {SYNTHESIZERS, NUCLEOTIDES} from '../const';
import {MonomerLibWrapper} from '../monomer-lib-utils/lib-wrapper';

export class FormatDetector {
  constructor (private sequence: string) {
    this.libWrapper = MonomerLibWrapper.getInstance();
  };

  private libWrapper: MonomerLibWrapper;

  getInvalidCodeIndex(format?: string): number {
    const encodingFormat = (format === undefined) ? this.getFormat() : format;
    if (!encodingFormat)
      return 0;
    return this.computeInvalidCharIdx(encodingFormat);
  }

  getFormat(): string | null {
    const possibleFormats = this.getListOfPossibleSynthesizersByFirstMatchedCode();
    if (possibleFormats.length === 0)
      return null;

    const outputIndices = Array(possibleFormats.length).fill(0);
    for (let i = 0; i < possibleFormats.length; ++i) {
      const format = possibleFormats[i];
      outputIndices[i] = this.computeInvalidCharIdx(format);
    }
    const formatIdx = (outputIndices.some((idx) => idx === -1)) ? -1 : Math.max(...outputIndices);
    return possibleFormats[outputIndices.indexOf(formatIdx)];
  }

  isValidSequence(format?: string): boolean {
    return this.getInvalidCodeIndex(format) === -1;
  }

  private computeInvalidCharIdx(format: string): number {
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
      indexOfFirstInvalidChar = -1
    return indexOfFirstInvalidChar;
  }

  private getListOfPossibleSynthesizersByFirstMatchedCode(): string[] {
    const sequence = this.sequence;
    let synthesizers: string[] = [];
    for (const synthesizer of Object.values(SYNTHESIZERS)) {
      let codes = sortByReverseLength(this.libWrapper.getCodesByFormat(synthesizer));
      let start = 0;
      for (let i = 0; i < sequence.length; i++) {
        if (sequence[i] === ')' && i !== sequence.length - 1) {
          start = i + 1;
          break;
        }
      }
      if (codes.some((s: string) => s === sequence.slice(start, start + s.length)))
        synthesizers.push(synthesizer);
    }
    return synthesizers;
  }
}
