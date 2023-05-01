/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {sortByReverseLength} from '../helpers';
import {SYNTHESIZERS, NUCLEOTIDES} from '../const';
import {MonomerLibWrapper} from '../monomer-lib-utils/lib-wrapper';

const libWrapper = MonomerLibWrapper.getInstance();

export class FormatDetector {
  constructor (private sequence: string) { };

  private invalidCodeIdx?: number;

  getInvalidCodeIndex(format?: string) {
    if (this.invalidCodeIdx)
      return this.invalidCodeIdx!;
    return this.computeInvalidCodeIndex(format);
  }

  private computeInvalidCodeIndex(format?: string): number {
    const possibleSynthesizers = (format === undefined) ?
      this.getListOfPossibleSynthesizersByFirstMatchedCode() :
      [format];
    if (possibleSynthesizers.length === 0)
      return 0;
    const result = this.compute(possibleSynthesizers);
    return result.indexOfFirstInvalidChar;
  }

  getFormat(): string | null {
    const possibleSynthesizers = this.getListOfPossibleSynthesizersByFirstMatchedCode();
    if (possibleSynthesizers.length === 0)
      return null;
    const result = this.compute(possibleSynthesizers);
    return result?.synthesizer[0];
  }

  // todo: rename + refactor + include this.invalidCodeIdx insead of
  // indexOfFirstInvalidChar
  private compute(possibleSynthesizers: string[]) {
    const outputIndices = Array(possibleSynthesizers.length).fill(0);

    const firstUniqueCharacters = ['r', 'd'];
    for (let i = 0; i < possibleSynthesizers.length; ++i) {
      const synthesizer = possibleSynthesizers[i];
      const codes = sortByReverseLength(
        libWrapper.getCodesByFromat(synthesizer)
      );
      while (outputIndices[i] < this.sequence.length) {
        const matchedCode = codes.find((c) => c === this.sequence.slice(outputIndices[i], outputIndices[i] + c.length));

        if (!matchedCode) break;

        if ( // for mistake pattern 'rAA'
          outputIndices[i] > 1 &&
          NUCLEOTIDES.includes(this.sequence[outputIndices[i]]) &&
          firstUniqueCharacters.includes(this.sequence[outputIndices[i] - 2])
        ) break;

        if ( // for mistake pattern 'ArA'
          firstUniqueCharacters.includes(this.sequence[outputIndices[i] + 1]) &&
          NUCLEOTIDES.includes(this.sequence[outputIndices[i]])
        ) {
          outputIndices[i]++;
          break;
        }
        outputIndices[i] += matchedCode.length;
      }
    }

    const outputIndex = Math.max(...outputIndices);
    const synthesizer = possibleSynthesizers[outputIndices.indexOf(outputIndex)];
    const indexOfFirstInvalidChar = (outputIndex === this.sequence.length) ? -1 : outputIndex;
    if (indexOfFirstInvalidChar !== -1) {
      return {
        indexOfFirstInvalidChar: indexOfFirstInvalidChar,
        synthesizer: [synthesizer],
      };
    }
    return {
      indexOfFirstInvalidChar: indexOfFirstInvalidChar,
      synthesizer: [synthesizer],
    }
  }

  private getListOfPossibleSynthesizersByFirstMatchedCode(): string[] {
    const sequence = this.sequence;
    let synthesizers: string[] = [];
    // iterate over all possible formats and 
    Object.values(SYNTHESIZERS).forEach((synthesizer: string) => {
      let codes = sortByReverseLength(libWrapper.getCodesByFromat(synthesizer));
      let start = 0;
      for (let i = 0; i < sequence.length; i++) {
        if (sequence[i] === ')' && i !== sequence.length - 1) {
          start = i + 1;
          break;
        }
      }
      if (codes.some((s: string) => s === sequence.slice(start, start + s.length)))
        synthesizers.push(synthesizer);
    });
    return synthesizers;
  }
}
