/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HardcodeTerminator} from '../hardcode-terminator';

import {DELIMITER} from '../const';
import {sortByStringLengthInDescendingOrder} from '../helpers';
import {LINKER_CODES} from './const';

const terminator = new HardcodeTerminator();

// todo: remove intersections with sequence-to-molfile
export class SequenceToSmilesConverter {
  constructor(
    private sequence: string, private invert: boolean = false, private format: string
  ) {
    this.codeToSmilesMap = terminator.getCodeToSmilesMap(this.sequence, this.format);
  };

  private codeToSmilesMap: Map<string, string>;

  convert(): string {
    const parsedRawCodes = this.parseRawSequence();
    const modifications = terminator.getModificationCodes();
    let resultingSmiles = '';
    for (let i = 0; i < parsedRawCodes.length; i++) {
      const code = parsedRawCodes[i];
      const phosphate = terminator.getPhosphateSmiles();
      if (modifications.includes(code)) {
        const modificationSmiles = (i > parsedRawCodes.length / 2) ?
          terminator.getModificationSmiles(code).right :
          terminator.getModificationSmiles(code).left;
        resultingSmiles += modificationSmiles + phosphate;
      } else {
        if (
          LINKER_CODES.includes(code) ||
          this.isMonomerAttachedToLink(code) ||
          (i < parsedRawCodes.length - 1 && LINKER_CODES.includes(parsedRawCodes[i + 1]))
        )
          resultingSmiles += this.codeToSmilesMap.get(code);
        else
          resultingSmiles += this.codeToSmilesMap.get(code) + terminator.getPhosphateSmiles();
      }
    }
    resultingSmiles = resultingSmiles.replace(/OO/g, 'O');

    const len = parsedRawCodes.length;
    const vagueLegacyPredicate = (LINKER_CODES.includes(parsedRawCodes[len - 1]) &&
      len > 1 && !this.isMonomerAttachedToLink(parsedRawCodes[len - 2])) ||
      modifications.includes(parsedRawCodes[len - 1]) ||
      this.isMonomerAttachedToLink(parsedRawCodes[len - 1]);

    return vagueLegacyPredicate ? resultingSmiles :
      resultingSmiles.slice(0, resultingSmiles.length - terminator.getPhosphateSmiles().length + 1);
  }

  private isMonomerAttachedToLink(code: string) {
    // todo: eliminate this legacy list, leads to bugs
    const legacyList = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
    return legacyList.includes(code);
  }

  private parseRawSequence(): string[] {
    const allCodesOfFormat = this.getAllCodesOfFormat();
    const parsedCodes = [];
    let i = 0;
    while (i < this.sequence.length) {
      const code = allCodesOfFormat.find(
        (s: string) => s === this.sequence.slice(i, i + s.length)
      )!;
      this.invert ? parsedCodes.unshift(code) : parsedCodes.push(code);
      i += code.length;
    }
    return parsedCodes;
  }

  private getAllCodesOfFormat(): string[] {
    let allCodesInTheFormat = Array.from(this.codeToSmilesMap.keys());
    const modifications = terminator.getModificationCodes();
    allCodesInTheFormat = allCodesInTheFormat.concat(modifications).concat(DELIMITER);
    return sortByStringLengthInDescendingOrder(allCodesInTheFormat);
  }
}
