/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HardcodeTerminator} from '../hardcode-terminator';

import {DELIMITER} from '../const';
import {LINKER_CODES, P_LINKAGE} from './const';

const terminator = new HardcodeTerminator();

/** Wrapper for parsing a strand and getting a sequence of monomer IDs (with
 * omitted linkers, if needed)  */
export class MonomerCodeParser {
  constructor(
    private sequence: string, private invert: boolean = false,
    private codeMap: Map<string, string>
  ) { }

  parseSequence(): string[] {
    const parsedCodes = this.parseRawSequence();
    return this.addLinkers(parsedCodes);
  }

  private addLinkers(parsedRawCodes: string[]) {
    const monomerIdSequence: string[] = [];
    for (let i = 0; i < parsedRawCodes.length; i++) {
      const code = parsedRawCodes[i];
      const monomerId = this.codeMap.get(code);
      // todo: port this validation to more appropriate place
      if (monomerId === undefined)
        throw new Error('SequenceTranslator: monomer is absent in the code map');
      monomerIdSequence.push(monomerId);

      const isLinker = LINKER_CODES.includes(code);
      const attachedToLink = this.isMonomerAttachedToLink(code);
      const nextMonomerIsLinker = (i < parsedRawCodes.length - 1 && LINKER_CODES.includes(parsedRawCodes[i + 1]));

      // todo: refactor as molfile-specific
      if (!isLinker && !attachedToLink && !nextMonomerIsLinker) {
        // todo: replace by phosphate linkage ID
        monomerIdSequence.push(P_LINKAGE);
      }
    }
    return monomerIdSequence;
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
    let allCodesInTheFormat = Array.from(this.codeMap.keys());
    const modifications = terminator.getModificationCodes();
    allCodesInTheFormat = allCodesInTheFormat.concat(modifications).concat(DELIMITER);
    return reverseLengthSort(allCodesInTheFormat);
  }

  // todo: eliminate this strange legacy condition, leads to bugs
  private isMonomerAttachedToLink(code: string) {
    const legacyList = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
    return legacyList.includes(code);
  }
}

function reverseLengthSort(array: string[]): string[] {
  return array.sort((a, b) => b.length - a.length);
}
