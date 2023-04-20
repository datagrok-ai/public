/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DELIMITER} from '../const';
import {LINKER_CODES, P_LINKAGE} from './const';
import {MonomerLibWrapper} from './monomer-handler';

/** Wrapper for parsing a strand and getting a sequence of monomer IDs (with
 * omitted linkers, if needed)  */
export class MonomerSequenceParser {
  constructor(
    private sequence: string, private invert: boolean = false,
    // todo: remove from the list of parameters
    private codeMap: Map<string, string>
  ) {
    this.lib = new MonomerLibWrapper();
  }

  private lib: MonomerLibWrapper;

  /** Get sequence of parsed monomer symbols, which are unique short names for
   * the monomers within the Monomer Library */
  parseSequence(): string[] {
    const parsedRawCodes = this.parseRawSequence();
    return this.addLinkers(parsedRawCodes);
  }

  private addLinkers(parsedRawCodes: string[]) {
    const monomerSymbolSequence: string[] = [];
    for (let i = 0; i < parsedRawCodes.length; i++) {
      const code = parsedRawCodes[i];
      const monomerSymbol = this.codeMap.get(code);
      // todo: port this validation to more appropriate place
      if (monomerSymbol === undefined) {
        console.log(this.codeMap);
        throw new Error(`SequenceTranslator: monomer with code ${code} is absent in the code map`);
      }
      monomerSymbolSequence.push(monomerSymbol);

      const isLinker = LINKER_CODES.includes(code);
      const attachedToLink = isMonomerAttachedToLink(code);
      const nextMonomerIsLinker = (i < parsedRawCodes.length - 1 && LINKER_CODES.includes(parsedRawCodes[i + 1]));

      // todo: refactor as molfile-specific
      if (!isLinker && !attachedToLink && !nextMonomerIsLinker) {
        // todo: replace by phosphate linkage ID
        monomerSymbolSequence.push(P_LINKAGE);
      }
    }
    return monomerSymbolSequence;
  }

  private parseRawSequence(): string[] {
    const allCodesOfFormat = this.getAllCodesOfFormat();
    const parsedCodes = [];
    let i = 0;
    while (i < this.sequence.length) {
      const code = allCodesOfFormat.find(
        (s: string) => s === this.sequence.substring(i, i + s.length)
      )!;
      this.invert ? parsedCodes.unshift(code) : parsedCodes.push(code);
      i += code.length;
    }
    return parsedCodes;
  }

  // todo: port to monomer handler
  private getAllCodesOfFormat(): string[] {
    let allCodesInTheFormat = Array.from(this.codeMap.keys());
    const modifications = this.lib.getModificationCodes();
    allCodesInTheFormat = allCodesInTheFormat.concat(modifications).concat(DELIMITER);
    return reverseLengthSort(allCodesInTheFormat);
  }
}
// todo: eliminate this strange legacy condition, leads to bugs
function isMonomerAttachedToLink(code: string) {
  const legacyList = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  return legacyList.includes(code);
}

function reverseLengthSort(array: string[]): string[] {
  return array.sort((a, b) => b.length - a.length);
}
