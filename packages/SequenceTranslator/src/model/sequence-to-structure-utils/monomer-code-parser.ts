/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PHOSPHATE_SYMBOL} from './const';
import {sortByReverseLength} from '../helpers';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';

/** Wrapper for parsing a strand and getting a sequence of monomer IDs (with
 * omitted linkers, if needed)  */
export class MonomerSequenceParser {
  constructor(
    private sequence: string, private invert: boolean = false,
    // todo: remove from the list of parameters
    private codeMap: Map<string, string>
  ) {
    this.lib = MonomerLibWrapper.getInstance();
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
    parsedRawCodes.forEach((code, i) => {
      const monomerSymbol = this.getSymbolForCode(code);
      monomerSymbolSequence.push(monomerSymbol);

      // todo: to be deleted
      const LINKER_CODES = ['s', 'ps', '*', 'Rpn', 'Spn', 'Rps', 'Sps']

      const isPhosphate = LINKER_CODES.includes(code);
      const hasPhosphate = isMonomerWithPhosphate(code);
      const lastMonomer = i === parsedRawCodes.length - 1;
      const nextMonomerIsPhosphate = (i + 1 < parsedRawCodes.length && LINKER_CODES.includes(parsedRawCodes[i + 1]));

      // todo: refactor as molfile-specific
      if (!isPhosphate && !hasPhosphate && !nextMonomerIsPhosphate && !lastMonomer) {
        monomerSymbolSequence.push(PHOSPHATE_SYMBOL);
      }
    });
    // console.log('monomerSymbolSequence:', monomerSymbolSequence);
    return monomerSymbolSequence;
  }

  private getSymbolForCode(code: string): string {
    let monomerSymbol = this.codeMap.get(code);
    // todo: remove as a legacy workaround, codeMap must contain all the
    // symbols, and symbols are not codes
    monomerSymbol ??= code;
    return monomerSymbol;
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
    const modifications = this.lib.getModificationGCRSCodes();
    allCodesInTheFormat = allCodesInTheFormat.concat(modifications);
    return sortByReverseLength(allCodesInTheFormat);
  }
}

// todo: eliminate this strange legacy condition, leads to bugs
function isMonomerWithPhosphate(code: string) {
  const legacyList = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
  return legacyList.includes(code);
}
