/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HardcodeTerminator} from '../hardcode-terminator';

import {DELIMITER} from '../const';
import {LINKS, P_LINKAGE} from './const';

const terminator = new HardcodeTerminator();

export class MonomerCodeParser {
  constructor(
    private sequence: string, private invert: boolean = false, private format: string
  ) {
    this.codeToNameMap = terminator.getCodeToNameMap(this.sequence, this.format);
  }

  private codeToNameMap: Map<string, string>;

  parseSequence(): string[] {
    const parsedCodes = this.parseRawSequence();
    return this.addLinkers(parsedCodes);
  }

  private addLinkers(parsedCodes: string[]) {
    const monomerNamesWithLinks: string[] = [];
    for (let i = 0; i < parsedCodes.length; i++) {
      if (
        LINKS.includes(parsedCodes[i]) ||
        this.isMonomerAttachedToLink(parsedCodes[i]) ||
        (i < parsedCodes.length - 1 && LINKS.includes(parsedCodes[i + 1]))
      ) {
        const name = this.codeToNameMap.get(parsedCodes[i]);
        if (name !== undefined)
          monomerNamesWithLinks.push(name);
        else
          monomerNamesWithLinks.push(parsedCodes[i]);
      } else {
        const name = this.codeToNameMap.get(parsedCodes[i]);
        if (name !== undefined)
          monomerNamesWithLinks.push(name);
        else
          monomerNamesWithLinks.push(parsedCodes[i]);
        monomerNamesWithLinks.push(P_LINKAGE);
      }
    }
    return monomerNamesWithLinks;
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
    let allCodesInTheFormat = Array.from(this.codeToNameMap.keys());
    const modifications = terminator.getModificationCodes();
    allCodesInTheFormat = allCodesInTheFormat.concat(modifications).concat(DELIMITER);
    return this.reverseLengthSort(allCodesInTheFormat);
  }

  private reverseLengthSort(array: string[]): string[] {
    return array.sort((a, b) => b.length - a.length);
  }

  private isMonomerAttachedToLink(code: string) {
    // todo: eliminate this legacy list, leads to bugs
    const legacyList = ['e', 'h', /*'g',*/ 'f', 'i', 'l', 'k', 'j'];
    return legacyList.includes(code);
  }
}
