/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerSequenceParser} from './monomer-code-parser';
import {MonomerLibWrapper} from '../monomer-lib-utils/lib-wrapper';

// todo: merge with sequence-to-molfile
export class SequenceToSmilesConverter {
  constructor(
    sequence: string, invert: boolean = false, format: string
  ) {
    this.lib = MonomerLibWrapper.getInstance();
    const codeToNameMap = this.lib.getCodeToNameMap(format);
    this.parser = new MonomerSequenceParser(sequence, invert, codeToNameMap);
  };

  private parser: MonomerSequenceParser;
  private lib: MonomerLibWrapper;

  convert(): string {
    const parsedSequence = this.parser.parseSequence();
    const monomerSmilesArray: string[] = [];
    for (let idx = 0; idx < parsedSequence.length; idx++) {
      const monomerSymbol = parsedSequence[idx];
      const monomerSmiles = this.getMonomerSmiles(monomerSymbol, idx);
      monomerSmilesArray.push(monomerSmiles);
    }
    console.log('monomer Smiles array:', monomerSmilesArray);
    return this.getPolymerSmiles(monomerSmilesArray);
  }

  private getMonomerSmiles(monomerName: string, idx: number): string {
    if (this.lib.isModification(monomerName) && idx > 0)
      return this.lib.get3PrimeTerminalSmiles(monomerName);
    else
      return this.lib.getSmilesByName(monomerName);
  }

  private getPolymerSmiles(monomerSmiles: string[]): string {
    let resultingSmiles = monomerSmiles.join('');
    resultingSmiles = resultingSmiles.replace(/OO/g, 'O');
    return resultingSmiles;
  }

  // convert_old(): string {
  //   const parsedRawCodes = this.parseRawSequence();
  //   const modifications = terminator.getModificationCodes();
  //   let resultingSmiles = '';
  //   for (let i = 0; i < parsedRawCodes.length; i++) {
  //     const code = parsedRawCodes[i];
  //     const phosphate = terminator.getPhosphateSmiles();
  //     if (modifications.includes(code)) {
  //       const modificationSmiles = (i > parsedRawCodes.length / 2) ?
  //         terminator.getModificationSmiles(code).right :
  //         terminator.getModificationSmiles(code).left;
  //       resultingSmiles += modificationSmiles + phosphate;
  //     } else {
  //       if (
  //         LINKER_CODES.includes(code) ||
  //         this.isMonomerAttachedToLink(code) ||
  //         (i < parsedRawCodes.length - 1 && LINKER_CODES.includes(parsedRawCodes[i + 1]))
  //       )
  //         resultingSmiles += this.codeToSmilesMap.get(code);
  //       else
  //         resultingSmiles += this.codeToSmilesMap.get(code) + terminator.getPhosphateSmiles();
  //     }
  //   }
  //   resultingSmiles = resultingSmiles.replace(/OO/g, 'O');

  //   const len = parsedRawCodes.length;
  //   // seemingly, a condition for the trailing phosphate
  //   const vagueLegacyPredicate = (LINKER_CODES.includes(parsedRawCodes[len - 1]) &&
  //     len > 1 && !this.isMonomerAttachedToLink(parsedRawCodes[len - 2])) ||
  //     modifications.includes(parsedRawCodes[len - 1]) ||
  //     this.isMonomerAttachedToLink(parsedRawCodes[len - 1]);

  //   // cutting of the trailing phosphate if necessary
  //   return vagueLegacyPredicate ? resultingSmiles :
  //     resultingSmiles.slice(0, resultingSmiles.length - terminator.getPhosphateSmiles().length + 1);
  // }
}
