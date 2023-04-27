/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomerLib} from '../../package';

import {IMonomerLib} from '@datagrok-libraries/bio/src/types';

import {HardcodeTerminator} from '../hardcode-terminator';
const terminator = new HardcodeTerminator();

export class MonomerLibWrapper {
  constructor() {
    const lib = getMonomerLib();
    if (lib === null)
      throw new Error('SequenceTranslator: monomer library is null');
    this.lib = lib!;
  }
  private lib: IMonomerLib;

  getMolfileByName(monomerName: string): string {
    const monomer = this.lib?.getMonomer('RNA', monomerName);
    return monomer?.molfile!;
  }

  getSmilesByName(monomerName: string): string {
    return terminator.getSmilesByName(monomerName);
  }

  get3PrimeTerminalSmiles(modificationName: string): string {
    return terminator.get3PrimeTerminalSmiles(modificationName);
  }

  isModification(monomerName: string): boolean {
    const monomerMolfile = this.getMolfileByName(monomerName);
    if (monomerMolfile !== undefined) {
      return (monomerMolfile.includes('MODIFICATION')) ? true : false;
    } else {
      const modifications = terminator.getModificationCodes();
      return (modifications.includes(monomerName)) ? true : false;
    }
  }

  getCodeToNameMap(sequence: string, format: string): Map<string, string> {
    return (new HardcodeTerminator().getCodeToNameMap(sequence, format));
  }

  getModificationCodes(): string[] {
    return terminator.getModificationCodes();
  }
}
