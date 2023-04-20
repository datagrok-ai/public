/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MODIFICATIONS} from '../hardcode-to-be-eliminated/const';
import {map} from '../hardcode-to-be-eliminated/map';


import {SYNTHESIZERS, DELIMITER, TECHNOLOGIES} from './const';
import {isValidSequence} from './code-converter/conversion-validation-tools';

/** Auxiliary class wrapping legacy dependencies on hardcode, to be deleted after full hardcode elimination  */
export class HardcodeTerminator {
  getModificationCodes(): string[] {
    return Object.keys(MODIFICATIONS);
  }

  getModificationSmiles(code: string): {left: string, right: string} {
    return MODIFICATIONS[code];
  }

  getCodeToNameMap(sequence: string, format: string) {
    const codeToNameMap = new Map<string, string>;
    const NAME = 'name';
    if (format == null) {
      for (const synthesizer of Object.keys(map)) {
        for (const technology of Object.keys(map[synthesizer])) {
          for (const code of Object.keys(map[synthesizer][technology]))
            codeToNameMap.set(code, map[synthesizer][technology][code][NAME]!);
        }
      }
    } else {
      for (const technology of Object.keys(map[format])) {
        for (const code of Object.keys(map[format][technology]))
          codeToNameMap.set(code, map[format][technology][code][NAME]!);
      }
    }
    // what is the reason of the following condition?
    codeToNameMap.set(DELIMITER, '');
    const output = isValidSequence(sequence, format);
    const G = 'g';
    if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12)) {
      const value = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA][G][NAME]!;
      codeToNameMap.set(G, value);
    } else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS)) {
      const value = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA][G][NAME]!;
      codeToNameMap.set(G, value);
    }
    return codeToNameMap;
  }

  getSmilesByName(monomerName: string): string {
    if (Object.keys(MODIFICATIONS).includes(monomerName)) {
      return MODIFICATIONS[monomerName].left;
    } else {
      const nameToSmilesMap = this.getNameToSmilesMap();
      const smiles = nameToSmilesMap.get(monomerName);
      if (smiles === undefined)
        throw new Error(`ST: cannot find smiles for monomerName ${monomerName}`);
      return nameToSmilesMap.get(monomerName)!;
    }
  }

  get3PrimeTerminalSmiles(modificationName: string): string {
    if (!Object.keys(MODIFICATIONS).includes(modificationName))
      throw new Error(`ST: ${modificationName} is not in the list of modifications`);
    return MODIFICATIONS[modificationName].right;
  }

  private getNameToSmilesMap(): Map<string, string> {
    const nameToSmilesMap = new Map<string, string>();
    const NAME = 'name';
    const SMILES = 'SMILES';
    for (const synthesizer of Object.keys(map)) {
      for (const technology of Object.keys(map[synthesizer])) {
        for (const code of Object.keys(map[synthesizer][technology])) {
          const name = map[synthesizer][technology][code][NAME]!;
          const smiles = map[synthesizer][technology][code][SMILES]!;
          nameToSmilesMap.set(name, smiles);
        }
      }
    }
    return nameToSmilesMap;
  }
}
