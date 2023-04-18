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
  getModificationCodes() {
    return Object.keys(MODIFICATIONS);
  }

  getModificationSmiles(code: string): {left: string, right: string} {
    return MODIFICATIONS[code];
  }

  // todo: port to helmLib
  getPhosphateSmiles(): string {
    return 'OP(=O)(O)O';
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

  getCodeToSmilesMap(sequence: string, format: string) {
    const codeToSmilesMap = new Map<string, string>();
    const SMILES = 'SMILES';
    if (format == null) {
      for (const synthesizer of Object.keys(map)) {
        for (const technology of Object.keys(map[synthesizer])) {
          for (const code of Object.keys(map[synthesizer][technology]))
            codeToSmilesMap.set(code, map[synthesizer][technology][code].SMILES);
        }
      }
    } else {
      for (const technology of Object.keys(map[format])) {
        for (const code of Object.keys(map[format][technology]))
          codeToSmilesMap.set(code, map[format][technology][code].SMILES);
      }
    }
    codeToSmilesMap.set(DELIMITER, '');
    const output = isValidSequence(sequence, format);
    const G = 'g';
    if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12)) {
      const value = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA][G][SMILES]!;
      codeToSmilesMap.set(G, value);
    } else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS)) {
      const value = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA][G][SMILES]!;
      codeToSmilesMap.set(G, value);
    }
    return codeToSmilesMap;
  }
}
