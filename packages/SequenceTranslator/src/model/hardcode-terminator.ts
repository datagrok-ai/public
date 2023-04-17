/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MODIFICATIONS} from '../hardcode-to-be-eliminated/const';
import {map} from '../hardcode-to-be-eliminated/map';


import {SYNTHESIZERS, DELIMITER, TECHNOLOGIES} from './const';
import {isValidSequence} from './code-converter/conversion-validation-tools';

/** Auxiliary class, to be deleted after full hardcode elimination  */
export class HardcodeTerminator {
  getModifications() {
    return Object.keys(MODIFICATIONS);
  }

  getCodeToNameMap(sequence: string, format: string) {
    const obj: { [code: string]: string } = {};
    const NAME = 'name';
    if (format == null) {
      for (const synthesizer of Object.keys(map)) {
        for (const technology of Object.keys(map[synthesizer])) {
          for (const code of Object.keys(map[synthesizer][technology]))
            obj[code] = map[synthesizer][technology][code][NAME]!;
        }
      }
    } else {
      for (const technology of Object.keys(map[format])) {
        for (const code of Object.keys(map[format][technology]))
          obj[code] = map[format][technology][code][NAME]!;
      }
    }
    obj[DELIMITER] = '';
    const output = isValidSequence(sequence, format);
    if (output.synthesizer!.includes(SYNTHESIZERS.MERMADE_12))
      obj['g'] = map[SYNTHESIZERS.MERMADE_12][TECHNOLOGIES.SI_RNA]['g'][NAME]!;
    else if (output.synthesizer!.includes(SYNTHESIZERS.AXOLABS))
      obj['g'] = map[SYNTHESIZERS.AXOLABS][TECHNOLOGIES.SI_RNA]['g'][NAME]!;
    return obj;
  }
}
