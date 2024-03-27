/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {sortByReverseLength} from '../helpers';
import {DEFAULT_FORMATS} from '../const';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {CODES_TO_HELM_DICT} from '../data-loader/json-loader';
import {SequenceValidator} from './sequence-validator';

export class FormatDetector {
  constructor(private sequence: string) {
    this.libWrapper = MonomerLibWrapper.getInstance();
    this.formats = Object.keys(CODES_TO_HELM_DICT);
  };

  private libWrapper: MonomerLibWrapper;
  private formats: string[];

  getFormat(): string | null {
    // todo: reliable criterion
    if (this.sequence.startsWith('RNA'))
      return DEFAULT_FORMATS.HELM;
    const possibleFormats = this.getListOfPossibleSynthesizersByFirstMatchedCode();
    if (possibleFormats.length === 0)
      return null;

    const validator = new SequenceValidator(this.sequence);
    const outputIndices = Array(possibleFormats.length).fill(0);
    for (let i = 0; i < possibleFormats.length; ++i) {
      const format = possibleFormats[i];
      outputIndices[i] = validator.getInvalidCodeIndex(format);
    }
    const formatIdx = (outputIndices.some((idx) => idx === -1)) ? -1 : Math.max(...outputIndices);
    return possibleFormats[outputIndices.indexOf(formatIdx)];
  }

  // todo: rename
  private getListOfPossibleSynthesizersByFirstMatchedCode(): string[] {
    const sequence = this.sequence;
    const synthesizers: string[] = [];
    for (const format of this.formats) {
      const codes = sortByReverseLength(this.libWrapper.getCodesByFormat(format));
      let start = 0;
      for (let i = 0; i < sequence.length; i++) {
        if (sequence[i] === ')' && i !== sequence.length - 1) {
          start = i + 1;
          break;
        }
      }
      if (codes.some((s: string) => s === sequence.slice(start, start + s.length)))
        synthesizers.push(format);
    }
    return synthesizers;
  }
}
