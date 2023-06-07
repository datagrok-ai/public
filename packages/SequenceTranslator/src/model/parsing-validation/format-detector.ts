/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {sortByReverseLength} from '../helpers';
import {FORMAT} from '../const';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {SequenceValidator} from './sequence-validator';

export class FormatDetector {
  constructor (private sequence: string) {
    this.libWrapper = MonomerLibWrapper.getInstance();
  };

  private libWrapper: MonomerLibWrapper;

  getFormat(): FORMAT | null {
    // todo: reliable criterion
    if (this.sequence.startsWith('RNA'))
      return FORMAT.HELM;
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

  private getListOfPossibleSynthesizersByFirstMatchedCode(): FORMAT[] {
    const sequence = this.sequence;
    let synthesizers: FORMAT[] = [];
    for (const synthesizer of Object.values(FORMAT)) {
      let codes = sortByReverseLength(this.libWrapper.getCodesByFormat(synthesizer));
      let start = 0;
      for (let i = 0; i < sequence.length; i++) {
        if (sequence[i] === ')' && i !== sequence.length - 1) {
          start = i + 1;
          break;
        }
      }
      if (codes.some((s: string) => s === sequence.slice(start, start + s.length)))
        synthesizers.push(synthesizer);
    }
    return synthesizers;
  }
}
