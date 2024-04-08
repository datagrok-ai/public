import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GAP_SYMBOL, INotationProvider, ISeqSplitted, SeqSplittedBase, SplitterFunc} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {getSplitterWithSeparator, StringListSeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {GapOriginals} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

export class CyclizedNotationProvider implements INotationProvider {
  private readonly separatorSplitter: SplitterFunc;
  public readonly splitter: SplitterFunc;

  constructor(
    public readonly separator: string
  ) {
    this.separatorSplitter = getSplitterWithSeparator(this.separator);
    this.splitter = this._splitter.bind(this);
  }

  private _splitter(seq: string): ISeqSplitted {
    const baseSS: ISeqSplitted = this.separatorSplitter(seq);
    return new CyclizedSeqSplitted(baseSS.originals, GapOriginals[NOTATION.SEPARATOR]);
  };
}

/** Gets canonical monomers for original ones with cyclization marks */
export class CyclizedSeqSplitted extends StringListSeqSplitted {
  private readonly seqCList: (string | null)[];

  private _canonicals: string[] | null = null;
  override get canonicals(): SeqSplittedBase {
    if (!this._canonicals) {
      const len = this.length;
      this._canonicals = new Array<string>(len);
      for (let posIdx = 0; posIdx < len; ++posIdx) {
        this._canonicals[posIdx] = this.getCanonical(posIdx);
      }
    }
    return this._canonicals;
  }

  override getCanonical(posIdx: number): string {
    if (this.isGap(posIdx)) return GAP_SYMBOL;

    let cmRes: string | null = this.seqCList[posIdx];
    if (cmRes === null) {
      const om = this.getOriginal(posIdx);
      cmRes = om;
      if (om[om.length - 1] === ')')
        cmRes = this.seqCList[posIdx] = om.replace(/\(\d+\)$/, '');
    }
    return cmRes;
  }

  constructor(seqOList: SeqSplittedBase, gapOriginalMonomer: string) {
    super(seqOList, gapOriginalMonomer);
    this.seqCList = new Array<string | null>(this.length).fill(null);
  }
}
