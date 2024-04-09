import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getSplitterWithSeparator} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {
  INotationProvider, ISeqSplitted, SeqSplittedBase, SplitterFunc
} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

//TODO: This is the rough draft yet
export class DimerizedNotationProvider implements INotationProvider {
  private readonly separatorSplitter: SplitterFunc;
  public readonly splitter: SplitterFunc;

  constructor(
    public readonly separator: string
  ) {
    this.separatorSplitter = getSplitterWithSeparator(this.separator);
    this.splitter = this._splitter.bind(this);
  }

  private _splitter(seq: string): ISeqSplitted {
    return new DimerizedSeqSplitted(seq);
  }
}

export class DimerizedSeqSplitted implements ISeqSplitted {
  private mList: string[];

  constructor(seq: string) {
    this.mList = seq.split('-');
  }

  get length(): number { return this.mList.length;}

  get canonicals(): SeqSplittedBase { return this.mList; }

  getCanonical(posIdx: number): string { return this.mList[posIdx]; }

  getOriginal(posIdx: number): string { return this.mList[posIdx]; }

  isGap(posIdx: number): boolean { return this.getOriginal(posIdx) === ''; }

  get originals(): SeqSplittedBase { return this.mList; }
}