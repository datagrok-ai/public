import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

/* eslint-disable max-len */
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {INotationProvider, ISeqSplitted, SeqSplittedBase, SplitterFunc} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {getSplitterWithSeparator} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GAP_SYMBOL, GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {monomerToShort, StringListSeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {Chain} from '../polytool/conversion/pt-chain';

import {_package} from '../package';
import {CyclizedCellRendererBack} from './cell-renderer-cyclized';

/* eslint-enable max-len */

export class CyclizedNotationProvider implements INotationProvider {
  private readonly separatorSplitter: SplitterFunc;
  public readonly splitter: SplitterFunc;

  get defaultGapOriginal(): string { return ''; }

  constructor(
    public readonly separator: string,
    protected readonly helmHelper: IHelmHelper
  ) {
    this.separatorSplitter = getSplitterWithSeparator(this.separator);
    this.splitter = this._splitter.bind(this);
  }

  setUnits(): void {}

  private _splitter(seq: string): ISeqSplitted {
    const baseSS: ISeqSplitted = this.separatorSplitter(seq);
    return new CyclizedSeqSplitted(
      wu.count(0).take(baseSS.length).map((p) => baseSS.getOriginal(p)).toArray(),
      GapOriginals[NOTATION.SEPARATOR]);
  }

  public getHelm(seq: string, options?: any): string {
    const seqChain = Chain.fromSeparator(seq, this.helmHelper);
    const resPseudoHelm = seqChain.getHelm();
    return resPseudoHelm;
  }

  public createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>):
  CellRendererBackBase<string> {
    const maxLengthOfMonomer: number = 4;
    // (_package.bioProperties ? _package.bioProperties.maxMonomerLength : 4) ?? 50;
    const back = new CyclizedCellRendererBack(gridCol, tableCol,
      maxLengthOfMonomer, this.helmHelper.seqHelper);

    back.init().then(() => {});
    return back;
  }
}

/** Gets canonical monomers for original ones with cyclization marks */
export class CyclizedSeqSplitted extends StringListSeqSplitted {
  override getCanonical(posIdx: number): string {
    if (this.isGap(posIdx)) return GAP_SYMBOL;

    const om = this.getOriginal(posIdx);
    let cmRes = om;
    if (om.startsWith('{'))
      cmRes = om.slice(1);
    else if (om.endsWith('}'))
      cmRes = om.slice(0, -1);
    else if (om.startsWith('('))
      cmRes = om.replace(/^\(.\d+\)/, '');
    else if (om.endsWith(')'))
      cmRes = om.replace(/\(\d+\)$/, '');
    return cmRes;
  }

  constructor(seqOList: SeqSplittedBase, gapOriginalMonomer: string) {
    super(seqOList, gapOriginalMonomer);
  }
}
