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
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {Chain} from '../polytool/pt-conversion';

import {_package} from '../package';

/* eslint-enable max-len */

export class CyclizedNotationProvider implements INotationProvider {
  private readonly separatorSplitter: SplitterFunc;
  public readonly splitter: SplitterFunc;

  constructor(
    public readonly separator: string,
    protected readonly seqHelper: ISeqHelper
  ) {
    this.separatorSplitter = getSplitterWithSeparator(this.separator);
    this.splitter = this._splitter.bind(this);
  }

  private _splitter(seq: string): ISeqSplitted {
    const baseSS: ISeqSplitted = this.separatorSplitter(seq);
    return new CyclizedSeqSplitted(
      wu.count(0).take(baseSS.length).map((p) => baseSS.getOriginal(p)).toArray(),
      GapOriginals[NOTATION.SEPARATOR]);
  }

  public async getHelm(seq: string, options?: any): Promise<string> {
    const helmHelper = await getHelmHelper();
    const seqChain = await Chain.parseNotation(seq, helmHelper);
    const resPseudoHelm = seqChain.getNotationHelm();
    return resPseudoHelm;
  }

  public createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string> {
    let maxLengthOfMonomer: number = 4; // (_package.bioProperties ? _package.bioProperties.maxMonomerLength : 4) ?? 50;
    const back = new MonomerPlacer(gridCol, tableCol, _package.logger, maxLengthOfMonomer,
      () => {
        const sh = this.seqHelper.getSeqHandler(tableCol);
        return {
          seqHandler: sh,
          monomerCharWidth: 7,
          separatorWidth: 11,
          monomerToShort: monomerToShort,
        };
      });
    back.init().then(() => {});
    return back;
  }
}

/** Gets canonical monomers for original ones with cyclization marks */
export class CyclizedSeqSplitted extends StringListSeqSplitted {
  private readonly seqCList: (string | null)[];

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
