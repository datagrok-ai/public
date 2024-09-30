import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {INotationProvider, ISeqSplitted, SeqSplittedBase, SplitterFunc}
  from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {getSplitterWithSeparator, StringListSeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {GAP_SYMBOL, GapOriginals} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {monomerToShortFunction} from './cell-renderer';

import {_package} from '../package';

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
    return new CyclizedSeqSplitted(
      wu.count(0).take(baseSS.length).map((p) => baseSS.getOriginal(p)).toArray(),
      GapOriginals[NOTATION.SEPARATOR]);
  }

  public async getHelm(seqCol: DG.Column<string>, options?: any): Promise<DG.Column<string>> {
    const polyToolPackageName: string = 'SequenceTranslator';

    const funcList = DG.Func.find({package: polyToolPackageName, name: 'polyToolConvert2'});
    if (funcList.length == 0)
      throw new Error(`Package '${polyToolPackageName}' must be installed for Cyclized notation provider.`);
    const func = funcList[0];

    const ptConvertCall = await func.prepare({table: seqCol.dataFrame, seqCol: seqCol, ...options});

    const editorFunc = DG.Func.find({package: polyToolPackageName, name: 'getPolyToolConvertEditor'})[0];
    const resHelmCol = (await editorFunc.prepare({call: ptConvertCall}).call()).getOutputParamValue() as DG.Column<string>;
    return resHelmCol;
  }

  public createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string> {
    let maxLengthOfMonomer: number = (_package.properties ? _package.properties.maxMonomerLength : 4) ?? 50;
    return new MonomerPlacer(gridCol, tableCol, _package.logger, maxLengthOfMonomer,
      () => {
        const sh = SeqHandler.forColumn(tableCol);
        return {
          seqHandler: sh,
          monomerCharWidth: 7,
          separatorWidth: 11,
          monomerToShort: monomerToShortFunction,
        };
      });
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
