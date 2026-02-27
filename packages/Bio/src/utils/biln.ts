/* eslint-disable max-len */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* eslint-disable max-len */
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {INotationProvider, NotationProviderBase, SplitterFunc} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {monomerToShort, splitterAsBiln} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {_package} from '../package';
/* eslint-enable max-len */

export class BilnNotationProvider extends NotationProviderBase implements INotationProvider {
  public readonly splitter: SplitterFunc;

  get defaultGapOriginal(): string { return ''; }

  static override get notationName(): string { return NOTATION.BILN; }

  static override get implementsFromHelm(): boolean { return false; }

  static override convertFromHelm(helm: string, options: any): string {
    throw new Error('Canonical way of converting from helm to biln must be used');
  }
  constructor(
    public readonly separator: string,
    public readonly seqHelper: ISeqHelper,
    public readonly seqCol: DG.Column
  ) {
    super();
    this.splitter = splitterAsBiln.bind(this);
  }

  setUnits(): void {}

  public getHelm(seq: string, _options?: any): string {
    // return resPseudoHelm;
    // generate helm from biln
    const seqSplitted = this.splitter(seq);
    const sh = this.seqHelper.getSeqHandler(this.seqCol);
    return sh.getJoiner({notation: NOTATION.HELM})(seqSplitted);
  }

  public createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>):
  CellRendererBackBase<string> {
    const maxLengthOfMonomer = _package.properties.maxMonomerLength || 4;
    // (_package.bioProperties ? _package.bioProperties.maxMonomerLength : 4) ?? 50;
    const back = new BilnCellRendererBack(gridCol, tableCol,
      maxLengthOfMonomer, this.seqHelper);

    back.init().then(() => {});
    return back;
  }
}

export class BilnCellRendererBack extends MonomerPlacer {
  constructor(
    gridCol: DG.GridColumn | null, tableCol: DG.Column,
    maxLengthOfMonomer: number, seqHelper: ISeqHelper
  ) {
    super(gridCol, tableCol, _package.logger, maxLengthOfMonomer, () => {
      const sh = seqHelper.getSeqHandler(tableCol);
      const {font, fontWidth} = MonomerPlacer.getFontSettings(tableCol);
      return {
        seqHandler: sh,
        font: font,
        fontCharWidth: fontWidth,
        separatorWidth: 0,
        monomerToShort: monomerToShort,
      };
    });
  }
}
