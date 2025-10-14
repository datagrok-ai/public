import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';

export abstract class CellRendererWithMonomerLibBackBase extends CellRendererBackBase<string> {
  protected monomerLib: IMonomerLibBase | null = null;
  protected seqHelper: ISeqHelper | null = null;

  protected constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
  ) {
    super(gridCol, tableCol, _package.logger);

    getMonomerLibHelper().then((libHelper) => {
      this.monomerLib = libHelper.getMonomerLib();
      this.dirty = true;
      this.gridCol?.grid?.invalidate();
      this.subs.push(this.monomerLib.onChanged.subscribe(() => {
        this.dirty = true;
        this.gridCol?.grid?.invalidate();
      }));
    });
    getSeqHelper().then((seqHelper) => {
      this.seqHelper = seqHelper;
      this.dirty = true;
      this.gridCol?.grid?.invalidate();
    });
  }
}
