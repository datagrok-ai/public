import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';

export class CyclizedCellRendererBack extends MonomerPlacer {

  constructor(
    gridCol: DG.GridColumn | null, tableCol: DG.Column,
    maxLengthOfMonomer: number, seqHelper: ISeqHelper
  ) {
    super(gridCol, tableCol, _package.logger, maxLengthOfMonomer, () => {
      const sh = seqHelper.getSeqHandler(tableCol);
      return {
        seqHandler: sh,
        monomerCharWidth: 7,
        separatorWidth: 11,
        monomerToShort: monomerToShort,
      };
    });
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
    const gridCellBounds: DG.Rect = gridCell.bounds;
    const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCellBounds.x);
    const left: number | null = this.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width);
    if(left != null) {

    }
    super.onMouseMove(gridCell, e);
  }
}
