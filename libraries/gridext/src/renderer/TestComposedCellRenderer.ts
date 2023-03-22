import {AbstractComposedCellRenderer} from "./AbstructComposedCellRenderer";
import * as DG from 'datagrok-api/dg';

export class TestComposedCellRenderer extends AbstractComposedCellRenderer {
  constructor(renderer: DG.GridCellRenderer) {
    super(renderer);
  }

  fillLabelsAndUI(cell: DG.GridCell, arTextLabels : string[], arTextFonts: string[], arTextColors : string[], arBackColors: string[]) : void {
    arTextLabels.push('test 1');
    arTextLabels.push('test 2');
    arTextFonts.push('Roboto 12px');
    arTextFonts.push('Roboto 12px');
    arTextColors.push('black');
    arTextColors.push('black');
    arBackColors.push('white');
    arBackColors.push('white');
  }
}
