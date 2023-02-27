import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as fit from './fit-data';
import {FitResult} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';

export class FitGridCellHandler extends DG.ObjectHandler {
  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType == 'fit';
  }

  renderProperties(gridCell: DG.GridCell, context: any = null): HTMLElement {
    const acc = ui.accordion();
    acc.addPane('Options', () => ui.input.form(fit.getColumnChartOptions(gridCell.gridColumn),
      fit.fitChartDataProperties));
    acc.addPane('Results', () => {
      const fitResultColumn = gridCell.cell.dataFrame.col(`~fit:${gridCell.gridColumn.name}`)!;
      const fitResult: FitResult = fitResultColumn.get(gridCell.cell.row.idx);
      return ui.input.form(fitResult, fit.fitResultProperties);
    });
    return acc.root;
  }
}
