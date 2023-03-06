import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {GridCell} from 'datagrok-api/dg';
import {divV} from 'datagrok-api/ui';

import * as fit from './fit-data';
import {fitSeries} from './fit-data';
import {fitResultProperties} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';


function addStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property): void {
  const grid = chartColumn.grid;
  grid.dataFrame.columns
    .addNew(p.name, p.propertyType as DG.ColumnType)
    .init((i) => {
      const chartData = fit.getChartData(
        GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i)));
      const fitResult = fitSeries(chartData.series![0], true);
      return p.get(fitResult);
    });
}


export class FitGridCellHandler extends DG.ObjectHandler {
  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType == 'fit';
  }

  renderProperties(gridCell: DG.GridCell, context: any = null): HTMLElement {
    const acc = ui.accordion();
    const chartData = fit.getChartData(gridCell);

    acc.addPane('Options', () => ui.input.form(chartData.chartOptions, fit.fitChartDataProperties));
    acc.addPane('Results', () => {
      const host = divV([]);

      for (let i = 0; i < chartData.series!.length; i++) {
        const series = chartData.series![i];
        const fitResult = fitSeries(series, true);
        host.appendChild(ui.panel([
          ui.h1(series.name ?? 'series ' + i),
          ui.input.form(fitResult, fitResultProperties, {
            onCreated: (input) => input.root.appendChild(
              ui.iconFA('plus',
                () => addStatisticsColumn(gridCell.gridColumn, input.property),
                `Calculate ${input.property.name} for the whole column`))
          })
        ]));
      }

      return host;
    });
    return acc.root;
  }
}
