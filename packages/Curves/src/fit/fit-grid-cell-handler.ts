import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {
  fitSeries,
  getChartData,
  getColumnChartOptions,
  getSeriesStatistics,
  getSeriesFitFunction,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, fitSeriesProperties, fitChartDataProperties, FIT_CELL_TYPE} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {MultiCurveViewer} from './multi-curve-viewer';


function addStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property): void {
  const grid = chartColumn.grid;
  grid.dataFrame.columns
    .addNew(p.name, p.propertyType as DG.ColumnType)
    .init((i) => {
      const chartData = getChartData(
        DG.GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i)));
      const fitResult = getSeriesStatistics(chartData.series![0], getSeriesFitFunction(chartData.series![0]));
      return p.get(fitResult);
    });
}


export class FitGridCellHandler extends DG.ObjectHandler {
  get type(): string {
    return 'GridCell';
  }

  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType === FIT_CELL_TYPE;
  }

  renderProperties(gridCell: DG.GridCell, context: any = null): HTMLElement {
    const acc = ui.accordion();
    const chartData = getChartData(gridCell);
    const columnChartOptions = getColumnChartOptions(gridCell.gridColumn);
    const refresh = {onValueChanged: (_: any) => gridCell.grid.invalidate()};

    acc.addPane('Options', () => ui.divV([
        ui.input.form(columnChartOptions.seriesOptions, fitSeriesProperties, refresh),
        ui.input.form(columnChartOptions.chartOptions, fitChartDataProperties, refresh),
      ]));

    acc.addPane('Results', () => {
      const host = ui.divV([]);

      for (let i = 0; i < chartData.series!.length; i++) {
        const series = chartData.series![i];
        const fitFunction = getSeriesFitFunction(chartData.series![i]);
        if (!series.parameters)
          series.parameters = fitSeries(series, fitFunction).parameters;
        const seriesStatistics = getSeriesStatistics(series, fitFunction);
        host.appendChild(ui.panel([
          ui.h1(series.name ?? 'series ' + i),
          ui.input.form(seriesStatistics, statisticsProperties, {
            onCreated: (input) => input.root.appendChild(
              ui.iconFA('plus',
                () => addStatisticsColumn(gridCell.gridColumn, input.property),
                `Calculate ${input.property.name} for the whole column`))
          })
        ]));
      }

      return host;
    });

    acc.addPane('Chart', () => MultiCurveViewer.fromChartData(chartData).root);

    return acc.root;
  }
}
