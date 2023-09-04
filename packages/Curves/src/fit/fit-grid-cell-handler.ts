import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {
  fitSeries,
  getColumnChartOptions,
  getSeriesStatistics,
  getSeriesFitFunction,
  getDataFrameChartOptions,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, fitSeriesProperties, fitChartDataProperties, FIT_CELL_TYPE, TAG_FIT, IFitChartData, IFitSeries, IFitChartOptions, FIT_SEM_TYPE, IFitSeriesOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX, getChartData, mergeProperties} from './fit-renderer';
import {MultiCurveViewer} from './multi-curve-viewer';
import {convertXMLToIFitChartData} from './fit-parser';


const SOURCE_COLUMN_TAG = '.sourceColumn';
const SERIES_NUMBER_TAG = '.seriesNumber';
const STATISTICS_TAG = '.statistics';
const CHART_OPTIONS = 'chartOptions';
const SERIES_OPTIONS = 'seriesOptions';
enum MANIPULATION_LEVEL {
  DATAFRAME = 'Dataframe',
  COLUMN = 'Column',
  PLOT = 'Plot'
};


function addStatisticsColumn(chartColumn: DG.GridColumn, p: DG.Property, seriesNumber: number): void {
  const grid = chartColumn.grid;
  const column = DG.Column.float(p.name, chartColumn.column?.length);
  column.tags[SOURCE_COLUMN_TAG] = chartColumn.name;
  column.tags[SERIES_NUMBER_TAG] = seriesNumber;
  column.tags[STATISTICS_TAG] = p.name;

  column
    .init((i) => {
      const chartData = getChartData(
        DG.GridCell.fromColumnRow(grid, chartColumn.name, grid.tableRowToGrid(i)));
      const fitResult = getSeriesStatistics(chartData.series![0], getSeriesFitFunction(chartData.series![0]));
      return p.get(fitResult);
    });
  grid.dataFrame.columns.add(column);
}

function changePlotOptions(chartData: IFitChartData, inputBase: DG.InputBase, options: string): void {
  if (options === CHART_OPTIONS) {
    if (chartData.chartOptions === undefined) return;
    (chartData.chartOptions[inputBase.property.caption as keyof IFitChartOptions] as any) = inputBase.value;
  }
  else if (options === SERIES_OPTIONS) {
    if (chartData.series === undefined) return;
    for (let i = 0; i < chartData.series.length; i++)
      (chartData.series[i][inputBase.property.caption as keyof IFitSeries] as any) = inputBase.value;
  }
}

function changeColumnsCurvesOptions(columns: DG.Column[], inputBase: DG.InputBase, options: string): void {
  for (let i = 0; i < columns.length; i++) {
    for (let j = 0; j < columns[i].length; j++) {
      const value = columns[i].get(j);
      if (value === '') continue;
      const chartData: IFitChartData = columns[i].getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
        convertXMLToIFitChartData(value) : JSON.parse(value ?? '{}') ?? {};
      changePlotOptions(chartData, inputBase, options);
      columns[i].set(j, JSON.stringify(chartData));
    }
    if (columns[i].getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX)
      columns[i].setTag(TAG_FIT_CHART_FORMAT, '');
  }
}

function convertJnJColumnToJSON(column: DG.Column) {
  for (let i = 0; i < column.length; i++) {
    const value = column.get(i);
    if (value === '') continue;
    const chartData: IFitChartData = convertXMLToIFitChartData(value);
    column.set(i, JSON.stringify(chartData));
  }
  column.setTag(TAG_FIT_CHART_FORMAT, '');
}

function changeCurvesOptions(gridCell: DG.GridCell, inputBase: DG.InputBase, options: string, manipulationLevel: string): void {
  const chartOptions = manipulationLevel === MANIPULATION_LEVEL.DATAFRAME ?
    getDataFrameChartOptions(gridCell.cell.dataFrame) : getColumnChartOptions(gridCell.cell.column);
  if (options === CHART_OPTIONS)
    (chartOptions.chartOptions![inputBase.property.caption as keyof IFitChartOptions] as any) = inputBase.value;
  else if (options === SERIES_OPTIONS)
    (chartOptions.seriesOptions![inputBase.property.caption as keyof IFitSeriesOptions] as any) = inputBase.value;

  if (manipulationLevel === MANIPULATION_LEVEL.DATAFRAME) {
    gridCell.cell.dataFrame.tags[TAG_FIT] = JSON.stringify(chartOptions);
    const fitColumns = gridCell.cell.dataFrame.columns.bySemTypeAll(FIT_SEM_TYPE);
    changeColumnsCurvesOptions(fitColumns, inputBase, options);
  }
  else if (manipulationLevel === MANIPULATION_LEVEL.COLUMN) {
    gridCell.cell.column.tags[TAG_FIT] = JSON.stringify(chartOptions);
    changeColumnsCurvesOptions([gridCell.cell.column], inputBase, options);
  }
  else {
    const value = gridCell.cell.value;
    if (value === '') return;
    if (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX)
      convertJnJColumnToJSON(gridCell.cell.column);
    const chartData: IFitChartData = JSON.parse(value ?? '{}') ?? {};
    changePlotOptions(chartData, inputBase, options);
    gridCell.cell.value = JSON.stringify(chartData);
  }
  gridCell.grid.invalidate();
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
    // TODO: make just the base ui.choiceInput after nullable option is added
    const switchProperty = DG.Property.js('level', DG.TYPE.STRING, {description: 'Switch manipulation level',
      defaultValue: 'Column', choices: ['Dataframe', 'Column', 'Plot'], nullable: false});
    const switchLevelInput = DG.InputBase.forProperty(switchProperty);
    const chartData = getChartData(gridCell);
    const columnChartOptions = getColumnChartOptions(gridCell.cell.column);
    const dfChartOptions = getDataFrameChartOptions(gridCell.cell.dataFrame);

    const seriesOptionsRefresh = {onValueChanged: (inputBase: DG.InputBase) => 
      changeCurvesOptions(gridCell, inputBase, SERIES_OPTIONS, switchLevelInput.value)};
    const chartOptionsRefresh = {onValueChanged: (inputBase: DG.InputBase) =>
      changeCurvesOptions(gridCell, inputBase, CHART_OPTIONS, switchLevelInput.value)};

    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, chartData.seriesOptions ? chartData.seriesOptions :
      chartData.series ? chartData.series[0] ?? {} : {});
    mergeProperties(fitSeriesProperties, dfChartOptions.seriesOptions, chartData.seriesOptions ? chartData.seriesOptions :
      chartData.series ? chartData.series[0] ?? {} : {});
    mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions,
      chartData.chartOptions ? chartData.chartOptions : {});
    mergeProperties(fitChartDataProperties, dfChartOptions.chartOptions,
      chartData.chartOptions ? chartData.chartOptions : {});

    acc.addPane('Options', () => ui.divV([
      switchLevelInput.root,
      ui.h3('Series options'),
      ui.input.form(chartData.seriesOptions ? chartData.seriesOptions : chartData.series![0], fitSeriesProperties, seriesOptionsRefresh),
      ui.h3('Chart options'),
      ui.input.form(chartData.chartOptions, fitChartDataProperties, chartOptionsRefresh),
    ]));

    acc.addPane('Fit', () => {
      const host = ui.divV([]);

      for (let i = 0; i < chartData.series!.length; i++) {
        const series = chartData.series![i];
        const fitFunction = getSeriesFitFunction(chartData.series![i]);
        if (!series.parameters)
          series.parameters = fitSeries(series, fitFunction).parameters;
        const seriesStatistics = getSeriesStatistics(series, fitFunction);
        const color = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
          series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        host.appendChild(ui.panel([
          ui.h1(series.name ?? 'series ' + i, {style: {color: color}}),
          ui.input.form(seriesStatistics, statisticsProperties, {
            onCreated: (input) => input.root.appendChild(
              ui.iconFA('plus',
                () => addStatisticsColumn(gridCell.gridColumn, input.property, i),
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
