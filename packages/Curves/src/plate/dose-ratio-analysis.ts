/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from './plate';
import {FIT_FUNCTION_4PL_REGRESSION, IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';

export class PlateDoseRatioAnalysis {
  static createDoseRatioGrid(plate: Plate): HTMLElement | null {
    const requiredColumns = ['Agonist_Concentration_M', 'Antagonist_Concentration_M', 'Percent_Inhibition'];
    if (!requiredColumns.every((col) => plate.data.columns.contains(col))) {
      console.warn(`Dose Ratio Analysis: Missing one of the required columns: ${requiredColumns.join(', ')}`);
      return null;
    }

    const antagonistCol = plate.data.col('Antagonist_Concentration_M')!;
    const agonistCol = plate.data.col('Agonist_Concentration_M')!;
    const responseCol = plate.data.col('Percent_Inhibition')!;
    const uniqueConcentrations = Array.from(new Set(antagonistCol.toList()));
    const antagonistConcentrations = (uniqueConcentrations as number[]).sort((a, b) => a - b);

    const series: IFitSeries[] = [];

    for (let i = 0; i < antagonistConcentrations.length; i++) {
      const antagonistConc = antagonistConcentrations[i];
      const color = DG.Color.toHtml(DG.Color.getCategoricalColor(i));
      const currentSeries: IFitSeries = {
        name: `Antagonist: ${antagonistConc === 0 ? 'None' : `${antagonistConc.toExponential(1)} M`}`,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        points: [],
        pointColor: color,
        fitLineColor: color,
        droplines: ['IC50'],
        showPoints: 'points',
      };

      for (let j = 0; j < plate.data.rowCount; j++) {
        if (antagonistCol.get(j) === antagonistConc) {
          currentSeries.points.push({
            x: agonistCol.get(j),
            y: responseCol.get(j),
            outlier: false,
          });
        }
      }

      currentSeries.points.sort((a, b) => a.x - b.x);

      if (currentSeries.points.length > 1)
        series.push(currentSeries);
    }

    if (series.length === 0)
      return ui.divText('No data available to plot dose-ratio curves.');


    const chartData: IFitChartData = {
      chartOptions: {
        title: `Dose Ratio Analysis for ${plate.barcode}`,
        xAxisName: 'Agonist Concentration (M)',
        yAxisName: 'Percent Inhibition',
        logX: true,
        minY: -10,
        maxY: 110,
      },
      series: series,
    };

    const chartDf = DG.DataFrame.create(1);
    const curveCol = chartDf.columns.addNewString('Dose Ratio Curve');
    curveCol.init(() => JSON.stringify(chartData));
    curveCol.semType = 'fit';

    const grid = chartDf.plot.grid();

    const container = ui.div([grid.root], 'drc-grid-container');

    container.style.width = '100%';
    container.style.height = '100%';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';

    ui.tools.handleResize(container, (w: number, h: number) => {
      if (w > 20 && h > 20) {
        grid.col(curveCol.name)!.width = w - 20;
        grid.props.rowHeight = h - 20;
      }
    });

    return container;
  }
}
