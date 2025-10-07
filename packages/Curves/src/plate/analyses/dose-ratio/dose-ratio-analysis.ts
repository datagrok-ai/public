/* eslint-disable max-len */
// src/plate/analyses/dose-ratio/dose-ratio-analysis.ts

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Plate} from '../../plate';
import {AbstractPlateAnalysis, IAnalysisProperty} from '../base-analysis';
import {AnalysisRequiredFields} from '../../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {FIT_FUNCTION_4PL_REGRESSION, IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {createAnalysisRun, saveAnalysisResult, saveAnalysisRunParameter} from '../../../plates/plates-crud';

export class DoseRatioAnalysis extends AbstractPlateAnalysis {
  readonly name: string = 'Dose-Ratio';
  readonly friendlyName: string = 'Dose Ratio';

  private _lastRunSeries: IFitSeries[] = [];
  private _lastRunConcentrations: number[] = [];

  parameters: IAnalysisProperty[] = [];

  outputs: IAnalysisProperty[] = [
    {name: 'Curve', type: DG.TYPE.STRING, category: 'Output'},
    {name: 'Series Name', type: DG.TYPE.STRING, category: 'Output'},
    {name: 'Min Value', type: DG.TYPE.FLOAT, category: 'Output'},
    {name: 'Max Value', type: DG.TYPE.FLOAT, category: 'Output'},
  ];

  getRequiredFields(): AnalysisRequiredFields[] {
    return [
      {name: 'Agonist_Concentration_M', required: true, description: 'Agonist concentration values'},
      {name: 'Antagonist_Concentration_M', required: true, description: 'Antagonist concentration values'},
      {name: 'Percent_Inhibition', required: true, description: 'Response values as percent inhibition'},
      {name: 'Agonist_ID', required: false, description: 'Optional: Agonist identifier'},
      {name: 'Antagonist_ID', required: false, description: 'Optional: Antagonist identifier'},
    ];
  }

  async saveResults(
    plate: Plate,
    resultsDf: DG.DataFrame,
    params: Record<string, any>,
    currentMappings: Map<string, string>
  ): Promise<void> {
    if (!plate.id) {
      grok.shell.warning('Please save the plate first using the "CREATE" button before saving analysis results.');
      return;
    }
    if (this._lastRunSeries.length === 0) {
      grok.shell.error('No analysis results available to save.');
      return;
    }

    const pi = DG.TaskBarProgressIndicator.create('Saving Dose-Ratio results...');
    try {
      const groups = this._lastRunConcentrations.map(String);
      const runId = await createAnalysisRun(plate.id, this.name, groups);

      for (const [requiredField, actualColumn] of currentMappings.entries()) {
        const mappingPropName = `Mapping: ${requiredField}`;
        await saveAnalysisRunParameter({
          runId: runId,
          propertyName: mappingPropName,
          propertyType: DG.TYPE.STRING,
          value: actualColumn,
        });
      }

      for (let i = 0; i < this._lastRunSeries.length; i++) {
        const currentSeries = this._lastRunSeries[i];
        const antagonistConc = this._lastRunConcentrations[i];
        const groupCombination = [String(antagonistConc)];

        const yValues = (currentSeries.points || []).filter((p) => !p.outlier).map((p) => p.y!);
        const minValue = yValues.length > 0 ? Math.min(...yValues) : null;
        const maxValue = yValues.length > 0 ? Math.max(...yValues) : null;

        const chartData: IFitChartData = {
          chartOptions: {
            title: currentSeries.name,
            xAxisName: 'Agonist Concentration (M)',
            yAxisName: 'Percent Inhibition',
            logX: true,
          },
          series: [currentSeries],
        };

        const singleSeriesJson = JSON.stringify(chartData);

        await this.saveAnalysisProperty(runId, 'Curve', singleSeriesJson, groupCombination);
        await this.saveAnalysisProperty(runId, 'Series Name', currentSeries.name, groupCombination);
        if (minValue !== null)
          await this.saveAnalysisProperty(runId, 'Min Value', minValue, groupCombination);
        if (maxValue !== null)
          await this.saveAnalysisProperty(runId, 'Max Value', maxValue, groupCombination);
      }
      grok.shell.info(`Saved ${this._lastRunSeries.length} dose-ratio curves for plate ${plate.barcode}.`);
    } catch (e: any) {
      grok.shell.error(`Failed to save dose-ratio results: ${e.message}`);
      console.error(e);
    } finally {
      pi.close();
    }
  }


  private async saveAnalysisProperty(runId: number, propName: string, value: any, group: string[]): Promise<void> {
    const prop = this._outputProperties.get(propName);
    if (prop && value !== null && value !== undefined)
      await saveAnalysisResult({runId, propertyId: prop.id, propertyName: prop.name, propertyType: prop.type, value, groupCombination: group});
  }

  createView(plate: Plate, plateWidget: DG.Widget, currentMappings: Map<string, string>, onMap: (t: string, s: string) => void, onUndo: (t: string) => void, onRerender?:()=>void): HTMLElement {
    return this._createStandardMappingView(plate, currentMappings, onMap, onUndo, (p, m) => {
      return this._createResultsView(p, m);
    });
  }


  private _createResultsView(plate: Plate, mappings: Map<string, string>): HTMLElement | null {
    const {series, antagonistConcentrations} = this._runAnalysis(plate, mappings);
    this._lastRunSeries = series;
    this._lastRunConcentrations = antagonistConcentrations;

    if (series.length === 0)
      return ui.divText('No data available to plot dose-ratio curves.');

    const chartData: IFitChartData = {
      chartOptions: {
        title: `Dose Ratio Analysis for ${plate.barcode}`,
        xAxisName: 'Agonist Concentration (M)', yAxisName: 'Percent Inhibition',
        logX: true, minY: -10, maxY: 110,
      },
      series: series,
    };

    const chartDf = DG.DataFrame.fromObjects([{'curve': JSON.stringify(chartData)}])!;
    chartDf.col('curve')!.semType = 'fit';
    const grid = chartDf.plot.grid();

    // The call here needs to be updated to pass the mappings
    const saveButton = ui.button('SAVE RESULTS', async () => {
      await this.saveResults(plate, chartDf, {}, mappings); // <-- Pass mappings here
    });
    const container = ui.divV([grid.root, ui.div([saveButton], 'ui-box')], 'drc-grid-container');

    container.style.cssText = 'display: flex; flex-direction: column; width: 100%; height: 100%;';
    grid.root.style.flexGrow = '1';
    ui.tools.handleResize(container, (w: number, h: number) => {
      if (grid.col('curve')) {
        grid.col('curve')!.width = w - 20;
        grid.props.rowHeight = h - 60;
      }
    });

    return container;
  }

  private _runAnalysis(plate: Plate, mappings: Map<string, string>): {series: IFitSeries[], antagonistConcentrations: number[]} {
    const agonistColName = mappings.get('Agonist_Concentration_M')!;
    const antagonistColName = mappings.get('Antagonist_Concentration_M')!;
    const responseColName = mappings.get('Percent_Inhibition')!;
    const antagonistCol = plate.data.col(antagonistColName)!;
    const uniqueConcs = Array.from(new Set(antagonistCol.toList())).filter((c) => c !== null) as number[];
    const antagonistConcentrations = uniqueConcs.sort((a, b) => a - b);
    const series: IFitSeries[] = [];
    for (let i = 0; i < antagonistConcentrations.length; i++) {
      const antagConc = antagonistConcentrations[i];
      const color = DG.Color.toHtml(DG.Color.getCategoricalColor(i));
      const currentSeries: IFitSeries = {
        name: `Antagonist: ${antagConc === 0 ? 'None' : `${antagConc.toExponential(1)} M`}`,
        fitFunction: FIT_FUNCTION_4PL_REGRESSION,
        points: [],
        pointColor: color, fitLineColor: color,
        droplines: ['IC50'], showPoints: 'points',
      };
      for (let j = 0; j < plate.data.rowCount; j++) {
        if (antagonistCol.get(j) === antagConc) {
          currentSeries.points.push({
            x: plate.data.get(agonistColName, j),
            y: plate.data.get(responseColName, j),
          });
        }
      }
      currentSeries.points.sort((a, b) => (a.x || 0) - (b.x || 0));
      if (currentSeries.points.length > 1)
        series.push(currentSeries);
    }
    return {series, antagonistConcentrations};
  }

  protected _getGroups(resultsDf: DG.DataFrame): { groupColumn: string; groups: string[]; } {
    throw new Error('Method not implemented for DoseRatioAnalysis.');
  }
  override formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame {
    console.log('--- Dose Ratio: formatResultsForGrid START (Final Version) ---');
    if (rawResults.rowCount === 0) return rawResults;

    const runs = new Map<number, any[]>();
    const propsCol = rawResults.col('properties');

    if (!propsCol) {
      console.error('Dose Ratio: CRITICAL - `properties` column not found!');
      return DG.DataFrame.create(0);
    }

    // Group rows by run_id, extracting clean data from the 'properties' JSON
    for (let i = 0; i < rawResults.rowCount; i++) {
      const runId = rawResults.get('run_id', i);
      const propsJson = propsCol.get(i);
      if (!propsJson) continue;

      if (!runs.has(runId))
        runs.set(runId, []);

      const properties = JSON.parse(propsJson);

        runs.get(runId)!.push({
          run_id: runId,
          plate_id: rawResults.get('plate_id', i),
          barcode: rawResults.get('barcode', i),
          Curve: properties.Curve, // This is already a JSON string of a single-series chart
        });
    }


    const finalRows: any[] = [];

    for (const [runId, runRows] of runs.entries()) {
      const allSeries: IFitSeries[] = [];
      const firstRow = runRows[0]; // For common info like barcode

      runRows.forEach((rowObj) => {
        // The 'Curve' column contains a stringified IFitChartData with a *single* series.
        // We need to extract that single series object.
        const curveJson = rowObj.Curve;
        if (curveJson && typeof curveJson === 'string') {
          try {
            const chartData: IFitChartData = JSON.parse(curveJson);
            if (chartData.series && chartData.series.length > 0)
              allSeries.push(chartData.series[0]);
          } catch (e) {
            console.error(`Failed to parse curve JSON for run ${runId}:`, curveJson, e);
          }
        }
      });

      if (allSeries.length > 0) {
        // Construct the new, aggregated chart data containing all series for the run
        const combinedChartData: IFitChartData = {
          chartOptions: {
            title: `Dose Ratio for Plate ${firstRow.barcode}`,
            xAxisName: 'Agonist Concentration (M)',
            yAxisName: 'Percent Inhibition',
            logX: true,
          },
          series: allSeries,
        };

        // Create the final aggregated row for our new DataFrame
        finalRows.push({
          'run_id': runId,
          'plate_id': firstRow.plate_id,
          'barcode': firstRow.barcode,
          'Curve': JSON.stringify(combinedChartData),
        });
      }
    }

    if (finalRows.length === 0)
      return DG.DataFrame.create(0);

    const resultDf = DG.DataFrame.fromObjects(finalRows)!;
    const curveCol = resultDf.col('Curve');
    if (curveCol)
      curveCol.semType = 'fit';
    console.log('--- Dose Ratio: formatResultsForGrid END ---');
    console.log('Final Dose Ratio DataFrame structure:\n', resultDf.toCsv());

    return resultDf;
  }
}
