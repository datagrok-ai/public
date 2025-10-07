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
    // CHANGE a: Change type from FLOAT to INT to fix the filter UI component crash.
    {name: 'Min Value', type: DG.TYPE.INT, category: 'Output'},
    {name: 'Max Value', type: DG.TYPE.INT, category: 'Output'},
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
  // dose-ratio-analysis.ts

  // ... (inside the DoseRatioAnalysis class)

  override formatResultsForGrid(rawResults: DG.DataFrame): DG.DataFrame {
    console.log(`--- Dose Ratio Debug V2: START ---. Received ${rawResults.rowCount} raw rows.`);
    if (rawResults.rowCount === 0) {
      console.log('--- Dose Ratio Debug V2: END (no rows to process) ---');
      return rawResults;
    }

    const propsCol = rawResults.col('properties');
    if (!propsCol) {
      console.error('Dose Ratio Debug V2: CRITICAL - `properties` column not found!');
      return DG.DataFrame.create(0);
    }

    console.log('Debug V2 Step 1: Grouping rows by run_id...');
    const runs = new Map<number, { firstRow: DG.Row, seriesList: IFitSeries[] }>();
    try {
      for (let i = 0; i < rawResults.rowCount; i++) {
        const row = rawResults.row(i);
        const runId = row.get('run_id');
        const propsJson = propsCol.get(i);
        console.log(`--- Processing raw row ${i} for runId ${runId} ---`);

        if (!runs.has(runId))
          runs.set(runId, {firstRow: row, seriesList: []});

        if (!propsJson || typeof propsJson !== 'string') {
          console.warn(`Row ${i} FAILED! Properties data is not a valid string. Skipping.`);
          continue;
        }

        const properties = JSON.parse(propsJson);
        const singleSeriesChartJson = properties.Curve;
        console.log(`Row ${i}: Check 1 - 'properties.Curve' is of type [${typeof singleSeriesChartJson}]. Is it a non-empty string? ${typeof singleSeriesChartJson === 'string' && singleSeriesChartJson.length > 0}`);

        if (singleSeriesChartJson && typeof singleSeriesChartJson === 'string') {
          const chartData = JSON.parse(singleSeriesChartJson);
          const hasSeries = chartData.hasOwnProperty('series');
          const seriesIsArray = Array.isArray(chartData.series);
          console.log(`Row ${i}: Check 2 - Parsed chartData. Does it have a 'series' property? ${hasSeries}. Is 'series' an array? ${seriesIsArray}`);

          if (hasSeries && seriesIsArray && chartData.series.length > 0) {
            console.log(`Row ${i}: SUCCESS! Pushing series to list for runId ${runId}.`);
                    runs.get(runId)!.seriesList.push(chartData.series[0]);
          } else {
            console.error(`Row ${i}: FAILED! chartData.series is not a valid, non-empty array. Actual value:`, chartData.series);
          }
        } else {
          console.error(`Row ${i}: FAILED! properties.Curve was not a string or was empty.`);
        }
      }
    } catch (e) {
      console.error('Debug V2 ERROR during Step 1 (Grouping):', e);
      return DG.DataFrame.create(0);
    }

    console.log(`Debug V2 Step 1 Complete. Found ${runs.size} unique runs.`);
    // NEW: Log a summary of what was collected before proceeding
    runs.forEach((value, key) => {
      console.log(`Summary for Run ID ${key}: Collected ${value.seriesList.length} series.`);
    });


    // The rest of the function remains the same as our V2 debug version
    console.log('Debug V2 Step 2: Aggregating runs into final rows...');
    const finalRows: any[] = [];
    try {
      for (const [runId, runData] of runs.entries()) {
        if (runData.seriesList.length === 0) continue;

        let runMin = Infinity; let runMax = -Infinity;
        for (const series of runData.seriesList) {
          if (!series.points || !Array.isArray(series.points)) continue;
          for (const point of series.points) {
            if (!point.outlier && point.y != null) {
              if (point.y < runMin) runMin = point.y;
              if (point.y > runMax) runMax = point.y;
            }
          }
        }

        const combinedChartData: IFitChartData = {
          chartOptions: {
            title: `Dose Ratio for Plate ${runData.firstRow.get('barcode')}`,
            xAxisName: 'Agonist Concentration (M)', yAxisName: 'Percent Inhibition', logX: true,
          },
          series: runData.seriesList,
        };

        finalRows.push({
          'run_id': runId, 'plate_id': runData.firstRow.get('plate_id'),
          'barcode': runData.firstRow.get('barcode'), 'Curve': JSON.stringify(combinedChartData),
          'Min Value': runMin === Infinity ? null : Math.round(runMin),
          'Max Value': runMax === -Infinity ? null : Math.round(runMax),
        });
      }
    } catch (e) {
      console.error('Debug V2 ERROR during Step 2 (Aggregation):', e);
      return DG.DataFrame.create(0);
    }
    console.log(`Debug V2 Step 2 Complete. Created ${finalRows.length} final rows.`);

    console.log('Debug V2 Step 3: Creating final DataFrame.');
    if (finalRows.length === 0) {
      console.log('--- Dose Ratio Debug V2: END (no valid final rows to display) ---');
      return DG.DataFrame.create(0);
    }

    const resultDf = DG.DataFrame.fromObjects(finalRows)!;
    const curveCol = resultDf.col('Curve');
    if (curveCol) curveCol.semType = 'fit';

    console.log('--- Dose Ratio Debug V2: END (SUCCESS) ---');
    return resultDf;
  }
}
