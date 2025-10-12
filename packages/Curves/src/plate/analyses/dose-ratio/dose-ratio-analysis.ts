/* eslint-disable max-len */
// src/plate/analyses/dose-ratio/dose-ratio-analysis.ts

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Plate} from '../../plate';
import {AnalysisBase, IAnalysisProperty} from '../base-analysis';
import {AnalysisRequiredFields} from '../../../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {FIT_FUNCTION_4PL_REGRESSION, IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {createAnalysisRun, saveAnalysisResult, saveAnalysisRunParameter} from '../../../plates/plates-crud';
import './../plate-analyses.css';

export class DoseRatioAnalysis extends AnalysisBase {
  readonly name: string = 'Dose-Ratio';
  readonly friendlyName: string = 'Dose Ratio';

  private _lastRunSeries: IFitSeries[] = [];
  private _lastRunConcentrations: number[] = [];

  parameters: IAnalysisProperty[] = [];

  outputs: IAnalysisProperty[] = [
    {name: 'Curve', type: DG.TYPE.STRING, category: 'Output'},
    {name: 'Series Name', type: DG.TYPE.STRING, category: 'Output'},
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

    const saveButton = ui.button('SAVE RESULTS', async () => {
      ui.setUpdateIndicator(saveButton, true);
      try {
        await this.saveResults(plate, chartDf, {}, mappings);
      } catch (e) {
      } finally {
        ui.setUpdateIndicator(saveButton, false);
      }
    });
    
    // Use the new scoped class name and remove inline styles
    const container = ui.divV([grid.root, ui.div([saveButton], 'ui-box')], 'assay-plates--analysis-grid-container');

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
    if (rawResults.rowCount === 0)
      return rawResults;

    const propsCol = rawResults.col('properties');
    if (!propsCol) {
      console.error('Dose Ratio: CRITICAL - `properties` column not found!');
      return DG.DataFrame.create(0);
    }

    const finalRows: any[] = [];
    for (let i = 0; i < rawResults.rowCount; i++) {
      const row = rawResults.row(i);
      const propsJson = propsCol.get(i);
      if (!propsJson) continue;

      try {
        const properties = JSON.parse(propsJson);
        const groupArray = row.get('group_combination');

        finalRows.push({
          'run_id': row.get('run_id'),
          'plate_id': row.get('plate_id'),
          'barcode': row.get('barcode'),
          'group_combination': groupArray,
          'Antagonist Conc.': Array.isArray(groupArray) && groupArray.length > 0 ?
            parseFloat(groupArray[0]).toExponential(1) + ' M' : null,
          'Curve': properties.Curve,
          'Min Value': properties['Min Value'],
          'Max Value': properties['Max Value'],
        });
      } catch (e) {
        console.error(`Failed to parse Dose Ratio properties JSON for row ${i}:`, propsJson, e);
      }
    }

    if (finalRows.length === 0)
      return DG.DataFrame.create(0);

    const resultDf = DG.DataFrame.fromObjects(finalRows)!;

    const curveCol = resultDf.col('Curve');
    if (curveCol)
      curveCol.semType = 'fit';

    return resultDf;
  }
}
