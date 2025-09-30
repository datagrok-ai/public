/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Plate} from './plate';
import {FIT_FUNCTION_4PL_REGRESSION, IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {AnalysisMappingPanel, AnalysisRequiredFields} from '../plates/views/components/analysis-mapping/analysis-mapping-panel';
import {BaseAnalysisView} from './base-analysis-view';
import {
  AnalysisProperty, PlateProperty, createAnalysisRun,
  getOrCreateProperty, saveAnalysisResult, savePlate
} from '../plates/plates-crud';

export const doseRatioAnalysisProperties: AnalysisProperty[] = [
  {name: 'Curve', type: DG.TYPE.STRING},
  {name: 'Series Name', type: DG.TYPE.STRING},
  {name: 'Min Value', type: DG.TYPE.FLOAT},
  {name: 'Max Value', type: DG.TYPE.FLOAT},
];
export class PlateDoseRatioAnalysis {
  private static REQUIRED_FIELDS: AnalysisRequiredFields[] = [
    {
      name: 'Agonist_Concentration_M',
      required: true,
      description: 'Agonist concentration values in molar units'
    },
    {
      name: 'Antagonist_Concentration_M',
      required: true,
      description: 'Antagonist concentration values in molar units'
    },
    {
      name: 'Percent_Inhibition',
      required: true,
      description: 'Response values as percent inhibition'
    },
    {
      name: 'Agonist_ID',
      required: false,
      description: 'Optional: Identifier for the agonist compound'
    },
    {
      name: 'Antagonist_ID',
      required: false,
      description: 'Optional: Identifier for the antagonist compound'
    }
  ];

  // In PlateDoseRatioAnalysis
  static createAnalysisView(
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (target: string, source: string) => void,
    onUndo: (target: string) => void
  ): HTMLElement {
    const analysisView = new BaseAnalysisView(
      plate,
      {
        analysisName: 'Dose Ratio Analysis',
        requiredFields: this.getRequiredFields(),
        createResultsView: (plate, mappings) => {
          return this.createDoseRatioGrid(plate, mappings);
        }
      },
      currentMappings,
      onMap,
      onUndo
    );

    return analysisView.getRoot();
  }

  static async saveDoseRatioAnalysisResults(
    plate: Plate,
    series: IFitSeries[],
    antagonistConcentrations: number[],
    mappings: Map<string, string>
  ): Promise<void> {
    if (!plate.id) {
      grok.shell.error('Plate must be saved before saving analysis results.');
      return;
    }

    const pi = DG.TaskBarProgressIndicator.create('Saving Dose-Ratio results...');
    try {
      // 1. Groups are the unique antagonist concentrations.
      const groups = antagonistConcentrations.map(String);

      // 2. Create the analysis run.
      const runId = await createAnalysisRun(plate.id, 'Dose-Ratio', groups);

      // Optional: Save parameters like mappings.
      // ...

      // 3. Get or create the property definitions.
      const propertyMap = new Map<string, PlateProperty>();
      propertyMap.set('Curve', await getOrCreateProperty('Curve', DG.TYPE.STRING));
      propertyMap.set('Series Name', await getOrCreateProperty('Series Name', DG.TYPE.STRING));
      propertyMap.set('Min Value', await getOrCreateProperty('Min Value', DG.TYPE.FLOAT));
      propertyMap.set('Max Value', await getOrCreateProperty('Max Value', DG.TYPE.FLOAT));

      // 4. Iterate over each curve (series) and save its properties.
      for (let i = 0; i < series.length; i++) {
        const currentSeries = series[i];
        const antagonistConc = antagonistConcentrations[i];
        const groupCombination = [String(antagonistConc)]; // The group is the antagonist concentration.

        // ... (logic to calculate min/max values remains the same)
        const validPoints = (currentSeries.points || []).filter((p) => !p.outlier && typeof p.y === 'number' && !isNaN(p.y));
        const yValues = validPoints.map((p) => p.y);
        const minValue = yValues.length > 0 ? Math.min(...yValues) : null;
        const maxValue = yValues.length > 0 ? Math.max(...yValues) : null;
        const singleSeriesJson = JSON.stringify({chartOptions: { /* ... */ }, series: [currentSeries]});

        // 5. Save each property for the current group combination.
        const curveProp = propertyMap.get('Curve')!;
        await saveAnalysisResult({runId, propertyId: curveProp.id, propertyType: curveProp.type, value: singleSeriesJson, groupCombination});

        const seriesNameProp = propertyMap.get('Series Name')!;
        await saveAnalysisResult({runId, propertyId: seriesNameProp.id, propertyType: seriesNameProp.type, value: currentSeries.name || `Antagonist: ${antagonistConc}`, groupCombination});

        if (minValue !== null) {
          const minProp = propertyMap.get('Min Value')!;
          await saveAnalysisResult({runId, propertyId: minProp.id, propertyType: minProp.type, value: minValue, groupCombination});
        }
        if (maxValue !== null) {
          const maxProp = propertyMap.get('Max Value')!;
          await saveAnalysisResult({runId, propertyId: maxProp.id, propertyType: maxProp.type, value: maxValue, groupCombination});
        }
      }

      grok.shell.info(`Saved ${series.length} dose-ratio curves for plate ${plate.barcode}.`);
    } catch (e: any) {
      console.error('Save Dose-Ratio Results failed:', e);
      grok.shell.error(`Failed to save dose-ratio results: ${e?.message ?? e}`);
    } finally {
      pi.close();
    }
  }


  static updateAnalysisView(
    container: HTMLElement,
    plate: Plate,
    currentMappings: Map<string, string>,
    onMap: (sourceColumn: string, targetProperty: string) => void,
    onUndo: (targetProperty: string) => void
  ): void {
    const mappingPanel = (container as any).__mappingPanel as AnalysisMappingPanel;
    const resultsContainer = (container as any).__resultsContainer as HTMLElement;

    if (!mappingPanel || !resultsContainer) return;

    mappingPanel.updateConfig({
      sourceColumns: plate.data.columns.names(),
      currentMappings,
      onMap,
      onUndo
    });
    ui.empty(resultsContainer);

    const validationStatus = mappingPanel.getValidationStatus();

    if (validationStatus.isValid) {
      const chartElement = this.createDoseRatioGrid(plate, currentMappings);
      if (chartElement) {
        resultsContainer.appendChild(chartElement);
      } else {
        resultsContainer.appendChild(
          ui.divText('Unable to create dose-ratio curves with current data mapping.', 'warning-message')
        );
      }
    } else {
      const missingFieldsText = validationStatus.missingFields.join(', ');
      resultsContainer.appendChild(
        ui.divText(
          `Please map the following required fields to continue: ${missingFieldsText}`,
          'info-message'
        )
      );
    }
  }


  static createDoseRatioGrid(plate: Plate, mappings?: Map<string, string>): HTMLElement | null {
    const agonistColumn = mappings?.get('Agonist_Concentration_M') || 'Agonist_Concentration_M';
    const antagonistColumn = mappings?.get('Antagonist_Concentration_M') || 'Antagonist_Concentration_M';
    const responseColumn = mappings?.get('Percent_Inhibition') || 'Percent_Inhibition';

    if (!plate.data.columns.contains(agonistColumn) ||
      !plate.data.columns.contains(antagonistColumn) ||
      !plate.data.columns.contains(responseColumn))
      return null;

    const antagonistCol = plate.data.col(antagonistColumn)!;
    const agonistCol = plate.data.col(agonistColumn)!;
    const responseCol = plate.data.col(responseColumn)!;

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

    const saveButton = ui.button('SAVE RESULTS', async () => {
      // The check to ensure the plate is saved first is crucial.
      if (!plate.id) {
        grok.shell.warning('Please save the plate using the CREATE button first, then save the analysis results.');
        // The dialog logic to prompt for saving can be kept if desired.
        return;
      }

      // This now calls the NEWLY REFACTORED save function.
      // It has the same signature but saves to the new tables.
      await PlateDoseRatioAnalysis.saveDoseRatioAnalysisResults(
        plate,
        series,
        antagonistConcentrations,
        mappings || new Map()
      );
    });

    saveButton.style.marginTop = '8px';

    const container = ui.divV([
      grid.root,
      ui.div([saveButton], {style: {display: 'flex', justifyContent: 'flex-end', paddingRight: '4px'}})
    ], 'drc-grid-container');

    container.style.width = '100%';
    container.style.height = '100%';
    container.style.display = 'flex';
    container.style.flexDirection = 'column';

    grid.root.style.width = '100%';
    grid.root.style.flexGrow = '1';

    ui.tools.handleResize(container, (w: number, h: number) => {
      if (w > 20 && h > 20) {
      grid.col(curveCol.name)!.width = w - 20;
      grid.props.rowHeight = h - 20;
      }
    });

    return container;
  }

  static getRequiredFields(): AnalysisRequiredFields[] {
    return [...this.REQUIRED_FIELDS];
  }
}
